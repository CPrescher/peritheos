"""Source shock data, porous precursor and full experimental path coverage."""

import copy
import hashlib
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    convert_eosmat,
    eosmat_schema,
    get_eos_record,
    get_material,
    get_material_document,
    load_eosmat,
    save_eosmat,
)
from scripts.hu_2019_fe3c import (
    DATASET,
    IDENTIFIER,
    NUMERICAL,
    RHO0,
    ROOT,
    V0,
    catalog_dataset,
    catalog_record,
    observations,
    reproduce,
)
from scripts.standardize_reconstruction_metadata import assessment_for_record


def test_all_source_rows_absolute_error_tokens_and_quantity_origins():
    raw, rows = observations()
    assert len(rows) == 6
    assert [r["shot"] for r in rows] == list(range(1, 7))
    assert raw[0]["shock_velocity_km_s"] == "8.845(0.081)"
    assert rows[0]["shock_velocity_km_s_error_width"] == 0.081
    assert rows[-1]["pressure_gpa"] == 245
    assert rows[-1]["calculated_hugoniot_temperature_k"] == 5530
    assert rows[-1]["calculated_hugoniot_temperature_k_error_width"] == 190
    assert all(r["initial_density_g_cm3_error_width"] == 0.010 for r in rows)
    assert all(r["impact_velocity_km_s_error_width"] is None for r in rows)
    manifest = json.loads(
        (ROOT / "peritheos/data/datasets/cementite-hu-2019-source.json").read_text()
    )
    for resource_key, hash_key in (
        ("source_workbook", "source_workbook_sha256"),
        ("transcription_resource", "transcription_sha256"),
    ):
        path = ROOT / "peritheos/data" / manifest[resource_key]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == manifest[hash_key]
    assert manifest["source_cell_values"]["E8"] == raw[-1]["shock_velocity_km_s"]
    d = get_material_document("cementite_nonmagnetic", format_version=4)
    ds = next(x for x in d["datasets"] if x["identifier"] == DATASET)
    assert ds == catalog_dataset()
    assert ds["license"] == "CC-BY-4.0"
    assert "NOT measured" in ds["quantity_origins"]["calculated_hugoniot_temperature_k"]
    typed = get_material("cementite_nonmagnetic").get_dataset(DATASET)
    np.testing.assert_array_equal(typed["solid_subset"], [1, 1, 1, 0, 0, 0])
    np.testing.assert_array_equal(
        typed["initial_density_g_cm3"], [7.546, 7.524, 7.535, 7.504, 7.509, 7.478]
    )


def test_published_coefficients_preserved_and_conditional_fit_is_separate(
    assert_audit_close,
):
    result = reproduce()
    assert_audit_close(
        result, json.loads((ROOT / NUMERICAL).read_text()), rel=1e-12, abs=1e-12
    )
    assert result["source_equation_check"]["all_within_reported_Us_width"]
    assert result["source_equation_check"]["rms_km_s"] == pytest.approx(0.03468433578)
    refit = result["independent_diagnostic_refit"]
    np.testing.assert_allclose(refit["parameters_c0_s"], [4.62819684355, 1.61907248479])
    assert all(refit["within_published_coefficient_widths"])
    assert not refit["author_workflow_reproduced"]
    r = catalog_record()
    assert r["parameters"]["c0"] == {"value": 4.627, "error": 0.073}
    assert r["parameters"]["s"] == {"value": 1.614, "error": 0.028}
    assert r["parameter_error_confidence"] is None
    assert r["record_kind"] == "published"
    for s in result["rankine_hugoniot_replay"]["rows"]:
        row = observations()[1][s["shot"] - 1]
        assert abs(s["pressure_difference_gpa"]) < row["pressure_gpa_error_width"]
        assert abs(s["density_difference_g_cm3"]) < row["density_g_cm3_error_width"]


def test_saved_record_assessment_matches_ledger_and_existing_default(
    assert_record_replay_close,
):
    d = get_material_document("cementite_nonmagnetic", format_version=4)
    r = next(x for x in d["eos_records"] if x["identifier"] == IDENTIFIER)
    assert_record_replay_close(r, catalog_record())
    entry = next(
        x
        for x in json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())[
            "records"
        ]
        if x["record_identifier"] == IDENTIFIER
    )
    assert entry["observations"] == 6
    assert r["evidence"]["reproduction"]["assessment"] == assessment_for_record(
        r, entry
    )
    material = get_material("cementite_nonmagnetic")
    assert (
        material.default_equilibrium_record().identifier
        == "cementite_nonmagnetic_takahashi_2019_kbr_bm3_mgd"
    )
    assert material.default_hugoniot_record().identifier == IDENTIFIER
    assert r["source_fit_domain"]["pressure_gpa"] == [85.6, 245]
    assert r["source_fit_domain"]["installed_shots"] == [1, 2, 3, 4, 5, 6]
    assert r["validity"]["pressure_gpa"] == [85.6, 245]
    assert "Possible melting" in r["evidence"]["scientific_validation"]["warning"]


def test_porous_initial_state_and_full_experimental_path_execution():
    r = get_eos_record(IDENTIFIER)
    assert r.equation_kind == "hugoniot"
    assert r.reference_volume == pytest.approx(V0)
    assert V0 == pytest.approx(158.481407231785)
    assert r.eos.rho0 == RHO0
    state = r.state_from_particle_velocity(2.164)
    assert state.shock_velocity == pytest.approx(8.119696)
    assert state.pressure == pytest.approx(132.2219416336)
    assert state.density == pytest.approx(10.2592060441)
    for row in observations()[1]:
        assert np.isfinite(
            r.state_from_particle_velocity(row["particle_velocity_km_s"]).pressure
        )
    for p in [85.6, 100, 132.5, 171.6, 200, 227.5, 245]:
        v = r.volume(p, check_validity=True)
        assert r.pressure(v, check_validity=True) == pytest.approx(p, abs=1e-8)
    with pytest.raises(ValueError):
        r.state_from_particle_velocity(3.4)
    with pytest.raises(ValueError):
        r.volume(250, check_validity=True)
    with pytest.raises(ValueError):
        r.volume(100, temperature=1190)


@pytest.mark.parametrize("version", [3, 4])
def test_native_and_compatibility_roundtrip_preserve_porosity_rows_and_scope(
    tmp_path, version
):
    d = get_material_document("cementite_nonmagnetic", format_version=4)
    before = copy.deepcopy(d)
    out = convert_eosmat(d, version=version)
    Draft202012Validator(eosmat_schema(version)).validate(out)
    p = tmp_path / "fe3c.eosmat"
    save_eosmat(p, out)
    restored = convert_eosmat(load_eosmat(p), version=4)
    assert restored["eos_records"] == d["eos_records"]
    assert restored["datasets"] == d["datasets"]
    chosen = Material.from_eosmat(out, record_identifiers=[IDENTIFIER])
    assert len(chosen.get_dataset(DATASET)) == 6
    assert chosen.eos_records[0].initial_state.density_g_cm3 == RHO0
    assert d == before


def test_generic_us_up_report_replays_all_six_rows_and_both_error_coordinates():
    from scripts.reconstruction_plot_data import curve_comparison

    native = get_material_document("cementite_nonmagnetic", format_version=4)
    r = next(x for x in native["eos_records"] if x["identifier"] == IDENTIFIER)
    comparison, reason = curve_comparison(native, r)
    assert not reason
    panel = comparison.panels[0]
    assert len(panel.x) == 6
    assert panel.rms("source") == pytest.approx(
        reproduce()["source_equation_check"]["rms_km_s"]
    )
    assert panel.rms("refit") == pytest.approx(
        reproduce()["independent_diagnostic_refit"]["rms_km_s"]
    )
    np.testing.assert_array_equal(
        panel.sigma, [0.081, 0.034, 0.113, 0.036, 0.122, 0.055]
    )
    np.testing.assert_array_equal(
        panel.x_sigma, [0.016, 0.009, 0.016, 0.029, 0.029, 0.030]
    )
