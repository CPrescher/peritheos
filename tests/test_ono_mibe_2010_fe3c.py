"""Original experiment, finite scope and preservation of source distinctions."""

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
    resolve_dataset_observations,
    save_eosmat,
)
from scripts.ono_mibe_2010_fe3c import (
    DATASET,
    IDENTIFIER,
    NUMERICAL,
    ROOT,
    catalog_dataset,
    catalog_record,
    observations,
    pressure,
    reproduce,
)
from scripts.standardize_reconstruction_metadata import assessment_for_record


def test_original_parameters_adopted_anchor_and_strict_scope(
    assert_record_replay_close,
):
    d = get_material_document("cementite", format_version=4)
    r = next(r for r in d["eos_records"] if r["identifier"] == IDENTIFIER)
    expected = catalog_record()
    assert_record_replay_close(r, expected)
    assert r["parameters"] == {
        "V0": {"value": 155.4, "error": 0.2},
        "K0": {"value": 167, "error": None},
        "K0_prime": {"value": 6.7, "error": None},
    }
    assert r["temperature_ref"] == 300
    assert r["experimental_pressure_range_gpa"] == [0, 34.3]
    assert r["validity"] == catalog_record()["validity"]
    assert r["record_kind"] == "published"
    a = r["evidence"]["reproduction"]["assessment"]
    assert a["scope"] == "source_equation"
    assert a["status"] == "source_reconstruction"
    entry = next(
        row
        for row in json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())[
            "records"
        ]
        if row["record_identifier"] == IDENTIFIER
    )
    assert a == assessment_for_record(r, entry)
    assert r["evidence"]["reproduction"]["author_workflow_reproduced"] is False
    assert (
        get_material("cementite").default_equilibrium_record().identifier
        == "cementite_scott_2001_bm3_1"
    )


def test_all_fourteen_literal_rows_and_error_tokens_are_retained():
    raw, rows = observations()
    assert len(rows) == 14
    assert raw[4]["a_a"] == "4.901(10)"
    assert rows[4]["a_a_error_width"] == 0.010
    assert rows[0]["volume_a3"] == 155.4
    assert rows[-1]["volume_a3"] == 125.1
    assert rows[-1]["pressure_gpa"] == 66.5
    assert all(r["pressure_gpa_error_width"] is None for r in rows)
    ds = next(
        d
        for d in get_material_document("cementite", format_version=4)["datasets"]
        if d["identifier"] == DATASET
    )
    assert ds == catalog_dataset()
    typed = get_material("cementite").get_dataset(DATASET)
    np.testing.assert_array_equal(typed["source_row"], np.arange(1, 15))
    assert np.isnan(typed["pressure_gpa_error_width"]).all()
    obs = resolve_dataset_observations(
        get_material_document("cementite"), IDENTIFIER, DATASET
    )
    assert len(obs.row_indices) == 14
    np.testing.assert_array_equal(
        obs.selected_for_fit, [False] + [True] * 5 + [False] * 8
    )


def test_source_curve_check_and_separate_conditional_refit(assert_audit_close):
    r = reproduce()
    archived = json.loads((ROOT / NUMERICAL).read_text())
    assert_audit_close(
        r, archived, rel=1e-12, abs=1e-10, ignore_keys={"independent_diagnostic_refit"}
    )
    # Subtracting the published K-prime magnifies relative roundoff near zero;
    # bound the diagnostic in coefficient units, keeping its fixed anchor exact.
    assert_audit_close(
        r["independent_diagnostic_refit"],
        archived["independent_diagnostic_refit"],
        abs=2e-7,
        float_paths=(
            "report.parameters_K0_K0_prime*",
            "report.difference_from_published*",
            "report.rms_pressure_gpa",
        ),
    )
    c = r["source_equation_check"]
    assert c["passed"]
    assert c["pressure_gpa"] == [10, 20, 30]
    np.testing.assert_allclose(c["predicted_volume_a3"], c["read_volume_a3"], atol=0.41)
    assert r["source_observation_comparison"]["rms_pressure_gpa"] == pytest.approx(
        0.41015257
    )
    f = r["independent_diagnostic_refit"]
    np.testing.assert_allclose(
        f["parameters_K0_K0_prime"], [165.92506, 6.70657], atol=1e-5
    )
    assert not f["coefficient_parity_claimed"]
    assert not f["author_workflow_reproduced"]


def test_original_nacl_calibration_is_not_later_recalibration():
    r = catalog_record()
    assert (
        r["pressure_calibration"]["recalculation"]["status"]
        == "missing_calibrant_observations"
    )
    assert "NaCl" in r["pressure_calibration"]["notes"]
    assert all(
        "reference_eos_record" not in x for x in r["pressure_calibration"]["methods"]
    )
    assert catalog_dataset()["provenance"]["calibrant_checkpoint"]["NaCl_a_a"] == 2.872
    assert "0 K" in r["notes"]


def test_public_execution_independent_expression_inversion_and_validity():
    m = get_eos_record(IDENTIFIER)
    v = np.linspace(136, 155.4, 20)
    np.testing.assert_allclose(m.pressure(v), pressure(v), atol=1e-11)
    pp = np.linspace(0, 34.3, 20)
    np.testing.assert_allclose(m.pressure(m.volume(pp)), pp, atol=1e-8)
    assert m.eos.bulk_modulus(155.4) == pytest.approx(167)
    assert m.volume(30) == pytest.approx(137.7551, abs=0.41)
    assert np.isfinite(m.volume(34.3, 300, check_validity=True))
    with pytest.raises(ValueError):
        m.volume(40, 300, check_validity=True)
    with pytest.raises(ValueError):
        m.volume(20, 301, check_validity=True)


@pytest.mark.parametrize("version", [3, 4])
def test_native_and_compatibility_roundtrip_preserve_rows_and_record(tmp_path, version):
    d = get_material_document("cementite", format_version=4)
    out = convert_eosmat(d, version=version)
    Draft202012Validator(eosmat_schema(version)).validate(out)
    p = tmp_path / "cementite.eosmat"
    save_eosmat(p, out)
    restored = convert_eosmat(load_eosmat(p), version=4)
    assert restored["datasets"] == d["datasets"]
    assert restored["eos_records"] == d["eos_records"]
    assert (
        Material.from_eosmat(out, record_identifiers=[IDENTIFIER])
        .eos_records[0]
        .identifier
        == IDENTIFIER
    )
