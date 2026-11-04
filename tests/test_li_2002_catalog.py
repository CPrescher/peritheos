"""Published Li coefficients, source observations and qualified public metadata."""

import copy
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
from scripts.li_2002_catalog import DATASET, IDENTIFIER, catalog_dataset, ledger_outcome
from scripts.reproduce_li_2002_fe3c import ROOT, arrays, pressure, rows
from scripts.standardize_reconstruction_metadata import assessment_for_record


def record(document):
    return next(r for r in document["eos_records"] if r["identifier"] == IDENTIFIER)


def test_published_coefficients_fixed_measured_error_and_source_scope(
    assert_assessment_close,
):
    d = get_material_document("cementite", format_version=4)
    r = record(d)
    assert r["parameters"] == {
        "V0": {"value": 155.28, "error": 0.11, "fixed": True},
        "K0": {"value": 174, "error": 6},
        "K0_prime": {"value": 4.8, "error": 0.8},
    }
    assert r["record_kind"] == "published"
    assert r["parameter_error_confidence"] is None
    assert r["temperature_ref"] == 300
    assert r["experimental_pressure_range_gpa"] == [0.0001, 30.5]
    assert (
        r["evidence"]["scientific_validation"]["status"] == "primary_source_validated"
    )
    assert r["evidence"]["reproduction"]["author_workflow_reproduced"] is False
    a = r["evidence"]["reproduction"]["assessment"]
    assert a["status"] == "parity" and a["scope"] == "parameter_agreement"
    assert "confidence" in " ".join(a["qualifications"])
    outcome = {
        **ledger_outcome(r),
        "material_identifier": "cementite",
        "record_identifier": IDENTIFIER,
    }
    assert_assessment_close(a, assessment_for_record(r, outcome))
    assert [p["parameter"] for p in a["parameters"]] == ["K0", "K0_prime"]
    assert all(p["within_combined_2sigma"] is None for p in a["parameters"])
    assert (
        get_material("cementite").default_equilibrium_record().identifier
        == "cementite_scott_2001_bm3_1"
    )


def test_installed_table_and_tokens_are_complete_not_a_calculated_curve():
    d = get_material_document("cementite", format_version=4)
    ds = next(x for x in d["datasets"] if x["identifier"] == DATASET)
    assert ds == catalog_dataset()
    assert ds["source_tokens"] == rows()[0]
    typed = get_material("cementite").get_dataset(DATASET)
    assert len(typed) == 14
    np.testing.assert_array_equal(typed["source_row"], np.arange(1, 15))
    assert typed["volume_a3"][0] == 155.40
    assert typed["volume_a3"][1] == 155.28
    assert typed["pressure_gpa_error_width"][1] == 1e-6
    assert np.isnan(typed["neon_a_a"][:4]).all()
    assert np.isfinite(typed["neon_a_a"][4:]).all()
    assert all(
        c.role == "uncertainty" for c in typed.columns if c.name.endswith("error_width")
    )


def test_observation_view_excludes_reference_card_and_marks_fixed_anchor():
    d = get_material_document("cementite")
    obs = resolve_dataset_observations(d, IDENTIFIER, DATASET)
    assert len(obs.row_indices) == 13
    np.testing.assert_array_equal(obs.row_indices, np.arange(1, 14))
    np.testing.assert_array_equal(obs.selected_for_fit, [False] + [True] * 12)
    _, (p, v, _, _) = arrays()
    np.testing.assert_array_equal(obs.pressure_gpa[obs.selected_for_fit], p)
    np.testing.assert_array_equal(obs.volume_a3[obs.selected_for_fit], v)


def test_original_calibration_is_not_replaced_by_a_later_scale():
    r = record(get_material_document("cementite", format_version=4))
    c = r["pressure_calibration"]
    assert c["recalculation"]["status"] == "missing_calibrant_observations"
    assert c["methods"][0]["parameters"] == {
        "A_gpa": 1904,
        "B": 7.715,
        "reference_wavelength_nm": 694.3,
    }
    assert c["methods"][1]["reference"]["year"] == 1989
    assert all("reference_eos_record" not in x for x in c["methods"])
    assert "average" in c["notes"]


def test_public_execution_native_inversion_and_source_benchmark():
    model = get_eos_record(IDENTIFIER)
    _, (_, v, _, _) = arrays()
    np.testing.assert_allclose(model.pressure(v), pressure(v), atol=2e-12)
    pp = np.linspace(0, 30.5, 30)
    np.testing.assert_allclose(model.pressure(model.volume(pp)), pp, atol=1e-8)
    assert model.pressure(155.28) == pytest.approx(0, abs=1e-12)
    assert model.eos.bulk_modulus(155.28) == pytest.approx(174)
    assert model.volume(30) == pytest.approx(136.73043478, abs=0.12)


@pytest.mark.parametrize("version", [3, 4])
def test_schema_and_roundtrip_preserve_every_source_token_and_metadata(
    tmp_path, version
):
    d = get_material_document("cementite", format_version=4)
    before = copy.deepcopy(d)
    out = convert_eosmat(d, version=version)
    Draft202012Validator(eosmat_schema(version=version)).validate(out)
    path = tmp_path / f"cementite-v{version}.eosmat"
    save_eosmat(path, out)
    restored = convert_eosmat(load_eosmat(path), version=4)
    assert restored["datasets"] == before["datasets"]
    assert restored["eos_records"] == before["eos_records"]
    selected = Material.from_eosmat(out, record_identifiers=[IDENTIFIER])
    assert len(selected.eos_records) == 1
    assert len(selected.get_dataset(DATASET)) == 14
    assert d == before


def test_ledger_and_both_curves_replay_the_same_twelve_observations(
    assert_parameter_replay_close,
):
    d = get_material_document("cementite")
    r = record(d)
    outcome = ledger_outcome(r)
    stored = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    saved = next(x for x in stored["records"] if x["record_identifier"] == IDENTIFIER)
    assert_parameter_replay_close(saved["parameters"], outcome["parameters"])
    assert saved["author_workflow_reproduced"] is False
    _, (p, v, _, _) = arrays()
    assert len(p) == 12
    source_pressure = get_eos_record(IDENTIFIER).pressure(v)
    np.testing.assert_allclose(source_pressure, pressure(v), atol=2e-12)
    assert np.sqrt(np.mean((source_pressure - p) ** 2)) == pytest.approx(
        0.2633934456555
    )
    refit_parameters = [x["refit"] for x in saved["parameters"]]
    refit_pressure = pressure(v, refit_parameters)
    assert np.sqrt(np.mean((refit_pressure - p) ** 2)) == pytest.approx(
        outcome["rmse_gpa"]
    )
