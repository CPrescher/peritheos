"""Keep source transcription, refit repeatability and graphical checks separate."""

import copy
import hashlib
import json

import numpy as np
import pytest

from peritheos import (
    convert_eosmat,
    get_material,
    load_eosmat,
    validate_eosmat_document,
)
from scripts import reproduce_bykova_2018_coesite_iv as audit
from scripts.eosmat_io import read_material
from scripts.generate_reconstruction_reports import report_page
from scripts.standardize_reconstruction_metadata import assessment_for_record


def document():
    return load_eosmat(audit.ROOT / "peritheos/data/materials/coesite_iv.eosmat")


def test_original_coefficients_and_source_table_are_preserved():
    native = document()
    record = native["eos_records"][0]
    assert record["parameters"] == {
        "V0": {"value": 438.6702, "error": None},
        "K0": {"value": 168.52, "error": None, "fixed": True},
        "K0_prime": {"value": 3.34, "error": None, "fixed": True},
    }
    assert hashlib.sha256(audit.TABLE.read_bytes()).hexdigest() == (
        "f0f2c318a6f54f757f818cbc5e5d7dbc5b1b02c3b53fb0a55598eb79064869e4"
    )
    source = np.genfromtxt(audit.TABLE, delimiter=",", names=True)
    np.testing.assert_array_equal(source["pressure_gpa"], [39, 41, 48])
    np.testing.assert_array_equal(
        source["volume_a3_conventional_cell"], [368.5011, 366.1221, 357.8274]
    )
    assert record["record_kind"] == "refit"
    assert record["fit_datasets"] == [audit.TABLE_ID]
    assert native["datasets"][1]["used_by_eos_records"] == []
    assert record["parameters"]["V0"]["error"] is None
    validate_eosmat_document(native)
    recovered = convert_eosmat(convert_eosmat(native, version=3), version=4)
    assert recovered["eos_records"] == native["eos_records"]
    assert recovered["datasets"] == native["datasets"]
    material = get_material("coesite_iv")
    figure = material.get_dataset(audit.FIGURE_ID)
    assert figure is not None


def test_repeatability_scope_and_covariance_do_not_claim_published_v0():
    native = document()
    legacy = read_material(audit.ROOT / "peritheos/data/materials/coesite_iv.eosmat")
    row = {
        "record_identifier": audit.RECORD_ID,
        "material_identifier": "coesite_iv",
        **audit.ledger_outcome(legacy["eos_records"][0]),
    }
    assessment = native["eos_records"][0]["evidence"]["reproduction"]["assessment"]
    assert assessment == assessment_for_record(native["eos_records"][0], row)
    assert assessment["scope"] == "independent_refit"
    assert assessment["basis_code"] == "registered_peritheos_refit_reproduced"
    assert assessment["status"] == "parity"
    assert "not a published coefficient" in row["comparison_scope"]
    assert row["parameters"][0]["within_combined_2sigma"] is None
    assert row["absolute_sigma"] is False
    assert "reduced_chi_square" not in row
    assert "conditional" in row["covariance_convention"]
    page = report_page(native, root=audit.ROOT)
    assert "stored **Peritheos refit**" in page
    assert "whole GPa" in page
    assert "coesite-iv-bykova-2018.md" in page
    changed = copy.deepcopy(legacy["eos_records"][0])
    changed["eos"]["parameters"]["K0"] += 1
    with pytest.raises(ValueError, match="coefficients differ"):
        audit.ledger_outcome(changed)


def test_table_rounding_diagnostic_is_not_a_statistical_confidence_interval():
    result = audit.reproduce()
    fit = result["fit"]
    assert fit["V0_a3"] == pytest.approx(438.6702185, abs=1e-6)
    assert fit["rmse_gpa"] == pytest.approx(0.0808114, abs=1e-7)
    assert fit["formal_V0_standard_error_a3"] == pytest.approx(0.0842556, abs=1e-6)
    rounding = result["rounding_diagnostic"]
    assert rounding["confidence_interval"] is False
    assert rounding["intersection_V0_interval_a3"] == pytest.approx(
        [438.0894721, 439.2672041], abs=1e-6
    )
    table = np.genfromtxt(audit.TABLE, names=True, delimiter=",")
    for v0 in rounding["intersection_V0_interval_a3"]:
        residuals = (
            audit.pressure(table["volume_a3_conventional_cell"], v0)
            - table["pressure_gpa"]
        )
        assert max(abs(residuals)) <= 0.500000001


def test_figure_readings_are_independent_checks_and_bounds_are_not_errors():
    assert audit.FIGURE.read_text() == audit.digitized_csv()
    checks = audit.reproduce()["figure_check"]
    assert checks["coefficients_fitted_to_figure"] is False
    assert checks["all_markers"]["points"] == 8
    extra = checks["additional_markers_only"]
    assert extra["points"] == 5
    assert extra["all_within_reading_bounds"] is True
    assert extra["pressure_rmse_gpa"] == pytest.approx(0.2259419, abs=1e-7)
    assert extra["max_absolute_pressure_residual_gpa"] == pytest.approx(
        0.3544235, abs=1e-7
    )
    assert checks["two_pixel_sensitivity"]["markers_outside_bounds"] == [3, 4]
    assert [r["marker"] for r in checks["points"] if r["table10_counterpart"]] == [
        5,
        7,
        8,
    ]
    figure = document()["datasets"][1]
    assert (
        hashlib.sha256(audit.FIGURE.read_bytes()).hexdigest()
        == figure["resource"]["sha256"]
    )
    assert all(
        c["role"] == "bound"
        for c in figure["columns"]
        if "digitization_bound" in c["name"]
    )
    # Public native evaluation must agree with the separate hand expression at
    # every independently read graphical state, without optimizing either one.
    model = get_material("coesite_iv").get_eos_record(audit.RECORD_ID).eos
    v = np.array([r["volume_a3_conventional_cell"] for r in checks["points"]])
    np.testing.assert_allclose(model.pressure(v), audit.pressure(v), atol=1e-10)


def test_bad_figure_readings_fail_validation_without_changing_the_fit(monkeypatch):
    before = audit.reproduce()
    rows = audit.digitized_rows()
    for row in rows:
        row["volume_a3_per_formula_unit"] += 0.2
        row["volume_a3_conventional_cell"] += 3.2
    monkeypatch.setattr(audit, "digitized_rows", lambda: rows)
    after = audit.reproduce()
    assert before["fit"] == after["fit"]
    assert (
        after["figure_check"]["additional_markers_only"]["all_within_reading_bounds"]
        is False
    )


def test_saved_numerical_audit_replays():
    from scripts.validate_primary_eos_refits import check_ledger_values

    check_ledger_values(json.loads(audit.RESULT.read_text()), audit.reproduce())
