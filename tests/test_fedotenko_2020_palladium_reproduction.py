"""Printed weights must survive an unreported ambient pressure uncertainty."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_fedotenko_2020_palladium import (
    DATA,
    RECORDS,
    fit,
    ledger_outcome,
    observations,
)
from scripts.validate_primary_eos_refits import _fit_record, check_ledger_values


def record(identifier):
    return next(
        r
        for r in get_material_document("palladium")["eos_records"]
        if r["identifier"] == identifier
    )


@pytest.mark.parametrize("identifier", sorted(RECORDS))
def test_fedotenko_retains_reported_widths_before_generic_fallback(
    identifier, monkeypatch
):
    def forbidden(*args):
        pytest.fail(
            "Fedotenko must not discard its pressure widths in the generic reader"
        )

    monkeypatch.setattr("scripts.validate_primary_eos_refits._series", forbidden)
    doc = get_material_document("palladium")
    result = _fit_record(doc, record(identifier), doc["datasets"][0])
    assert result["observations"] == 15
    assert result["weighting"]["positive_pressure_rows"] == 14
    assert result["weighting"]["pressure_sigma_gpa"] == 0.2
    assert result["weighting"]["missing_ambient_pressure_sigma"] == "not imputed"
    assert result["coefficient_comparison_status"] == "parity"
    assert result["author_workflow_reproduced"] is False


def test_absolute_coordinate_widths_control_covariance_and_chi_square():
    r = record("palladium_fedotenko_2020_bm2_2")
    rows = observations()
    original = fit(r, rows)
    rows["pressure_error_gpa"][1:] *= 5
    rows["volume_error_a3_conventional_cell"][1:] *= 5
    scaled = fit(r, rows)
    assert scaled.parameters["K0"] == pytest.approx(original.parameters["K0"], abs=1e-5)
    assert scaled.standard_errors["K0"] == pytest.approx(
        5 * original.standard_errors["K0"], rel=1e-5
    )
    assert scaled.reduced_chi_square == pytest.approx(original.reduced_chi_square / 25)
    assert original.degrees_of_freedom == 13
    assert original.reduced_chi_square == pytest.approx(25.0054834)


def test_poor_bm2_agreement_is_separate_from_coefficient_parity():
    outcome = ledger_outcome(record("palladium_fedotenko_2020_bm2_2"))
    rows = observations()
    # Independent literal BM2 expression, evaluated at measured coordinates.
    ratio = 58.868 / rows["volume_a3_conventional_cell"]
    residual = 1.5 * 203 * (ratio ** (7 / 3) - ratio ** (5 / 3)) - rows["pressure_gpa"]
    assert outcome["published_rmse_gpa"] == pytest.approx(
        np.sqrt(np.mean(residual**2)), abs=1e-10
    )
    assert outcome["published_max_abs_pressure_residual_gpa"] == pytest.approx(
        1.4385589886
    )
    assert outcome["rmse_gpa"] == pytest.approx(0.9337035754)
    assert outcome["status"] == "parity"
    assert outcome["reproduction_status"] == "parity"
    assert outcome["coefficient_comparison_status"] == "parity"
    assert outcome["data_agreement_status"] == "poor"
    assert outcome["usage_recommendation"] == "not_recommended_for_quantitative_use"


def test_saved_fedotenko_assessments_match_the_corrected_audit():
    root = Path(__file__).resolve().parents[1]
    rows = json.loads((root / "docs/data/primary-eos-refits.json").read_text())[
        "records"
    ]
    saved = {r["record_identifier"]: r for r in rows}
    for identifier in RECORDS:
        r = record(identifier)
        current = ledger_outcome(r)
        # validate_all appends the reviewed recommendation after fitting.
        reason = r["scientific_validation"].get("recommendation_reason")
        if reason:
            current["qualification"] += " " + reason
        check_ledger_values({k: saved[identifier][k] for k in current}, current)


def test_published_coefficients_errors_and_observations_are_preserved():
    assert (
        hashlib.sha256(DATA.read_bytes()).hexdigest()
        == "1711e9f0bae4905cb18fbb8abce1ee79f9c1127cc718962a95350fefe7e065b1"
    )
    bm2 = record("palladium_fedotenko_2020_bm2_2")
    assert bm2["eos"]["parameters"] == {"V0": 58.868, "K0": 203.0}
    assert bm2["parameter_errors"] == {"V0": 0.002, "K0": 3.0}
    bm3 = record("palladium_fedotenko_2020_bm3_1")
    assert bm3["eos"]["parameters"] == {"V0": 58.868, "K0": 157.0, "K0_prime": 9.9}
    assert bm3["parameter_errors"] == {"V0": 0.002, "K0": 3.0, "K0_prime": 0.4}
    assert np.isnan(observations()["pressure_error_gpa"][0])


def test_bm2_remains_available_only_by_explicit_diagnostic_selection():
    rid = "palladium_fedotenko_2020_bm2_2"
    doc = get_material_document("palladium")
    assert rid not in {r.identifier for r in Material.from_eosmat(doc).eos_records}
    material = Material.from_eosmat(doc, record_identifiers=[rid])
    raw = next(r for r in material.to_eosmat()["eos_records"] if r["identifier"] == rid)
    assert raw["scientific_validation"]["data_agreement_status"] == "poor"
    assert raw["scientific_validation"]["status"] == "primary_source_validated"
    assert "Do not use" in raw["scientific_validation"]["warning"]
    assert "poor fit; diagnostic only" in raw["label"]
    assert raw["reproduction"]["assessment"]["status"] == "parity"
    assert raw["eos"] == record(rid)["eos"]
    model = get_eos_record(rid)
    assert not model.is_default
    assert model.pressure(52.411, 300) == pytest.approx(29.7614410114)
    assert model.volume(model.pressure(52.411, 300), 300) == pytest.approx(52.411)
    report = (
        Path(__file__).resolve().parents[1] / "docs/reconstruction-reports/palladium.md"
    )
    assert (
        "**Outcome: Coefficient parity achieved; poor data agreement.**"
        in report.read_text()
    )
    assert "**Warning:** Do not use" in report.read_text()


def test_source_reaudit_preserves_parity_and_negative_use_recommendation():
    from scripts.apply_primary_source_audit import audit_record

    r = record("palladium_fedotenko_2020_bm2_2")
    assert audit_record(r, "palladium.eosmat") == r
