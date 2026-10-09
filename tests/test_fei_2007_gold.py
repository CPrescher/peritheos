"""Au replay checks protect measured selection, equations and uncertainty scope."""

import hashlib
import json

import numpy as np
import pytest

from scripts.reproduce_fei_2007_gold import (
    DATA,
    OUTPUT,
    PUBLISHED,
    check_saved,
    ledger_outcome,
    observations,
    pressure,
    read_csv,
    reproduce,
)


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_hirose_original_columns_and_missing_errors_are_preserved():
    rows = read_csv("gold-hirose-2006-table1.csv")
    assert len(rows) == 13
    assert rows[8]["gold_a_angstrom_printed"] == "3.7106(0)"
    assert rows[8]["gold_a_angstrom_uncertainty"] == ""
    assert rows[9]["gold_a_angstrom_printed"] == "3.7220(5)"
    selected = [r for r in rows if r["use_in_fei2007_thermal_replay"] == "true"]
    assert [r["mgo_speziale_pressure_gpa_printed"] for r in selected] == [
        "106.5(12)",
        "117.9(21)",
    ]
    assert all(r["gold_a_angstrom_uncertainty"] == "" for r in selected)
    for r in rows:
        assert float(r["volume_a3"]) == pytest.approx(float(r["gold_a_angstrom"]) ** 3)
        if r["gold_a_angstrom_uncertainty"]:
            assert float(r["volume_uncertainty_a3"]) == pytest.approx(
                3
                * float(r["gold_a_angstrom"]) ** 2
                * float(r["gold_a_angstrom_uncertainty"])
            )
    source = json.loads((DATA / "gold-hirose-2006-source.json").read_text())
    assert (
        source["csv_sha256"]
        == hashlib.sha256(
            (DATA / "gold-hirose-2006-table1.csv").read_bytes()
        ).hexdigest()
    )
    assert source["uncertainties"]["confidence"] is None


def test_hot_targets_use_reported_mgo_pressures():
    rows, data = observations()
    assert len(rows) == 28
    assert data["p"][-2:].tolist() == [106.5, 117.9]
    assert data["t"][-2:].tolist() == [1340.0, 2330.0]
    assert np.isnan(data["sv"][-2:]).all()
    assert data["v"][-2:].tolist() == pytest.approx(np.array([3.7237, 3.7254]) ** 3)


def test_independent_debye_quadrature_and_native_equation(report):
    assert report["independent_equation_states"] == 119
    assert report["independent_equation_max_difference_gpa"] < 1e-9
    assert report["independent_quadrature_max_difference_gpa"] < 1e-9
    assert pressure(PUBLISHED[0], 300) == pytest.approx(0, abs=1e-12)
    # The two laws coincide at q=0, including expanded volumes.
    c = PUBLISHED.copy()
    c[4] = 0
    v = np.array([50.0, 67.85, 70.0])
    assert pressure(v, 2000, c, law="integrated") == pytest.approx(
        pressure(v, 2000, c), abs=1e-10
    )
    assert report["published_error_confidence"] is None
    assert report["published_covariance"] is None


def test_curve_agreement_does_not_identify_law_or_enter_objective(report):
    curves = report["source_curve_checks"]
    assert curves["kind"] == "calculated_source_curves_never_fit_targets"
    assert sum(c["laws"]["printed"]["rows"] for c in curves["curves"]) == 111
    assert all(
        c["laws"][law]["within_one_pressure_direction_stroke"]
        for c in curves["curves"]
        for law in ("printed", "integrated")
    )
    assert max(c["max_law_difference_gpa"] for c in curves["curves"]) < 0.024
    assert report["fits"]["all28_equal"]["refit"]["rows"] == 28


def test_conditional_weight_sensitivity_and_calibration_gap(report):
    fits = report["fits"]
    assert fits["all28_equal"]["parameters"]["q"] == pytest.approx(0.63219425, abs=1e-5)
    assert fits["all28_equal"]["refit"]["rmse_gpa"] == pytest.approx(
        0.71142160, abs=1e-6
    )
    assert fits["all28_available_errors"]["parameters"]["q"] == pytest.approx(
        1.1310102, abs=1e-5
    )
    assert report["available_error_eiv"]["q"] == pytest.approx(1.120668, abs=1e-5)
    assert fits["fei2004_only_equal"]["parameters"]["q"] > 1.1
    assert report["points"][-1]["Au_volume_width_a3"] is None
    assert report["MgO_reduction"]["Hirose2006_difference"][
        "max_abs_residual_gpa"
    ] == pytest.approx(1.9521275, abs=1e-6)
    assert report["author_fit_reproduced"] is False
    assert report["published_coefficients_changed"] is False
    assert report["defaults_changed"] is False


def test_archived_au_report_reproduces(report):
    check_saved(json.loads(OUTPUT.read_text()), report)


def test_staged_parity_keeps_author_fit_and_confidence_unresolved(report):
    from peritheos import get_material_document

    document = get_material_document("gold")
    record = next(
        r for r in document["eos_records"] if r["identifier"] == "gold_fei_2007_vinet_2"
    )
    outcome = ledger_outcome(record)
    assert outcome["status"] == report["classification"]["status"] == "parity"
    assert outcome["parity_basis"] == "within_reported_parameter_errors"
    assert outcome["parameters"] == report["classification"]["parameters"]
    assert outcome["observations"] == 65
    assert outcome["author_fit_reproduced"] is False
    assert all(
        p["within_reported_error"] and p["within_combined_2sigma"] is None
        for p in outcome["parameters"]
    )
    assert [p["refit"] for p in outcome["parameters"]] == pytest.approx(
        [5.995326273, 0.619148632], abs=1e-6
    )
    sv = record["scientific_validation"]
    assert sv["pvt_check"]["status"] == "parity"
    assert sv["pvt_check"]["report"] == "docs/data/fei-2007-gold-reproduction.json"
    assert record["pressure_calibration"]["status"] == "partially_resolved"
    assert record["thermal"]["parameters"]["q"] == 0.6
    ledger = json.loads((OUTPUT.parent / "primary-eos-refits.json").read_text())
    saved = next(
        r for r in ledger["records"] if r["record_identifier"] == record["identifier"]
    )
    assert all(saved[k] == value for k, value in outcome.items())
