"""Matched pressure-recalculation refits remain separate from the source EOS."""

import csv
import json

import numpy as np
import pytest
from scipy.optimize import least_squares

from scripts.audit_sata_2010_fe3c import VR, bm3_fit, selected_rows
from scripts.refit_sata_2010_recalculated_pressure import INPUTS, REPORT, reproduce


@pytest.fixture(scope="module")
def report():
    return reproduce()


def independent_bm3(volume, parameters):
    k, kp = parameters
    ratio = VR / volume
    return (
        1.5
        * k
        * (ratio ** (7 / 3) - ratio ** (5 / 3))
        * (1 + 0.75 * (kp - 4) * (ratio ** (2 / 3) - 1))
    )


def test_only_sata_pressure_coordinates_change_and_li_points_stay_original(report):
    original = selected_rows()
    assert len(report["inputs"]) == len(original) == 12
    assert report["kind"] == "diagnostic_refit"
    assert report["catalog_integration"] is False
    assert report["temperature_k"] == 300
    for row, source in zip(report["inputs"], original):
        for name in (
            "source",
            "source_row",
            "run",
            "measurement_id",
            "volume_a3",
            "error_volume_a3",
            "temperature_k",
        ):
            assert row[name] == source[name]
        assert row["published_pressure_gpa"] == source["pressure_gpa"]
        assert row["published_pressure_error_width_gpa"] == source["error_pressure_gpa"]
        if row["source"] == 1:
            assert row["fit_pressure_gpa"] == source["pressure_gpa"]
            assert row["pressure_shift_gpa"] == 0
            assert row["pressure_basis"] == "unchanged_Li_original_ruby_neon"
            assert row["run"] is None and row["measurement_id"] is None
        else:
            assert -0.353 < row["pressure_shift_gpa"] < 0
    with INPUTS.open(newline="") as stream:
        archived = list(csv.DictReader(stream))
    assert [float(r["fit_pressure_gpa"]) for r in archived[9:]] == [27.2, 30.3, 30.5]
    assert all(r["run"] == r["measurement_id"] == "" for r in archived[9:])


def test_linear_refits_match_independent_nonlinear_objective_and_prior_diagnostic(
    report,
):
    v = np.array([r["volume_a3"] for r in report["inputs"]])
    for key, coordinate in [
        ("original_pressures", "published_pressure_gpa"),
        ("recalculated_Sata_unchanged_Li", "fit_pressure_gpa"),
    ]:
        observed = np.array([r[coordinate] for r in report["inputs"]])
        reference = least_squares(
            lambda values: independent_bm3(v, values) - observed,
            [280, 3.8],
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-10,
        )
        assert reference.success
        fit = report["fixed_V0_BM3"][key]
        np.testing.assert_allclose(fit["parameters"], reference.x, atol=2e-5)
        assert fit["degrees_of_freedom"] == 10
        assert np.all(np.linalg.eigvalsh(fit["conditional_covariance"]) > 0)
        assert fit["rmse_gpa"] == pytest.approx(
            np.sqrt(np.mean(reference.fun**2)), abs=1e-10
        )
    np.testing.assert_allclose(
        report["fixed_V0_BM3"]["original_pressures"]["parameters"],
        bm3_fit()["parameters"],
        atol=1e-10,
    )
    new = report["fixed_V0_BM3"]["recalculated_Sata_unchanged_Li"]
    np.testing.assert_allclose(
        new["parameters"], [290.524467979, 3.752551645], atol=1e-8
    )
    assert new["rmse_gpa"] == pytest.approx(1.07854111345)
    assert new["max_abs_gpa"] > 2.38


def test_both_pressure_coordinates_and_reference_form_are_compared_without_promotion(
    report,
):
    cross = report["cross_evaluation"]
    assert cross["published_BM3"]["on_original_pressures"]["rmse_gpa"] == pytest.approx(
        1.10673314549
    )
    assert cross["published_BM3"]["on_recalculated_Sata_unchanged_Li"][
        "rmse_gpa"
    ] == pytest.approx(1.08076738921)
    check = report["reference_mBM3_check"]
    np.testing.assert_allclose(
        check["recalculated_Sata_unchanged_Li"]["parameters_Pr_Kr_Kprime"],
        [0.033908277, 290.272067286, 3.756051272],
        atol=1e-5,
    )
    assert (
        check["original_pressures"]["dof"]
        == check["recalculated_Sata_unchanged_Li"]["dof"]
        == 9
    )
    assert len(report["residuals"]) == 12
    archived = json.loads(REPORT.read_text())
    assert report == archived
