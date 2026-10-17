"""Independent equation, normalization and fit comparisons for diagnostic MGD."""

import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq

from scripts.brucite_mgd_comparison import (
    NAMES,
    REPORT,
    check_implementation,
    data,
    debye3,
    pressure,
    reproduce,
    volumes,
)


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_debye_quadrature_against_adaptive_integral():
    ratios = np.array([1e-5, 0.1, 1, 4, 12, 30])
    expected = np.array(
        [
            3 / x**3 * quad(lambda z: z**3 / np.expm1(z), 0, x, epsabs=1e-12)[0]
            for x in ratios
        ]
    )
    np.testing.assert_allclose(debye3(ratios), expected, rtol=2e-13)


@pytest.mark.parametrize("q", [0, -0.2, 1.0, 5.0])
def test_mgd_unit_conversion_reference_and_runtime_parity(q):
    parameters = [40.93, 41.8, 6.6, 0.65, q, 1100]
    checked = check_implementation(parameters)
    assert checked["native_max_difference_gpa"] < 1e-10
    assert checked["python_max_difference_gpa"] < 1e-10
    assert checked["minimum_bulk_modulus_gpa"] > 0
    assert pressure("mgd", 40.93, 300, parameters) == pytest.approx(0, abs=1e-12)


def test_inversion_against_bracketed_roots(report):
    parameters = list(report["fits"]["mgd_volume"]["parameters"].values())
    p = np.array([-0.06, 6, 12, 22.56])
    t = np.array([300, 473, 673, 873])
    expected = [
        brentq(lambda v: pressure("mgd", v, temp, parameters) - target, 30, 45)
        for target, temp in zip(p, t)
    ]
    result = volumes("mgd", p, t, parameters, [41, 37, 35, 32])
    np.testing.assert_allclose(result, expected, atol=2e-11)


def test_matched_rows_objectives_and_conditional_covariance(report):
    import json

    saved = json.loads(REPORT.read_text())
    assert saved["source_file_sha256"] == report["source_file_sha256"]
    assert len(report["row_order"]) == 48
    assert report["row_order"][:37] == [
        {"table": "table4", "source_row": i} for i in range(1, 38)
    ]
    assert report["row_order"][37:] == [
        {"table": "table1", "source_row": i} for i in range(1, 12)
    ]
    dataset = data()
    for key, fit in report["fits"].items():
        assert fit["solver_success"]
        assert not fit["active_bounds"]
        assert fit["jacobian_rank"] == fit["free_parameters"]
        assert fit["observations"] == 48
        residual = np.array(
            fit[
                f"row_{fit['objective']}_residuals_{'gpa' if fit['objective'] == 'pressure' else 'a3'}"
            ]
        )
        assert fit["objective_sse"] == pytest.approx(residual @ residual, abs=1e-11)
        assert fit["objective_sse"] == pytest.approx(
            saved["fits"][key]["objective_sse"], abs=1e-10
        )
        covariance = np.array(fit["conditional_covariance"])
        errors = np.array(list(fit["conditional_standard_errors"].values()))
        np.testing.assert_allclose(covariance.diagonal(), errors**2)
        assert fit["parameter_order"] == NAMES[fit["model"]]
        p = pressure(
            fit["model"],
            dataset["volume"],
            dataset["temperature"],
            list(fit["parameters"].values()),
        )
        np.testing.assert_allclose(
            p - dataset["pressure"], fit["row_pressure_residuals_gpa"], atol=1e-12
        )
        for t, count in [(300, 13), (473, 8), (673, 8), (873, 8)]:
            assert fit["assessment"][f"table4_{t}k"]["observations"] == count


def test_mgd_improves_old_model_but_matched_empirical_fit_is_better(report):
    old = report["current_continuous"]["assessment"]["table4_all"]["volume_rmse_a3"]
    mgd = report["fits"]["mgd_volume"]["assessment"]["table4_all"]["volume_rmse_a3"]
    empirical = report["fits"]["empirical_volume"]["assessment"]["table4_all"][
        "volume_rmse_a3"
    ]
    assert old == pytest.approx(0.1218333762, abs=1e-8)
    assert mgd == pytest.approx(0.0834024441, abs=1e-7)
    assert empirical == pytest.approx(0.0767864456, abs=1e-7)
    assert empirical < mgd < old
    # Joint fitting trades a slightly worse 300 K curve for better hot slices.
    assert (
        report["fits"]["mgd_volume"]["assessment"]["table4_300k"]["volume_rmse_a3"]
        > report["current_continuous"]["assessment"]["table4_300k"]["volume_rmse_a3"]
    )
    pressure_theta = report["fits"]["mgd_pressure"]["parameters"]["theta0_k"]
    volume_theta = report["fits"]["mgd_volume"]["parameters"]["theta0_k"]
    assert pressure_theta > 2 * volume_theta
    assert len(report["theta_sensitivity_profiles"]["volume"]) == 9
