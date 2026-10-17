"""Verify quadratic temperature laws under the unchanged two-coordinate fit."""

import json

import numpy as np
import pytest

from scripts.brucite_joint_pv_errors import REPORT as LINEAR_REPORT
from scripts.brucite_quadratic_pv_errors import (
    REPORT,
    jacobian,
    moduli,
    polynomial_transform,
    read_data,
    reproduce,
    residual,
)
from scripts.brucite_reproduce import bm3


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_quadratic_analytic_jacobian_matches_numeric_derivatives():
    _, data = read_data()
    x = np.r_[[40, 30, 25, 7, 8, 10], 0.3 * np.sin(np.arange(37))]
    expected = []
    for i, value in enumerate(x):
        step = 2e-5 * max(abs(value), 1)
        delta = np.zeros_like(x)
        delta[i] = step
        expected.append(
            (residual(x + delta, data, 2, True) - residual(x - delta, data, 2, True))
            / (2 * step)
        )
    np.testing.assert_allclose(
        jacobian(x, data, 2, True), np.array(expected).T, rtol=2e-7, atol=1e-7
    )


def test_polynomial_conversion_uses_full_squared_term_and_fixed_ambient(report):
    linear = json.loads(LINEAR_REPORT.read_text())
    assert report["fixed_ambient_law"] == linear["fixed_ambient_law"]
    assert report["fits"]["linear_pv_errors"] == linear["fits"]["joint_pv_errors"]
    fit = report["fits"]["quadratic_pv_errors"]
    beta = np.array(fit["bernstein_coefficients"])
    coefficients = polynomial_transform() @ beta
    t = np.linspace(300, 873, 31)
    dt = t - 300
    k, kp = moduli(beta, t)
    np.testing.assert_allclose(
        k, coefficients[0] + coefficients[1] * dt + coefficients[2] * dt**2
    )
    np.testing.assert_allclose(
        kp, coefficients[3] + coefficients[4] * dt + coefficients[5] * dt**2
    )
    assert np.min(k) > 0 and np.min(kp) > 4
    ambient = report["fixed_ambient_law"]
    v0 = ambient["V300_a3"] * np.exp(ambient["alpha_per_k"] * dt)
    np.testing.assert_allclose(bm3(v0, v0, k, kp), 0, atol=1e-12)


def test_exact_chi_square_and_degrees_of_freedom_and_all_rows(report):
    fit = report["fits"]["quadratic_pv_errors"]
    assert fit["solver_success"]
    assert report["source_rows"] == list(range(1, 38))
    assert [group["observations"] for group in report["per_temperature"]] == [
        13,
        8,
        8,
        8,
    ]
    assert fit["global_parameters"] == 6
    assert (
        fit["residual_components"] - fit["latent_volumes"] - 6
        == fit["degrees_of_freedom"]
        == 31
    )
    assert fit["active_global_bounds"] == [0] * 6
    p = fit["projection"]
    rp, rv, distance = (
        np.array(p[key])
        for key in [
            "pressure_residual_sigma",
            "volume_correction_sigma",
            "signed_orthogonal_distance_sigma",
        ]
    )
    np.testing.assert_allclose(distance**2, rp**2 + rv**2, atol=1e-9)
    assert rp @ rp + rv @ rv == pytest.approx(fit["chi_square"], abs=1e-8)
    assert fit["chi_square"] == pytest.approx(98.0932030423, abs=1e-7)
    assert fit["reduced_chi_square"] == pytest.approx(fit["chi_square"] / 31)
    assert fit["reduced_chi_square"] > 3
    assert all(run["success"] for run in fit["multistart"])
    assert np.ptp([run["objective"] for run in fit["multistart"]]) < 1e-7


def test_profiled_covariance_matches_full_latent_inverse(report):
    _, data = read_data()
    fit = report["fits"]["quadratic_pv_errors"]
    x = np.r_[
        fit["bernstein_coefficients"], fit["projection"]["volume_correction_sigma"]
    ]
    j = jacobian(x, data, 2, True)
    full = np.linalg.inv(np.einsum("ni,nj->ij", j, j))
    np.testing.assert_allclose(
        fit["bernstein_covariance"], full[:6, :6], rtol=3e-6, atol=1e-9
    )
    transform = polynomial_transform()
    physical = np.einsum("ij,jk,lk->il", transform, full[:6, :6], transform)
    np.testing.assert_allclose(fit["covariance"], physical, rtol=3e-6, atol=1e-9)


def test_improvement_and_independent_checks_remain_qualified(report):
    quadratic = report["fits"]["quadratic_pv_errors"]
    linear = report["fits"]["linear_pv_errors"]
    assert quadratic["chi_square"] < linear["chi_square"]
    assert quadratic["reduced_chi_square"] < linear["reduced_chi_square"]
    assert quadratic["volume_rmse_a3"] == pytest.approx(0.0775231219, abs=1e-8)
    for group in report["per_temperature"]:
        assert (
            group["comparisons"]["quadratic_pv_errors"]["volume_rmse_a3"]
            < group["comparisons"]["linear_pv_errors"]["volume_rmse_a3"]
        )
    checked = report["verification"]
    assert abs(checked["independent_odrpack"]["objective_difference"]) < 1e-5
    assert checked["native_bm3"]["max_pressure_difference_gpa"] < 1e-10
    assert checked["native_bm3"]["max_volume_difference_a3"] < 1e-8
    assert checked["native_bm3"]["minimum_sampled_bulk_modulus_gpa"] > 0
    control = report["independent_isotherm_pv_error_control"]
    assert control["degrees_of_freedom"] == 29
    assert (
        sum(fit["chi_square"] for fit in control["isotherms"]) == control["chi_square"]
    )
    assert control["reduced_chi_square"] > 2
    assert control["chi_square"] < quadratic["chi_square"]
    assert report == json.loads(REPORT.read_text())
