"""Check exact two-coordinate residuals, profiled uncertainty and joint replay."""

import json

import numpy as np
import pytest

from scripts.brucite_joint_pv_errors import (
    ORDER,
    REPORT,
    SCALE,
    STAGED_REPORT,
    jacobian,
    read_data,
    reproduce,
    residual,
)
from scripts.brucite_replay_checks import check_report
from scripts.brucite_staged_continuous import runtime_model


@pytest.fixture(scope="module")
def report():
    return reproduce()


def test_analytic_jacobian_for_both_coordinate_errors():
    _, data = read_data()
    x = np.r_[[40, 7, -20, 1.5], 0.3 * np.sin(np.arange(37))]
    expected = []
    for i, value in enumerate(x):
        step = 2e-5 * max(abs(value), 1)
        delta = np.zeros_like(x)
        delta[i] = step
        expected.append(
            (residual(x + delta, data, True) - residual(x - delta, data, True))
            / (2 * step)
        )
    np.testing.assert_allclose(
        jacobian(x, data, True), np.array(expected).T, rtol=2e-7, atol=1e-7
    )


def test_all_rows_fixed_ambient_law_and_exact_chi_square(report):
    staged = json.loads(STAGED_REPORT.read_text(encoding="utf-8"))
    assert report["fixed_ambient_law"] == staged["stage1_ambient"]
    assert report["source_rows"] == list(range(1, 38))
    assert [group["observations"] for group in report["per_temperature"]] == [
        13,
        8,
        8,
        8,
    ]
    fit = report["fits"]["joint_pv_errors"]
    assert fit["solver_success"]
    assert fit["active_global_bounds"] == [0, 0, 0, 0]
    assert (
        fit["residual_components"] - fit["latent_volumes"] - fit["global_parameters"]
        == fit["degrees_of_freedom"]
        == 33
    )
    projection = fit["projection"]
    rp, rv, distance = (
        np.array(projection[key])
        for key in [
            "pressure_residual_sigma",
            "volume_correction_sigma",
            "signed_orthogonal_distance_sigma",
        ]
    )
    np.testing.assert_allclose(distance**2, rp**2 + rv**2, atol=1e-10)
    assert rp @ rp + rv @ rv == pytest.approx(fit["chi_square"], abs=1e-8)
    assert np.sum(distance**2) / 33 == pytest.approx(
        fit["reduced_chi_square"], abs=1e-9
    )
    assert fit["chi_square"] == pytest.approx(173.407569012, abs=1e-7)
    assert fit["reduced_chi_square"] > 5  # Preserve the unresolved fit discrepancy.
    assert all(run["success"] for run in fit["multistart"])
    assert np.ptp([run["objective"] for run in fit["multistart"]]) < 1e-7


def test_profiled_covariance_matches_full_latent_parameter_inverse(report):
    _, data = read_data()
    fit = report["fits"]["joint_pv_errors"]
    x = np.r_[
        np.array([fit["parameters"][key] for key in ORDER]) / SCALE,
        fit["projection"]["volume_correction_sigma"],
    ]
    j = jacobian(x, data, True)
    full_covariance = np.linalg.inv(np.einsum("ni,nj->ij", j, j))
    expected = full_covariance[:4, :4] * SCALE[:, None] * SCALE[None, :]
    np.testing.assert_allclose(fit["covariance"], expected, rtol=2e-6, atol=1e-10)
    np.testing.assert_allclose(
        np.diag(fit["covariance"]),
        np.array([fit["standard_errors"][key] for key in ORDER]) ** 2,
    )
    np.testing.assert_allclose(
        fit["covariance_scaled_by_observed_scatter"],
        np.array(fit["covariance"]) * fit["reduced_chi_square"],
    )


def test_public_native_fit_agrees_and_evaluation_preserves_ambient_volumes(report):
    checked = report["verification"]["native_public_fit"]
    assert checked["success"]
    assert abs(checked["objective_difference"]) < 2e-5
    assert checked["maximum_parameter_difference_in_native_standard_errors"] < 0.002
    assert report["verification"]["runtime_max_pressure_difference_gpa"] < 1e-10
    ambient = report["fixed_ambient_law"]
    parameters = report["fits"]["joint_pv_errors"]["parameters"]
    eos = runtime_model(
        {"fixed_ambient_law": ambient, "linear_moduli": {"parameters": parameters}}
    )
    temperatures = np.array([299.3, 300, 400, 473, 573, 673, 773, 873])
    v0 = ambient["V300_a3"] * np.exp(ambient["alpha_per_k"] * (temperatures - 300))
    np.testing.assert_allclose(eos.pressure(v0, temperatures), 0, atol=1e-10)
    for t in [300, 400, 573, 673, 773, 873]:
        volume = eos.volume(10, t)
        assert eos.pressure(volume, t) == pytest.approx(10, abs=1e-9)


def test_different_objectives_and_comparators_remain_distinct(report):
    weighted = report["fits"]["joint_pv_errors"]
    baseline = report["fits"]["unweighted_joint"]
    assert weighted["projection"]["chi_square"] < baseline["projection"]["chi_square"]
    # Minimizing the normalized two-coordinate distance is a different task.
    assert weighted["volume_rmse_a3"] > baseline["volume_rmse_a3"]
    for fit in report["fits"].values():
        assert len(fit["row_pressure_residuals_gpa"]) == 37
        assert len(fit["row_volume_residuals_a3"]) == 37
        assert len(fit["projection"]["adjusted_volume_a3"]) == 37
    check_report(report, json.loads(REPORT.read_text(encoding="utf-8")))
