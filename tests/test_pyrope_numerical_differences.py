"""Check native-kernel attribution and the accurate analytic replay prototype."""

import hashlib
import json

import numpy as np
import pytest

from scripts.audit_pyrope_author_inputs import DATA, NAMES
from scripts.audit_pyrope_numerical_differences import (
    PRECISE,
    SCALE,
    STEM,
    coarse_pressure,
    covariance,
    forward_sigma,
    inputs,
)
from scripts.refit_pyrope_angel_2022_subset import effective_sigma
from scripts.refit_pyrope_author_analytic import evaluate, fit


def test_native_capture_closes_covariance_and_gauss_newton_step():
    rows, d, native, evidence, original = inputs()
    assert len(native) == len(rows) == 194
    assert evidence["native_exit_status"] == 0
    assert evidence["debugged_final_log_parameters_match_uninstrumented_replay"]
    assert (
        hashlib.sha256((DATA / (STEM + "-native-kernel.csv")).read_bytes()).hexdigest()
        == evidence["native_kernel_csv_sha256"]
    )
    j = np.array([[r["dP_d_" + name] for name in NAMES] for r in native])
    w = np.array([r["native_weight"] for r in native])
    cov = covariance(j, w, evidence["reduced_chi_square"])
    ix = evidence["covariance_parameter_indices_zero_based"]
    np.testing.assert_allclose(
        np.sqrt(cov.diagonal()),
        np.array(evidence["errors59_after_step"])[ix],
        rtol=1e-7,
    )
    report = json.loads((DATA / (STEM + "-analysis.json")).read_text())
    assert (
        report["native_step_reconstruction"]["maximum_difference_in_float32_ulps"] < 2
    )
    assert abs(report["native_report_minus_reconstructed_cycle_chi_square"]) < 2e-6
    total_python = (
        sum(g["python_weighted_sum"] for g in report["dataset_contributions"]) / 189
    )
    total_native = (
        sum(g["native_weighted_sum"] for g in report["dataset_contributions"]) / 189
    )
    np.testing.assert_allclose(
        total_python, report["objective_component_swaps"][0]["reduced_sum_squares"]
    )
    np.testing.assert_allclose(
        total_native, report["objective_component_swaps"][-1]["reduced_sum_squares"]
    )
    assert original["independent_fits"][0]["degrees_of_freedom"] == 189


def test_finite_error_propagation_retains_curvature_and_differs_from_derivative():
    # P(Y,T)=Y^2+T^2 has exact derivatives and finite increments known by hand.
    d = {
        "Y": np.array([10.0]),
        "T": np.array([5.0]),
        "sY": np.array([1.0]),
        "sT": np.array([2.0]),
        "sP": np.array([3.0]),
    }

    def function(c, d, *, y=None, t=None):
        return (d["Y"] if y is None else y) ** 2 + (d["T"] if t is None else t) ** 2

    propagated = forward_sigma(None, d, function)[0]
    assert propagated == pytest.approx(np.sqrt(3**2 + 21**2 + 24**2))
    assert propagated != pytest.approx(np.sqrt(3**2 + 20**2 + 20**2))


def test_tightened_modulus_inversion_recovers_precise_pressure():
    _, d, _, evidence, _ = inputs()
    c = np.array(list(evidence["starting_parameters"].values()))
    precise = PRECISE(c, d)
    np.testing.assert_allclose(
        coarse_pressure(c, d, tight=True), precise, atol=1e-8, rtol=0
    )
    assert np.max(abs(coarse_pressure(c, d) - precise)) > 0.009


def test_analytic_parameter_and_coordinate_derivatives_match_five_point_checks():
    _, d, _, evidence, _ = inputs()
    c = np.array(list(evidence["starting_parameters"].values()))
    p, j, py, pt = evaluate(c, d)
    np.testing.assert_allclose(p, PRECISE(c, d), atol=0, rtol=0)
    for i, scale in enumerate(SCALE):
        h = scale * 1e-4
        values = []
        for multiple in [-2, -1, 1, 2]:
            trial = c.copy()
            trial[i] += multiple * h
            values.append(PRECISE(trial, d))
        check = (values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (12 * h)
        np.testing.assert_allclose(j[:, i], check, rtol=2e-6, atol=2e-9)
    for arg, key, h, derivative in [("y", "Y", 1e-3, py), ("t", "T", 0.1, pt)]:
        values = [
            PRECISE(c, d, **{arg: d[key] + multiple * h}) for multiple in [-2, -1, 1, 2]
        ]
        check = (values[0] - 8 * values[1] + 8 * values[2] - values[3]) / (12 * h)
        np.testing.assert_allclose(derivative, check, rtol=2e-7, atol=2e-10)
    sigma = np.sqrt(d["sP"] ** 2 + (py * d["sY"]) ** 2 + (pt * d["sT"]) ** 2)
    np.testing.assert_allclose(sigma, effective_sigma(c, d), rtol=1e-7)


def test_analytic_replay_preserves_parameters_errors_and_stationarity():
    rows, _, _, _, original = inputs()
    baseline = original["independent_fits"][0]
    result = fit(rows, [170, 4.5, 800, 1.2, 0.98])
    assert result["numerically_converged"] and not result["at_parameter_bound"]
    assert result["maximum_normalized_stationarity_score"] < 1e-6
    np.testing.assert_allclose(
        list(result["parameters"].values()),
        list(baseline["parameters"].values()),
        rtol=1e-7,
    )
    np.testing.assert_allclose(
        list(result["conditional_residual_scaled_widths"].values()),
        list(baseline["conditional_residual_scaled_widths"].values()),
        rtol=1e-7,
    )
    assert result["conditional_reduced_chi_square"] == pytest.approx(
        baseline["conditional_reduced_chi_square"], abs=1e-8
    )
