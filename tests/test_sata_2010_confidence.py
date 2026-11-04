"""Covariance geometry and mean-curve propagation for Sata diagnostics."""

import json

import numpy as np
import pytest
from scipy.stats import f as f_distribution
from scipy.stats import t as t_distribution

from scripts.audit_sata_2010_fe3c import VR
from scripts.plot_sata_2010_confidence import (
    FIT_REPORT,
    LEVELS,
    REPORT,
    VERSIONS,
    build,
    ellipse_points,
    fit_state,
    mean_and_gradient,
)


@pytest.fixture(scope="module")
def calculation():
    return build()


def test_joint_pair_ellipses_have_requested_mahalanobis_radius(calculation):
    report, _, _ = calculation
    for version in VERSIONS:
        for model, entry in report["fits"][version].items():
            cov = np.array(entry["pair_covariance"])
            center = np.array(entry["parameters"])[
                [0, 1] if model == "fixed_V0_BM3" else [1, 2]
            ]
            nu = entry["residual_degrees_of_freedom"]
            for level in LEVELS.values():
                _, points, q = ellipse_points(center, cov, level, nu)
                delta = points - center
                distance = np.einsum("ij,jk,ik->i", delta, np.linalg.inv(cov), delta)
                np.testing.assert_allclose(distance, q, atol=1e-10)
                assert f_distribution.cdf(q / 2, 2, nu) == pytest.approx(level)
                np.testing.assert_allclose(points[0], points[-1], atol=1e-12)
            assert entry["K_Kprime_correlation"] < -0.97


def test_mean_curve_bands_match_independent_linear_coefficient_propagation(calculation):
    _, _, bands = calculation
    source = json.loads(FIT_REPORT.read_text())
    measured_v = np.array([r["volume_a3"] for r in source["inputs"]])

    def design(v, model):
        f = ((VR / v) ** (2 / 3) - 1) / 2
        columns = [3 * f - 18 * f**2, 4.5 * f**2]
        if model == "reference_mBM3":
            columns.insert(0, 1 - 5 * f + 17.5 * f**2)
        return (1 + 2 * f)[:, None] ** 2.5 * np.column_stack(columns)

    for version in VERSIONS:
        observed = np.array(
            [
                r[
                    "published_pressure_gpa"
                    if version == "original_pressures"
                    else "fit_pressure_gpa"
                ]
                for r in source["inputs"]
            ]
        )
        for model in ("fixed_V0_BM3", "reference_mBM3"):
            rows = [
                r
                for r in bands
                if r["pressure_version"] == version and r["model"] == model
            ]
            v = np.array([r["volume_a3"] for r in rows])
            x = design(measured_v, model)
            beta = np.linalg.lstsq(x, observed, rcond=None)[0]
            dof = len(observed) - len(beta)
            s2 = np.sum((x @ beta - observed) ** 2) / dof
            covariance = s2 * np.linalg.inv(x.T @ x)
            grid = design(v, model)
            expected_se = np.sqrt(np.einsum("ij,jk,ik->i", grid, covariance, grid))
            se = np.array([r["mean_standard_error_pressure_gpa"] for r in rows])
            np.testing.assert_allclose(se, expected_se, rtol=1e-9)
            for level_name, level in LEVELS.items():
                width = np.array(
                    [
                        r[f"upper_{level_name}_pressure_gpa"] - r["mean_pressure_gpa"]
                        for r in rows
                    ]
                )
                np.testing.assert_allclose(
                    width,
                    t_distribution.ppf((1 + level) / 2, dof) * expected_se,
                    atol=1e-10,
                )
                # A future-observation interval would include s2; these do not.
                assert np.all(
                    width
                    < t_distribution.ppf((1 + level) / 2, dof)
                    * np.sqrt(expected_se**2 + s2)
                )


def test_full_physical_parameter_gradients_and_existing_fit_preservation(calculation):
    source = json.loads(FIT_REPORT.read_text())
    volumes = np.array([108.0, 120.0, 137.0])
    for model in ("fixed_V0_BM3", "reference_mBM3"):
        parameters, _, _ = fit_state(source, VERSIONS[1], model)
        _, gradient = mean_and_gradient(volumes, parameters, model)
        for i in range(len(parameters)):
            h = 1e-4 if i < len(parameters) - 1 else 1e-5
            step = np.zeros_like(parameters)
            step[i] = h
            finite = (
                mean_and_gradient(volumes, parameters + step, model)[0]
                - mean_and_gradient(volumes, parameters - step, model)[0]
            ) / (2 * h)
            np.testing.assert_allclose(gradient[:, i], finite, rtol=1e-7, atol=1e-8)
    report, _, _ = calculation
    assert report == json.loads(REPORT.read_text())
    assert "full 3x3 covariance" in report["band_convention"]
    assert (
        report["fits"][VERSIONS[1]]["reference_mBM3"]["residual_degrees_of_freedom"]
        == 9
    )
