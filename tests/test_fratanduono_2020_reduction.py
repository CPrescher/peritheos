"""Check the qualified reconstruction without promoting it to parameter parity."""

import json

import numpy as np
import pytest

from scripts.reproduce_fratanduono_2020_cu import (
    DATA,
    PREFIX,
    R_SPECIFIC,
    RHO0,
    SWITCH_DENSITY,
    THETA0,
    TR,
    debye_energy,
    kraus_gamma,
    kraus_temperature_factor,
    reconstruct,
    reproduce,
    thermal_correction,
)


def test_thermal_reference_and_thermodynamic_derivative():
    assert kraus_temperature_factor(RHO0) == 1
    assert kraus_gamma(RHO0) == 2
    np.testing.assert_array_equal(thermal_correction([RHO0], [2], [1]), [0])
    density = np.array([10.101, 12.195, SWITCH_DENSITY])
    step = 1e-5
    numerical_gamma = (
        np.log(kraus_temperature_factor(density * np.exp(step)))
        - np.log(kraus_temperature_factor(density * np.exp(-step)))
    ) / (2 * step)
    np.testing.assert_allclose(numerical_gamma, kraus_gamma(density), rtol=1e-8)
    # Dulong-Petit limit independently checks energy units and normalization.
    assert debye_energy(1e8, THETA0) / (3 * R_SPECIFIC * 1e8) == pytest.approx(
        1, rel=2e-6
    )
    factor = kraus_temperature_factor(density)
    correction = thermal_correction(density, kraus_gamma(density), factor)
    assert np.all(correction > 0)
    assert np.all(TR * factor > TR)


def test_reconstruction_stops_at_figure_coverage_and_converges():
    source, result = reconstruct()
    _, refined = reconstruct(grid_points=16385)
    assert source.shape == (100, 3)
    np.testing.assert_array_equal(source[0], [0, 0, 8.938])
    np.testing.assert_array_equal(source[-1], [2331.8, 173.8, 27.13])
    assert len(result["density"]) == 53
    assert result["low"].sum() == 7
    assert result["density"][-1] == 22.521
    assert result["density"][-1] <= result["hugoniot_max_density"] < source[54, 2]
    np.testing.assert_allclose(
        result["isotherm"], refined["isotherm"], atol=1e-5, rtol=0
    )
    assert np.all(result["isotherm"] < result["isentrope"])
    assert np.all(np.diff(result["isotherm"]) > 0)
    # Check the last fully specified Altshuler row against independent quadrature.
    assert result["correction"][6] == pytest.approx(1.76014977919467)


def test_audit_records_limits_and_does_not_claim_parameter_parity(assert_audit_close):
    result, rows = reproduce()
    saved = json.loads((DATA / f"{PREFIX}-thermal-reduction.json").read_text())
    assert_audit_close(result, saved)
    assert len(rows) == 53
    assert result["selection"]["unreconstructed_high_density_rows"] == 46
    assert result["figure_blue_check"]["negative_implied_thermal_correction_rows"] == 16
    assert result["curve_checks"][
        "candidate_combined_298k_rmse_vs_published_gpa"
    ] == pytest.approx(1.56222884, rel=1e-6)
    for name in ["kraus_branch_298k_diagnostic", "candidate_combined_298k_diagnostic"]:
        fit = result["fits"][name]
        assert fit["solver_success"]
        assert not any(fit["within_printed_parameter_errors"].values())
        assert fit["parameter_covariance"] is None
        assert fit["objective"] == "unweighted pressure residuals"
