"""Tests for density--velocity finite-strain models and fitting."""

import json

import numpy as np
import pytest

from peritheos.acoustics import EulerianFiniteStrainAcoustic
from peritheos.errors import EosValidationError, FitValidationError
from peritheos.fitting import fit_acoustic_finite_strain


def test_acoustic_model_recovers_reference_moduli_and_velocity_identity():
    model = EulerianFiniteStrainAcoustic(4.11, 252.0, 4.1, 175.0, 1.66)
    density = np.array([4.11, 4.30, 4.45])
    vp, vs = model.velocities(density)
    bulk, shear = model.moduli(density)

    assert bulk[0] == pytest.approx(252.0)
    assert shear[0] == pytest.approx(175.0)
    assert density * vs**2 == pytest.approx(shear)
    assert density * (vp**2 - 4.0 * vs**2 / 3.0) == pytest.approx(bulk)


def test_acoustic_fit_recovers_exact_synthetic_coefficients():
    expected = EulerianFiniteStrainAcoustic(4.11, 247.0, 4.5, 176.0, 1.6)
    density = np.linspace(4.13, 4.5, 12)
    vp, vs = expected.velocities(density)

    result = fit_acoustic_finite_strain(
        density,
        vp,
        vs,
        rho0=expected.rho0,
        initial={
            "K_S0": 240.0,
            "K_S0_prime": 4.0,
            "G0": 170.0,
            "G0_prime": 2.0,
        },
    )

    assert result.success
    assert result.parameters == pytest.approx(
        {
            "K_S0": expected.K_S0,
            "K_S0_prime": expected.K_S0_prime,
            "G0": expected.G0,
            "G0_prime": expected.G0_prime,
        },
        rel=2.0e-7,
    )
    assert result.adjusted_density == pytest.approx(density)
    assert result.compressional_velocity_residuals == pytest.approx(0.0, abs=1.0e-8)
    assert result.shear_velocity_residuals == pytest.approx(0.0, abs=1.0e-8)
    assert json.loads(result.to_json())["model"]["parameters"]["rho0"] == 4.11


def test_acoustic_fit_supports_density_eiv_and_correlated_covariance():
    expected = EulerianFiniteStrainAcoustic(4.11, 247.0, 4.5, 176.0, 1.6)
    density = np.linspace(4.13, 4.5, 10)
    vp, vs = expected.velocities(density)
    sigmas = np.array([0.06, 0.04, 0.004])
    correlation = np.array(
        [
            [1.0, 0.20, -0.10],
            [0.20, 1.0, 0.15],
            [-0.10, 0.15, 1.0],
        ]
    )
    covariance = correlation * np.outer(sigmas, sigmas)

    result = fit_acoustic_finite_strain(
        density,
        vp,
        vs,
        rho0=expected.rho0,
        initial={
            "K_S0": 245.0,
            "K_S0_prime": 4.4,
            "G0": 174.0,
            "G0_prime": 1.7,
        },
        observation_covariance=covariance,
        absolute_sigma=True,
    )

    assert result.success
    assert result.degrees_of_freedom == 16
    assert result.parameters == pytest.approx(
        {
            "K_S0": expected.K_S0,
            "K_S0_prime": expected.K_S0_prime,
            "G0": expected.G0,
            "G0_prime": expected.G0_prime,
        },
        rel=2.0e-7,
    )
    assert result.covariance.shape == (4, 4)
    assert np.all(np.isfinite(result.correlation))


@pytest.mark.parametrize(
    ("model_kwargs", "message"),
    [
        ({"rho0": 0.0}, "rho0 must be greater than zero"),
        ({"K_S0": -1.0}, "K_S0 must be greater than zero"),
    ],
)
def test_acoustic_model_rejects_invalid_parameters(model_kwargs, message):
    kwargs = {
        "rho0": 4.11,
        "K_S0": 247.0,
        "K_S0_prime": 4.5,
        "G0": 176.0,
        "G0_prime": 1.6,
        **model_kwargs,
    }
    with pytest.raises(EosValidationError, match=message):
        EulerianFiniteStrainAcoustic(**kwargs)


def test_acoustic_fit_rejects_incomplete_or_double_uncertainty_models():
    initial = {"K_S0": 247.0, "K_S0_prime": 4.5, "G0": 176.0}
    with pytest.raises(FitValidationError, match="Missing acoustic parameters"):
        fit_acoustic_finite_strain([4.2], [11.0], [6.6], rho0=4.11, initial=initial)

    initial["G0_prime"] = 1.6
    with pytest.raises(FitValidationError, match="cannot be combined"):
        fit_acoustic_finite_strain(
            [4.2],
            [11.0],
            [6.6],
            rho0=4.11,
            initial=initial,
            compressional_velocity_sigma=0.06,
            observation_covariance=np.eye(2),
        )
