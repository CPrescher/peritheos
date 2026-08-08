"""Tests for P-V and P-V-T EOS fitting."""

import numpy as np
import pytest

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenEinstein
from peritheos.fitting import fit_rt_eos, fit_thermal_eos


def test_fit_rt_eos_recovers_synthetic_parameters():
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(7.5, 11.0, 30)
    pressures = expected.pressure(volumes)
    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"V0": 9.8, "K0": 110.0, "K0_prime": 4.0},
        bounds={"V0": (8.0, 12.0), "K0": (50.0, 200.0)},
    )

    assert result.success
    assert np.isclose(result.parameters["V0"], expected.V0)
    assert np.isclose(result.parameters["K0"], expected.K0)
    assert np.isclose(result.parameters["K0_prime"], expected.K0_prime)
    assert result.covariance.shape == (3, 3)
    assert result.correlation.shape == (3, 3)
    assert result.degrees_of_freedom == volumes.size - 3


def test_fit_with_fixed_parameter_and_absolute_uncertainties():
    expected = BM3(10.0, 120.0, 4.3)
    volumes = np.linspace(8.0, 10.5, 20)
    pressures = expected.pressure(volumes)
    result = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        initial={"K0": 100.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        sigma=0.05,
        absolute_sigma=True,
    )

    assert result.parameters["V0"] == 10.0
    assert result.standard_errors["V0"] == 0.0
    assert np.all(np.diag(result.covariance) > 0.0)


def test_bm2_fit_matches_closed_form_weighted_least_squares():
    """Benchmark the optimizer against an independently solvable linear fit."""
    volumes = np.array([8.0, 8.5, 9.0, 9.5, 10.0])
    observed_pressure = np.array([35.2, 23.8, 14.5, 6.3, 0.2])
    pressure_sigma = np.array([0.4, 0.3, 0.25, 0.2, 0.2])
    unit_pressure = np.asarray(BM3(10.0, 1.0, 4.0).pressure(volumes))
    weights = pressure_sigma**-2
    expected_k0 = np.sum(weights * unit_pressure * observed_pressure) / np.sum(
        weights * unit_pressure**2
    )
    expected_variance = 1.0 / np.sum(weights * unit_pressure**2)

    result = fit_rt_eos(
        BM3,
        volumes,
        observed_pressure,
        initial={"K0": 90.0},
        fixed={"V0": 10.0, "K0_prime": 4.0},
        pressure_sigma=pressure_sigma,
        absolute_sigma=True,
    )

    assert np.isclose(result.parameters["K0"], expected_k0, rtol=1.0e-8)
    assert np.isclose(result.covariance[0, 0], expected_variance, rtol=1.0e-6)


def test_fit_rt_eos_handles_pressure_and_volume_uncertainties():
    expected = BM3(10.0, 120.0, 4.3)
    true_volumes = np.linspace(8.0, 10.5, 20)
    measured_volumes = true_volumes + 0.008 * np.sin(np.arange(true_volumes.size))
    pressures = expected.pressure(true_volumes)

    result = fit_rt_eos(
        BM3,
        measured_volumes,
        pressures,
        initial={"K0": 110.0, "K0_prime": 4.0},
        fixed={"V0": 10.0},
        pressure_sigma=0.002,
        volume_sigma=0.01,
        absolute_sigma=True,
    )

    assert result.success
    assert np.isclose(result.parameters["K0"], expected.K0, rtol=0.01)
    assert np.isclose(result.parameters["K0_prime"], expected.K0_prime, rtol=0.02)
    assert np.any(np.abs(result.volume_corrections) > 0.0)
    assert result.adjusted_temperature is None
    assert result.temperature_corrections is None
    assert result.weighted_residuals.size == 2 * pressures.size
    assert result.degrees_of_freedom == pressures.size - 2


def test_fit_thermal_eos_recovers_gamma():
    rt_eos = BM3(1.0, 160.0, 4.0)
    expected = MieGruneisenEinstein(rt_eos, 300.0, 800.0, 1.6, 1.0, 2.0)
    volumes = np.repeat(np.array([0.8, 0.9, 1.0]), 5)
    temperatures = np.tile(np.linspace(500.0, 2000.0, 5), 3)
    pressures = expected.pressure(volumes, temperatures)
    result = fit_thermal_eos(
        MieGruneisenEinstein,
        rt_eos,
        volumes,
        temperatures,
        pressures,
        initial={"gamma0": 1.3, "q": 0.8},
        fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
    )

    assert result.success
    assert np.isclose(result.parameters["gamma0"], 1.6)
    assert np.isclose(result.parameters["q"], 1.0)


def test_fit_thermal_eos_handles_uncertainties_in_all_observables():
    rt_eos = BM3(1.0, 160.0, 4.0)
    expected = MieGruneisenEinstein(rt_eos, 300.0, 800.0, 1.6, 1.0, 2.0)
    true_volumes = np.repeat(np.array([0.8, 0.9, 1.0]), 5)
    true_temperatures = np.tile(np.linspace(500.0, 2000.0, 5), 3)
    indices = np.arange(true_volumes.size)
    measured_volumes = true_volumes + 0.0003 * np.sin(indices)
    measured_temperatures = true_temperatures + 1.5 * np.cos(indices)
    pressures = expected.pressure(true_volumes, true_temperatures)

    result = fit_thermal_eos(
        MieGruneisenEinstein,
        rt_eos,
        measured_volumes,
        measured_temperatures,
        pressures,
        initial={"gamma0": 1.4, "q": 0.8},
        fixed={"Tr": 300.0, "theta0": 800.0, "n": 2.0},
        pressure_sigma=0.001,
        volume_sigma=0.0005,
        temperature_sigma=2.0,
        absolute_sigma=True,
    )

    assert result.success
    assert np.isclose(result.parameters["gamma0"], 1.6, rtol=0.01)
    assert np.isclose(result.parameters["q"], 1.0, rtol=0.04)
    assert np.any(np.abs(result.volume_corrections) > 0.0)
    assert np.any(np.abs(result.temperature_corrections) > 0.0)
    assert result.adjusted_temperature is not None
    assert result.weighted_residuals.size == 3 * pressures.size
    assert result.degrees_of_freedom == pressures.size - 2


def test_fit_input_validation():
    volumes = np.linspace(8.0, 10.0, 5)
    pressures = BM3(10.0, 120.0, 4.0).pressure(volumes)
    with pytest.raises(ValueError, match="both initial and fixed"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"V0": 10.0},
            fixed={"V0": 10.0, "K0": 120.0, "K0_prime": 4.0},
        )
    with pytest.raises(ValueError, match="sigma"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            sigma=0.0,
        )
    with pytest.raises(ValueError, match="either pressure_sigma or sigma"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            pressure_sigma=0.1,
            sigma=0.1,
        )
    with pytest.raises(ValueError, match="volume_sigma"):
        fit_rt_eos(
            BM3,
            volumes,
            pressures,
            initial={"K0": 110.0},
            fixed={"V0": 10.0, "K0_prime": 4.0},
            volume_sigma=-0.1,
        )
