"""Tests for the natural-strain equations of state."""

import numpy as np
import pytest

from peritheos.eos.rt import NaturalStrain2, NaturalStrain3, NaturalStrain4
from peritheos.utils import derivative


@pytest.mark.parametrize(
    "eos",
    [
        NaturalStrain2(10.0, 100.0),
        NaturalStrain3(10.0, 100.0, 4.0),
        NaturalStrain4(10.0, 100.0, 4.0, -0.01),
    ],
)
def test_reference_state_and_array_support(eos):
    volumes = eos.V0 * np.array([0.8, 1.0, 1.2])

    assert eos.pressure(eos.V0) == 0.0
    assert eos.bulk_modulus(eos.V0) == eos.K0
    assert eos.pressure(volumes).shape == volumes.shape
    assert eos.pressure(volumes)[0] > 0.0
    assert eos.pressure(volumes)[2] < 0.0


@pytest.mark.parametrize(
    "eos",
    [
        NaturalStrain2(10.0, 100.0),
        NaturalStrain3(10.0, 100.0, 4.0),
        NaturalStrain4(10.0, 100.0, 4.0, -0.01),
    ],
)
def test_bulk_modulus_matches_pressure_derivative(eos):
    volumes = eos.V0 * np.array([0.8, 0.9, 1.0, 1.1])

    assert np.allclose(
        eos.bulk_modulus(volumes),
        -volumes * derivative(eos.pressure, volumes),
    )


def test_second_order_implies_k_prime_of_two():
    eos = NaturalStrain2(10.0, 100.0)
    step = 1.0e-4
    volumes = eos.volume(np.array([-step, step]))
    numerical = np.diff(eos.bulk_modulus(volumes)) / (2.0 * step)

    assert np.isclose(numerical[0], 2.0, rtol=1.0e-6)


def test_fourth_order_recovers_both_bulk_modulus_derivatives():
    eos = NaturalStrain4(10.0, 100.0, 4.2, -0.015)
    pressures = np.linspace(-0.02, 0.02, 21)
    moduli = eos.bulk_modulus(eos.volume(pressures))
    coefficients = np.polynomial.polynomial.polyfit(pressures, moduli, 3)

    assert np.isclose(coefficients[1], eos.K0_prime, rtol=1.0e-7)
    assert np.isclose(2.0 * coefficients[2], eos.K0_double_prime, rtol=1.0e-5)
