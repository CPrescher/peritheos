"""Tests for the Murnaghan equation of state."""

import numpy as np
import pytest

from peritheos.eos.rt import Murnaghan
from peritheos.utils import derivative


@pytest.fixture
def eos():
    return Murnaghan(V0=10.0, K0=100.0, K0_prime=4.0)


def test_parameters(eos):
    assert eos.V0 == 10.0
    assert eos.K0 == 100.0
    assert eos.K0_prime == 4.0


def test_reference_state(eos):
    assert eos.pressure(eos.V0) == 0.0
    assert eos.bulk_modulus(eos.V0) == eos.K0


def test_scalar_and_array_pressure(eos):
    volumes = eos.V0 * np.array([0.8, 1.0, 1.2])
    pressures = eos.pressure(volumes)

    assert pressures.shape == volumes.shape
    assert pressures[0] > 0.0
    assert pressures[1] == 0.0
    assert pressures[2] < 0.0


def test_known_value(eos):
    volume = 0.8 * eos.V0
    expected = eos.K0 / eos.K0_prime * ((eos.V0 / volume) ** 4.0 - 1.0)

    assert np.isclose(eos.pressure(volume), expected)
    assert np.isclose(eos.bulk_modulus(volume), eos.K0 * (eos.V0 / volume) ** 4.0)


def test_bulk_modulus_matches_pressure_derivative(eos):
    volumes = eos.V0 * np.array([0.8, 0.9, 1.0, 1.1])

    assert np.allclose(
        eos.bulk_modulus(volumes),
        -volumes * derivative(eos.pressure, volumes),
    )


def test_zero_pressure_derivative_limit():
    eos = Murnaghan(V0=10.0, K0=100.0, K0_prime=0.0)
    volumes = np.array([8.0, 10.0, 12.0])

    assert np.allclose(eos.pressure(volumes), eos.K0 * np.log(eos.V0 / volumes))
    assert np.all(eos.bulk_modulus(volumes) == eos.K0)
