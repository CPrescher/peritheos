"""Tests for the modified Tait equation of state."""

import numpy as np
import pytest

from peritheos.eos.rt import ModifiedTait, Murnaghan
from peritheos.utils import derivative


@pytest.fixture
def eos():
    return ModifiedTait(
        V0=10.0,
        K0=100.0,
        K0_prime=4.0,
        K0_double_prime=-0.01,
    )


def test_parameters(eos):
    assert eos.V0 == 10.0
    assert eos.K0 == 100.0
    assert eos.K0_prime == 4.0
    assert eos.K0_double_prime == -0.01


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


def test_bulk_modulus_matches_pressure_derivative(eos):
    volumes = eos.V0 * np.array([0.8, 0.9, 1.0, 1.1])

    assert np.allclose(
        eos.bulk_modulus(volumes),
        -volumes * derivative(eos.pressure, volumes),
    )


def test_zero_second_derivative_reduces_to_murnaghan():
    tait = ModifiedTait(10.0, 100.0, 4.0, 0.0)
    murnaghan = Murnaghan(10.0, 100.0, 4.0)
    volumes = np.linspace(7.0, 12.0, 6)

    assert np.allclose(tait.pressure(volumes), murnaghan.pressure(volumes))
    assert np.allclose(tait.bulk_modulus(volumes), murnaghan.bulk_modulus(volumes))


@pytest.mark.parametrize(
    "K0_prime,K0_double_prime",
    [(-1.0, -0.01), (4.0, -0.05), (4.0, 0.2)],
)
def test_singular_parameter_sets_are_rejected(K0_prime, K0_double_prime):
    with pytest.raises(ValueError, match="singular modified Tait EOS"):
        ModifiedTait(10.0, 100.0, K0_prime, K0_double_prime)


def test_volume_outside_domain_is_rejected():
    eos = ModifiedTait(10.0, 100.0, 4.0, 0.01)

    with pytest.raises(ValueError, match="outside the modified Tait EOS domain"):
        eos.pressure(1.0)
