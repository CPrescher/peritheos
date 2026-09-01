"""
Tests for the Holzapfel equation of state
"""

import numpy as np
import pytest

from peritheos.eos.rt.holzapfel import Holzapfel, bulk_modulus_derivative_analytical
from peritheos.utils import derivative

# Example values for diamond taken from Sokolova et al. 2016
V0 = 0.34141  # molar volume under standard conditions in Jbar^-1
K0 = 441.5  # isothermal bulk modulus in GPa
K0_prime = 3.9  # pressure derivative of bulk modulus
n = 1.0  # number of atoms in a chemical formula
Z = 6.0  # atomic number of the formula unit


@pytest.fixture
def holzapfel_eos():
    """Fixture that returns a Holzapfel EOS instance with test parameters"""
    return Holzapfel(V0=V0, K0=K0, K0_prime=K0_prime, n=n, Z=Z)


def test_init(holzapfel_eos):
    """Test initialization of Holzapfel EOS"""
    assert holzapfel_eos.V0 == V0
    assert holzapfel_eos.K0 == K0
    assert holzapfel_eos.K0_prime == K0_prime
    assert holzapfel_eos.n == n
    assert holzapfel_eos.Z == Z


def test_pressure_at_v0(holzapfel_eos):
    """Test pressure at reference volume is zero"""
    assert np.isclose(holzapfel_eos.pressure(V0), 0.0)


def test_pressure_compression(holzapfel_eos):
    """Test pressure under compression is positive"""
    assert holzapfel_eos.pressure(0.8 * V0) > 0


def test_pressure_expansion(holzapfel_eos):
    """Test pressure under expansion is negative"""
    assert holzapfel_eos.pressure(1.2 * V0) < 0


def test_pressure_array(holzapfel_eos):
    """Test pressure calculation with array input"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    pressures = holzapfel_eos.pressure(volumes)
    assert isinstance(pressures, np.ndarray)
    assert len(pressures) == len(volumes)
    assert np.all(pressures[:2] > 0)  # compression
    assert np.isclose(pressures[2], 0.0)  # reference volume
    assert np.all(pressures[3:] < 0)  # expansion


def test_bulk_modulus_at_v0(holzapfel_eos):
    """Test bulk modulus at reference volume equals K0"""
    assert np.isclose(holzapfel_eos.bulk_modulus(V0), K0)


def test_bulk_modulus_compression(holzapfel_eos):
    """Test bulk modulus increases under compression"""
    assert holzapfel_eos.bulk_modulus(0.8 * V0) > K0


def test_bulk_modulus_array(holzapfel_eos):
    """Test bulk modulus calculation with array input"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    moduli = holzapfel_eos.bulk_modulus(volumes)
    assert isinstance(moduli, np.ndarray)
    assert len(moduli) == len(volumes)
    assert np.isclose(moduli[2], K0)  # reference volume


def test_known_values_for_pressure(holzapfel_eos):
    """Test against pre-calculated values"""
    # These values were calculated using the original excel spreadsheet from Sokolova et al. 2016
    test_volumes = [
        0.307269,
        0.2901985,
        0.273128,
    ]
    test_pressures = [
        57.0207,
        98.0832,
        151.0189,
    ]
    for v, p in zip(test_volumes, test_pressures):
        assert np.isclose(holzapfel_eos.pressure(v), p)


def test_bulk_modulus_with_derivative(holzapfel_eos):
    """Test bulk modulus calculation with direct derivative as comparison"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    moduli = holzapfel_eos.bulk_modulus(volumes)
    moduli_from_derivative = -volumes * derivative(holzapfel_eos.pressure, volumes)

    assert np.allclose(moduli, moduli_from_derivative)


def test_bulk_modulus_derivative_at_V(holzapfel_eos):
    """Test the native derivative against an independent pressure derivative."""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    expected = derivative(holzapfel_eos.bulk_modulus, volumes) / derivative(
        holzapfel_eos.pressure, volumes
    )
    actual = holzapfel_eos.bulk_modulus_derivative(volumes)
    assert np.allclose(actual, expected, rtol=1e-5, atol=1e-7)
    assert np.isclose(holzapfel_eos.bulk_modulus_derivative(V0), K0_prime)


def test_legacy_analytical_derivative_helper_remains_compatible(holzapfel_eos):
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    fermigas_pressure = 1003.6 * (Z * n / (V0 * 10)) ** (5 / 3)
    c0 = -np.log(3 * K0 / fermigas_pressure)
    c2 = 1.5 * (K0_prime - 3) - c0
    actual = bulk_modulus_derivative_analytical(
        V0, volumes, holzapfel_eos.bulk_modulus(volumes), K0, c0, c2
    )
    assert np.allclose(actual, holzapfel_eos.bulk_modulus_derivative(volumes))
