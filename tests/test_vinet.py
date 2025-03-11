"""
Tests for the Vinet equation of state
"""

import numpy as np
import pytest
from peritheos.eos.vinet import Vinet
from peritheos.utils import derivative

# Example values for a material
V0 = 10.0  # reference volume in cubic angstroms
K0 = 100.0  # bulk modulus in GPa
K0_prime = 4.0  # pressure derivative of bulk modulus


@pytest.fixture
def vinet_eos():
    """Fixture that returns a Vinet EOS instance with test parameters"""
    return Vinet(V0=V0, K0=K0, K0_prime=K0_prime)


def test_init(vinet_eos):
    """Test initialization of Vinet EOS"""
    assert vinet_eos.V0 == V0
    assert vinet_eos.K0 == K0
    assert vinet_eos.K0_prime == K0_prime


def test_pressure_at_v0(vinet_eos):
    """Test pressure at reference volume is zero"""
    assert np.isclose(vinet_eos.pressure(V0), 0.0)


def test_pressure_compression(vinet_eos):
    """Test pressure under compression is positive"""
    assert vinet_eos.pressure(0.8 * V0) > 0


def test_pressure_expansion(vinet_eos):
    """Test pressure under expansion is negative"""
    assert vinet_eos.pressure(1.2 * V0) < 0


def test_pressure_array(vinet_eos):
    """Test pressure calculation with array input"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    pressures = vinet_eos.pressure(volumes)
    assert isinstance(pressures, np.ndarray)
    assert len(pressures) == len(volumes)
    assert np.all(pressures[:2] > 0)  # compression
    assert np.isclose(pressures[2], 0.0)  # reference volume
    assert np.all(pressures[3:] < 0)  # expansion


def test_bulk_modulus_at_v0(vinet_eos):
    """Test bulk modulus at reference volume equals K0"""
    assert np.isclose(vinet_eos.bulk_modulus(V0), K0)


def test_bulk_modulus_compression(vinet_eos):
    """Test bulk modulus increases under compression"""
    assert vinet_eos.bulk_modulus(0.8 * V0) > K0


def test_bulk_modulus_array(vinet_eos):
    """Test bulk modulus calculation with array input"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    moduli = vinet_eos.bulk_modulus(volumes)
    assert isinstance(moduli, np.ndarray)
    assert len(moduli) == len(volumes)
    assert np.isclose(moduli[2], K0)  # reference volume


def test_known_values(vinet_eos):
    """Test against pre-calculated values"""
    # These values were calculated using the equations
    test_v = 0.9 * V0
    f = (test_v / V0) ** (1 / 3)
    expected_p = 3 * K0 * (1 - f) / f**2 * np.exp(3 / 2 * (K0_prime - 1) * (1 - f))
    expected_k = (
        K0
        * f**-2
        * (1 + (3 / 2 * (K0_prime - 1) * f + 1) * (1 - f))
        * np.exp(3 / 2 * (K0_prime - 1) * (1 - f))
    )

    assert np.isclose(vinet_eos.pressure(test_v), expected_p)
    assert np.isclose(vinet_eos.bulk_modulus(test_v), expected_k)


def test_compare_bulk_modulus_with_derivative(vinet_eos):
    """Test bulk modulus calculation from equation of state with direct derivative as comparison"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    moduli = vinet_eos.bulk_modulus(volumes)
    moduli_from_derivative = -volumes * derivative(vinet_eos.pressure, volumes)
    assert np.allclose(moduli, moduli_from_derivative)


