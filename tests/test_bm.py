"""
Tests for the Birch-Murnaghan equation of state implementations.

This module contains tests for the 2nd, 3rd and 4th order Birch-Murnaghan
equations of state, verifying both pressure and bulk modulus calculations.
"""

import numpy as np
import pytest

from peritheos.eos.rt.bm import BM2, BM3, BM4
from peritheos.utils import derivative

V0 = 100
K0 = 100
K0_prime = 5
K0_double_prime = 10


@pytest.fixture
def bm2_eos():
    return BM2(V0, K0)


def test_bm2_init(bm2_eos):
    """Test initialization of BM2 EOS"""
    assert bm2_eos.V0 == V0
    assert bm2_eos.K0 == K0


def test_bm2_pressure_at_v0(bm2_eos):
    """Test pressure at reference volume is zero"""
    assert np.isclose(bm2_eos.pressure(V0), 0.0)


def test_bm2_pressure_compression(bm2_eos):
    """Test pressure under compression is positive"""
    assert bm2_eos.pressure(0.8 * V0) > 0


def test_bm2_pressure_expansion(bm2_eos):
    """Test pressure under expansion is negative"""
    assert bm2_eos.pressure(1.2 * V0) < 0


def test_bm2_pressure_array(bm2_eos):
    """Test pressure calculation with array input"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    pressures = bm2_eos.pressure(volumes)
    assert isinstance(pressures, np.ndarray)


def test_bm2_bulk_modulus_at_v0(bm2_eos):
    """Test bulk modulus at reference volume equals K0"""
    assert np.isclose(bm2_eos.bulk_modulus(V0), K0)


def test_bm2_bulk_modulus_compression(bm2_eos):
    """Test bulk modulus under compression is greater than K0"""
    assert bm2_eos.bulk_modulus(0.8 * V0) > K0


def test_bm2_bulk_modulus_expansion(bm2_eos):
    """Test bulk modulus under expansion is less than K0"""
    assert bm2_eos.bulk_modulus(1.2 * V0) < K0


def test_bm2_K_with_derivative(bm2_eos):
    """Test bulk modulus calculation from equation of state with direct derivative as comparison"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    moduli = bm2_eos.bulk_modulus(volumes)
    moduli_from_derivative = -volumes * derivative(bm2_eos.pressure, volumes)
    assert np.allclose(moduli, moduli_from_derivative)


def test_bm2_calculate_volume(bm2_eos):
    """Test volume calculation at a given pressure"""
    pressure = 100
    volume = bm2_eos.calculate_volume(pressure)
    pressure_2 = bm2_eos.pressure(volume)

    assert np.isclose(pressure, pressure_2)


def test_bm2_calculate_volume_array(bm2_eos):
    """Test volume calculation at a given pressure"""
    pressures = [100, 200, 300, 400, 500]
    volumes = bm2_eos.calculate_volume(pressures)
    pressures_2 = bm2_eos.pressure(volumes)
    assert np.allclose(pressures, pressures_2)


@pytest.fixture
def bm3_eos():
    return BM3(V0, K0, K0_prime)


def test_bm3_init(bm3_eos):
    """Test initialization of BM3 EOS"""
    assert bm3_eos.V0 == V0
    assert bm3_eos.K0 == K0
    assert bm3_eos.K0_prime == K0_prime


def test_bm3_pressure_at_v0(bm3_eos):
    """Test pressure at reference volume is zero"""
    assert np.isclose(bm3_eos.pressure(V0), 0.0)


def test_bm3_pressure_compression(bm3_eos):
    """Test pressure under compression is positive"""
    assert bm3_eos.pressure(0.8 * V0) > 0


def test_bm3_pressure_expansion(bm3_eos):
    """Test pressure under expansion is negative"""
    assert bm3_eos.pressure(1.2 * V0) < 0


def test_bm3_pressure_array(bm3_eos):
    """Test pressure calculation with array input"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    pressures = bm3_eos.pressure(volumes)
    assert isinstance(pressures, np.ndarray)


def test_bm3_bulk_modulus_at_v0(bm3_eos):
    """Test bulk modulus at reference volume equals K0"""
    assert np.isclose(bm3_eos.bulk_modulus(V0), K0)


def test_bm3_bulk_modulus_compression(bm3_eos):
    """Test bulk modulus under compression is greater than K0"""
    assert bm3_eos.bulk_modulus(0.8 * V0) > K0


def test_bm3_K_with_derivative(bm3_eos):
    """Test bulk modulus calculation from equation of state with direct derivative as comparison"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    moduli = bm3_eos.bulk_modulus(volumes)
    moduli_from_derivative = -volumes * derivative(bm3_eos.pressure, volumes)
    assert np.allclose(moduli, moduli_from_derivative)


@pytest.fixture
def bm4_eos():
    return BM4(V0, K0, K0_prime, K0_double_prime)


def test_bm4_init(bm4_eos):
    """Test initialization of BM4 EOS"""
    assert bm4_eos.V0 == V0
    assert bm4_eos.K0 == K0
    assert bm4_eos.K0_prime == K0_prime
    assert bm4_eos.K0_double_prime == K0_double_prime


def test_bm4_pressure_at_v0(bm4_eos):
    """Test pressure at reference volume is zero"""
    assert np.isclose(bm4_eos.pressure(V0), 0.0)


def test_bm4_pressure_compression(bm4_eos):
    """Test pressure under compression is positive"""
    assert bm4_eos.pressure(0.8 * V0) > 0


def test_bm4_pressure_expansion(bm4_eos):
    """Test pressure under expansion is negative"""
    assert bm4_eos.pressure(1.2 * V0) < 0


def test_bm4_pressure_array(bm4_eos):
    """Test pressure calculation with array input"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    pressures = bm4_eos.pressure(volumes)
    assert isinstance(pressures, np.ndarray)


def test_bm4_bulk_modulus_at_v0(bm4_eos):
    """Test bulk modulus at reference volume equals K0"""
    assert np.isclose(bm4_eos.bulk_modulus(V0), K0)


def test_bm4_bulk_modulus_compression(bm4_eos):
    """Test bulk modulus under compression is greater than K0"""
    assert bm4_eos.bulk_modulus(0.8 * V0) > K0


def test_bm4_K_with_derivative(bm4_eos):
    """Test bulk modulus calculation from equation of state with direct derivative as comparison"""
    volumes = np.array([0.8, 0.9, 1.0, 1.1, 1.2]) * V0
    moduli = bm4_eos.bulk_modulus(volumes)
    moduli_from_derivative = -volumes * derivative(bm4_eos.pressure, volumes)
    assert np.allclose(moduli, moduli_from_derivative)
