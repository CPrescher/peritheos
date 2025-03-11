"""
Tests for the ideal gas equation of state implementation
"""

import pytest
import numpy as np
from peritheos.eos import ideal_gas
from peritheos.constants import R


def test_calculate_pressure():
    """Test the ideal gas pressure calculation"""
    # Test with standard values
    temperature = 298.15  # K
    volume = 0.0224  # m³
    moles = 1.0  # mol
    
    expected_pressure = moles * R * temperature / volume
    calculated_pressure = ideal_gas.calculate_pressure(temperature, volume, moles)
    
    assert np.isclose(calculated_pressure, expected_pressure)
    assert np.isclose(calculated_pressure, 110667.7, rtol=1e-5)


def test_calculate_volume():
    """Test the ideal gas volume calculation"""
    # Test with standard values
    temperature = 298.15  # K
    pressure = 101325.0  # Pa
    moles = 1.0  # mol
    
    expected_volume = moles * R * temperature / pressure
    calculated_volume = ideal_gas.calculate_volume(temperature, pressure, moles)
    
    assert np.isclose(calculated_volume, expected_volume)
    assert np.isclose(calculated_volume, 0.0245, atol=1e-3)


def test_calculate_temperature():
    """Test the ideal gas temperature calculation"""
    # Test with standard values
    pressure = 101325.0  # Pa
    volume = 0.0224  # m³
    moles = 1.0  # mol
    
    expected_temperature = pressure * volume / (moles * R)
    calculated_temperature = ideal_gas.calculate_temperature(pressure, volume, moles)
    
    assert np.isclose(calculated_temperature, expected_temperature)
    assert np.isclose(calculated_temperature, 273.15, rtol=1e-2)


def test_calculate_moles():
    """Test the ideal gas moles calculation"""
    # Test with standard values
    pressure = 101325.0  # Pa
    volume = 0.0224  # m³
    temperature = 273.15  # K
    
    expected_moles = pressure * volume / (R * temperature)
    calculated_moles = ideal_gas.calculate_moles(pressure, volume, temperature)
    
    assert np.isclose(calculated_moles, expected_moles)
    assert np.isclose(calculated_moles, 1.0, rtol=1e-2)


def test_ideal_gas_law_consistency():
    """Test that all functions are consistent with the ideal gas law"""
    # Define a set of values
    temperature = 300.0  # K
    pressure = 200000.0  # Pa
    volume = 0.012  # m³
    moles = 0.8  # mol
    
    # Calculate pressure and check consistency
    calc_pressure = ideal_gas.calculate_pressure(temperature, volume, moles)
    assert np.isclose(calc_pressure * volume, moles * R * temperature)
    
    # Calculate volume and check consistency
    calc_volume = ideal_gas.calculate_volume(temperature, pressure, moles)
    assert np.isclose(pressure * calc_volume, moles * R * temperature)
    
    # Calculate temperature and check consistency
    calc_temperature = ideal_gas.calculate_temperature(pressure, volume, moles)
    assert np.isclose(pressure * volume, moles * R * calc_temperature)
    
    # Calculate moles and check consistency
    calc_moles = ideal_gas.calculate_moles(pressure, volume, temperature)
    assert np.isclose(pressure * volume, calc_moles * R * temperature) 