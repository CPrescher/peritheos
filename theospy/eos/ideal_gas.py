"""
Ideal gas equation of state implementation
"""

import numpy as np
from theospy.constants import R


def calculate_pressure(temperature, volume, moles):
    """
    Calculate pressure using the ideal gas equation of state: P = nRT/V
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin
    volume : float
        Volume in cubic meters
    moles : float
        Number of moles
        
    Returns
    -------
    float
        Pressure in Pascal
    """
    return moles * R * temperature / volume


def calculate_volume(temperature, pressure, moles):
    """
    Calculate volume using the ideal gas equation of state: V = nRT/P
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin
    pressure : float
        Pressure in Pascal
    moles : float
        Number of moles
        
    Returns
    -------
    float
        Volume in cubic meters
    """
    return moles * R * temperature / pressure


def calculate_temperature(pressure, volume, moles):
    """
    Calculate temperature using the ideal gas equation of state: T = PV/nR
    
    Parameters
    ----------
    pressure : float
        Pressure in Pascal
    volume : float
        Volume in cubic meters
    moles : float
        Number of moles
        
    Returns
    -------
    float
        Temperature in Kelvin
    """
    return pressure * volume / (moles * R)


def calculate_moles(pressure, volume, temperature):
    """
    Calculate number of moles using the ideal gas equation of state: n = PV/RT
    
    Parameters
    ----------
    pressure : float
        Pressure in Pascal
    volume : float
        Volume in cubic meters
    temperature : float
        Temperature in Kelvin
        
    Returns
    -------
    float
        Number of moles
    """
    return pressure * volume / (R * temperature) 