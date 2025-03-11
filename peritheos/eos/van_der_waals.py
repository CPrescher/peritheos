"""
Van der Waals equation of state implementation
"""

import numpy as np
from peritheos.constants import R


def calculate_pressure(temperature, volume, moles, a, b):
    """
    Calculate pressure using the Van der Waals equation of state:
    P = nRT/(V-nb) - a(n/V)²
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin
    volume : float
        Volume in cubic meters
    moles : float
        Number of moles
    a : float
        Attraction parameter in Pa·m^6/mol^2
    b : float
        Volume parameter in m^3/mol
        
    Returns
    -------
    float
        Pressure in Pascal
    """
    return (moles * R * temperature) / (volume - moles * b) - a * (moles / volume) ** 2


def calculate_volume_numerical(temperature, pressure, moles, a, b, initial_guess=None, max_iterations=100, tolerance=1e-6):
    """
    Calculate volume using the Van der Waals equation of state numerically
    using Newton-Raphson method
    
    Parameters
    ----------
    temperature : float
        Temperature in Kelvin
    pressure : float
        Pressure in Pascal
    moles : float
        Number of moles
    a : float
        Attraction parameter in Pa·m^6/mol^2
    b : float
        Volume parameter in m^3/mol
    initial_guess : float, optional
        Initial guess for volume in cubic meters
    max_iterations : int, optional
        Maximum number of iterations for numerical solution
    tolerance : float, optional
        Convergence tolerance
        
    Returns
    -------
    float
        Volume in cubic meters
    """
    # Initial guess (use ideal gas as starting point if not provided)
    if initial_guess is None:
        initial_guess = moles * R * temperature / pressure
    
    v = initial_guess
    
    for i in range(max_iterations):
        # Van der Waals equation in the form f(V) = 0
        f_v = pressure - (moles * R * temperature) / (v - moles * b) + a * (moles / v) ** 2
        
        # Derivative of f(V)
        df_v = (moles * R * temperature) / ((v - moles * b) ** 2) - 2 * a * (moles ** 2) / (v ** 3)
        
        # Newton-Raphson update
        v_new = v - f_v / df_v
        
        # Check for convergence
        if abs(v_new - v) < tolerance:
            return v_new
        
        v = v_new
    
    # If we reach here, we didn't converge
    raise ValueError(f"Failed to converge after {max_iterations} iterations")


def critical_constants(a, b):
    """
    Calculate critical constants for a substance using Van der Waals parameters
    
    Parameters
    ----------
    a : float
        Attraction parameter in Pa·m^6/mol^2
    b : float
        Volume parameter in m^3/mol
        
    Returns
    -------
    dict
        Dictionary containing critical temperature, pressure, and volume
    """
    # Critical temperature: Tc = 8a/27Rb
    t_c = 8 * a / (27 * R * b)
    
    # Critical pressure: Pc = a/27b²
    p_c = a / (27 * b ** 2)
    
    # Critical volume: Vc = 3nb
    v_c = 3 * b
    
    return {
        "temperature": t_c,
        "pressure": p_c,
        "volume": v_c
    } 