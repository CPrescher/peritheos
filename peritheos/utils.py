"""
Utility functions for thermodynamic calculations
"""

from .constants import R, N_A, k_B


def convert_pressure(value, from_unit, to_unit):
    """
    Convert pressure between different units
    
    Parameters
    ----------
    value : float
        Pressure value to convert
    from_unit : str
        Original unit ('pa', 'bar', 'atm', 'torr', 'psi')
    to_unit : str
        Target unit ('pa', 'bar', 'atm', 'torr', 'psi')
        
    Returns
    -------
    float
        Converted pressure value
    """
    # Conversion factors to Pa
    to_pa = {
        'pa': 1.0,
        'bar': 1e5,
        'atm': 101325.0,
        'torr': 133.322,
        'psi': 6894.76
    }
    
    # Convert to Pa first
    pa_value = value * to_pa[from_unit.lower()]
    
    # Convert from Pa to target unit
    return pa_value / to_pa[to_unit.lower()]


def convert_temperature(value, from_unit, to_unit):
    """
    Convert temperature between different units
    
    Parameters
    ----------
    value : float
        Temperature value to convert
    from_unit : str
        Original unit ('k', 'c', 'f')
    to_unit : str
        Target unit ('k', 'c', 'f')
        
    Returns
    -------
    float
        Converted temperature value
    """
    from_unit = from_unit.lower()
    to_unit = to_unit.lower()
    
    # Convert to Kelvin first
    if from_unit == 'k':
        kelvin = value
    elif from_unit == 'c':
        kelvin = value + 273.15
    elif from_unit == 'f':
        kelvin = (value - 32) * 5/9 + 273.15
    else:
        raise ValueError(f"Unsupported temperature unit: {from_unit}")
    
    # Convert from Kelvin to target unit
    if to_unit == 'k':
        return kelvin
    elif to_unit == 'c':
        return kelvin - 273.15
    elif to_unit == 'f':
        return (kelvin - 273.15) * 9/5 + 32
    else:
        raise ValueError(f"Unsupported temperature unit: {to_unit}")


def compressibility_factor(pressure, volume, temperature, moles):
    """
    Calculate the compressibility factor Z = PV/nRT
    
    Parameters
    ----------
    pressure : float
        Pressure in Pascal
    volume : float
        Volume in cubic meters
    temperature : float
        Temperature in Kelvin
    moles : float
        Number of moles
        
    Returns
    -------
    float
        Compressibility factor (dimensionless)
    """
    return pressure * volume / (moles * R * temperature) 

def derivative(f, x, dx=1e-6):
    """Compute the derivative of f at x using finite differences"""
    return (f(x + dx) - f(x - dx)) / (2 * dx)
