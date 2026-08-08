"""
Physical constants used in thermodynamic calculations
"""

from scipy import constants

# Universal gas constant (J/(mol·K))
R = constants.R

# Boltzmann constant (J/K)
k_B = constants.Boltzmann

# Avogadro's number (1/mol)
N_A = constants.Avogadro

# Standard temperature and pressure
STP_TEMPERATURE = constants.zero_Celsius  # K (273.15)
STP_PRESSURE = constants.atm  # Pa (101325.0)

# Other useful constants
STANDARD_GRAVITY = constants.g  # m/s²
STEFAN_BOLTZMANN = constants.Stefan_Boltzmann  # W/(m²·K⁴)
