"""
This module contains the thermal equations of state (EOS) implementations.
"""

from .double_debye_helmholtz import DoubleDebyeHelmholtz
from .holland_powell import HollandPowell2011, ThermalModifiedTait
from .linear import (
    LinearThermalPressure,
    LogVolumeThermalPressure,
    ThermalReferenceStateEOS,
)
from .mie_gruneisen import (
    MieGruneisenDebye,
    MieGruneisenEinstein,
    Tange2009Debye,
)
from .multi_oscillator import MultiOscillatorGruneisenThermalEOS
from .sokolova2016 import Sokolova2016

__all__ = [
    "DoubleDebyeHelmholtz",
    "HollandPowell2011",
    "LinearThermalPressure",
    "LogVolumeThermalPressure",
    "MieGruneisenDebye",
    "MieGruneisenEinstein",
    "MultiOscillatorGruneisenThermalEOS",
    "Sokolova2016",
    "Tange2009Debye",
    "ThermalModifiedTait",
    "ThermalReferenceStateEOS",
]
