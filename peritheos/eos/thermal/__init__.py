"""
This module contains the thermal equations of state (EOS) implementations.
"""

from .dorogokupets_oganov import DorogokupetsOganov2007
from .double_debye_helmholtz import (
    DoubleDebyeHelmholtz,
    DoubleDebyeLogMomentHelmholtz,
)
from .holland_powell import HollandPowell2011, ThermalModifiedTait
from .linear import (
    LinearThermalPressure,
    LogVolumeThermalPressure,
    SecondOrderTaylorThermalPressure,
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
    "DoubleDebyeLogMomentHelmholtz",
    "DorogokupetsOganov2007",
    "HollandPowell2011",
    "LinearThermalPressure",
    "LogVolumeThermalPressure",
    "MieGruneisenDebye",
    "MieGruneisenEinstein",
    "MultiOscillatorGruneisenThermalEOS",
    "Sokolova2016",
    "SecondOrderTaylorThermalPressure",
    "Tange2009Debye",
    "ThermalModifiedTait",
    "ThermalReferenceStateEOS",
]
