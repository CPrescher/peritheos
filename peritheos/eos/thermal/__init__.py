"""
This module contains the thermal equations of state (EOS) implementations.
"""

from .dewaele import Dewaele2006
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
    AsymptoticPowerLawMieGruneisenDebyeExcess,
    MieGruneisenDebye,
    MieGruneisenEinstein,
    Tange2009Debye,
)
from .multi_oscillator import MultiOscillatorGruneisenThermalEOS
from .sokolova2016 import Sokolova2016
from .sound_velocity_debye import SoundVelocityDebyeHelmholtz

__all__ = [
    "DoubleDebyeHelmholtz",
    "DoubleDebyeLogMomentHelmholtz",
    "Dewaele2006",
    "DorogokupetsOganov2007",
    "HollandPowell2011",
    "LinearThermalPressure",
    "LogVolumeThermalPressure",
    "MieGruneisenDebye",
    "AsymptoticPowerLawMieGruneisenDebyeExcess",
    "MieGruneisenEinstein",
    "MultiOscillatorGruneisenThermalEOS",
    "Sokolova2016",
    "SecondOrderTaylorThermalPressure",
    "SoundVelocityDebyeHelmholtz",
    "Tange2009Debye",
    "ThermalModifiedTait",
    "ThermalReferenceStateEOS",
]
