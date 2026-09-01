"""Compatibility imports for the former paper-named thermal EOS module."""

from .multi_oscillator import (
    I_gamV,
    MultiOscillatorGruneisenThermalEOS,
    f_gamV,
)

Sokolova2016 = MultiOscillatorGruneisenThermalEOS

__all__ = [
    "I_gamV",
    "MultiOscillatorGruneisenThermalEOS",
    "Sokolova2016",
    "f_gamV",
]
