"""
This module contains the room temperature equations of state (EOS) implementations.
"""

from .bm import BM2, BM3, BM4
from .holzapfel import Holzapfel
from .murnaghan import Murnaghan
from .natural_strain import NaturalStrain2, NaturalStrain3, NaturalStrain4
from .rydberg_stacey import RydbergStacey
from .sun_morse import Morse3, SunMorse3, SunMorse4
from .tait import ModifiedTait
from .vinet import Vinet

__all__ = [
    "BM2",
    "BM3",
    "BM4",
    "Holzapfel",
    "ModifiedTait",
    "Murnaghan",
    "NaturalStrain2",
    "NaturalStrain3",
    "NaturalStrain4",
    "RydbergStacey",
    "Morse3",
    "SunMorse3",
    "SunMorse4",
    "Vinet",
]
