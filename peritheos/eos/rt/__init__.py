"""
This module contains the room temperature equations of state (EOS) implementations.
"""

from .baonza import Baonza
from .bm import BM2, BM3, BM4
from .density_polynomial import DensityPolynomial3
from .holzapfel import Holzapfel
from .murnaghan import Murnaghan
from .natural_strain import NaturalStrain2, NaturalStrain3, NaturalStrain4
from .odd_inverse_power import OddInversePower
from .rydberg_stacey import RydbergStacey
from .second_order_murnaghan import SecondOrderMurnaghan
from .sun_morse import Morse3, SunMorse3, SunMorse4
from .tait import ModifiedTait
from .vinet import Vinet
from .vinet3 import Vinet3

__all__ = [
    "BM2",
    "BM3",
    "BM4",
    "Baonza",
    "DensityPolynomial3",
    "Holzapfel",
    "ModifiedTait",
    "Murnaghan",
    "SecondOrderMurnaghan",
    "OddInversePower",
    "NaturalStrain2",
    "NaturalStrain3",
    "NaturalStrain4",
    "RydbergStacey",
    "Morse3",
    "SunMorse3",
    "SunMorse4",
    "Vinet",
    "Vinet3",
]
