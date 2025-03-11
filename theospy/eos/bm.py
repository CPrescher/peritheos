"""
Non-thermal equations of state module

This module implements various non-thermal equations of state commonly used in
high-pressure physics and materials science, including:

1. Birch-Murnaghan equation of state (2nd, 3rd, and 4th order)
2. Object-oriented implementation with a base class and specific EOS implementations

These equations describe the relationship between pressure and volume (or density)
at constant temperature, and are particularly useful for modeling the behavior of
solids under compression.

References:
-----------
1. Birch, F. (1947). Finite elastic strain of cubic crystals. Physical Review, 71(11), 809.
2. Birch, F. (1978). Finite strain isotherm and velocities for single-crystal and
   polycrystalline NaCl at high pressures and 300 K. Journal of Geophysical Research,
   83(B3), 1257-1268.
"""

import numpy as np


class EosBase:
    """
    Base class for equation of state implementations.

    This abstract class defines the interface that all equation of state
    implementations should follow.
    """

    def pressure(self, V):
        """
        Calculate pressure at a given volume.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in cubic angstroms or any consistent unit)

        Returns
        -------
        float or numpy.ndarray
            Pressure (in the same units as K0)
        """
        raise NotImplementedError("Subclasses must implement the pressure method.")

    def bulk_modulus(self):
        """
        Calculate the bulk modulus.

        Returns
        -------
        float
            Bulk modulus (in the same units as K0)
        """
        raise NotImplementedError("Subclasses must implement the bulk_modulus method.")


class BM2(EosBase):
    """
    2nd-order Birch-Murnaghan equation of state.

    The 2nd-order Birch-Murnaghan EOS is defined as:
    P(V) = (3/2) * K0 * [(V0/V)^(7/3) - (V0/V)^(5/3)]

    or in terms of the Eulerian strain f = [(V0/V)^(2/3) - 1]/2:
    P(V) = 3 * K0 * f * (1 + 2f)^(5/2)

    Equations from:
    Poirier, J.P., 2000, Introduction to the Physics of the Earth's
    Interior, Cambridge University Press, Cambridge, UK.
    """

    def __init__(self, V0, K0):
        """
        Initialize the 2nd-order Birch-Murnaghan EOS.

        Parameters
        ----------
        V0 : float
            Reference volume (in cubic angstroms or any consistent unit)
        K0 : float
            Bulk modulus at reference volume (in GPa or any consistent unit)
        """
        self.V0 = V0
        self.K0 = K0

    def pressure(self, V):
        """
        Calculate pressure using the 2nd-order Birch-Murnaghan EOS.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in the same units as V0)

        Returns
        -------
        float or numpy.ndarray
            Pressure (in the same units as K0)
        """
        f = ((self.V0 / V) ** (2 / 3) - 1) / 2
        return 3 * self.K0 * f * (1 + 2 * f) ** (5 / 2)

    def bulk_modulus(self, V):
        """
        Return the bulk modulus at the reference volume.

        Returns
        -------
        float
            Bulk modulus (in the same units as K0)
        """
        eta = ((self.V0 / V) ** (2 / 3) - 1) / 2
        return self.K0 * (1 + 7 * eta) * (1 + 2 * eta) ** (5 / 2)


class BM3(EosBase):
    """
    3rd-order Birch-Murnaghan equation of state.

    The 3rd-order Birch-Murnaghan EOS is defined as:
    P(V) = (3/2) * K0 * [(V0/V)^(7/3) - (V0/V)^(5/3)] *
           {1 + (3/4) * (K0' - 4) * [(V0/V)^(2/3) - 1]}

    or in terms of the Eulerian strain f = [(V0/V)^(2/3) - 1]/2:
    P(V) = 3 * K0 * f * (1 + 2f)^(5/2) * [1 + (3/2) * (K0' - 4) * f]

    Equations from:
    Anderson, O.L., Equations of State of Solids for Geophysics and Ceramic Science,
    Oxford University Press, Oxford, UK, 2000.
    """

    def __init__(self, V0, K0, K0_prime):
        """
        Initialize the 3rd-order Birch-Murnaghan EOS.

        Parameters
        ----------
        V0 : float
            Reference volume (in cubic angstroms or any consistent unit)
        K0 : float
            Bulk modulus at reference volume (in GPa or any consistent unit)
        K0_prime : float
            Pressure derivative of the bulk modulus at reference volume (dimensionless)
        """
        self.V0 = V0
        self.K0 = K0
        self.K0_prime = K0_prime

    def pressure(self, V):
        """
        Calculate pressure using the 3rd-order Birch-Murnaghan EOS.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in the same units as V0)

        Returns
        -------
        float or numpy.ndarray
            Pressure (in the same units as K0)
        """
        f = ((self.V0 / V) ** (2 / 3) - 1) / 2
        return (
            3
            * self.K0
            * f
            * (1 + 2 * f) ** (5 / 2)
            * (1 + (3 / 2) * (self.K0_prime - 4) * f)
        )

    def bulk_modulus(self, V):
        """
        Calculate the bulk modulus at the reference volume.

        For the 3rd-order Birch-Murnaghan EOS, the bulk modulus at V0 is:


        Returns
        -------
        float
            Bulk modulus (in the same units as K0)
        """
        f = ((self.V0 / V) ** (2 / 3) - 1) / 2
        return (1.0 + 2.0 * f) ** (5.0 / 2.0) * (
            self.K0
            + (3.0 * self.K0 * self.K0_prime - 5 * self.K0) * f
            + 27.0 / 2.0 * (self.K0 * self.K0_prime - 4.0 * self.K0) * f * f
        )
