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
3. Anderson, O.L. (1995) Equations of State of Solids for Geophysics and Ceramic Science,
   Oxford University Press, Oxford, UK.
4. Poirier, J.P., (2000), Introduction to the Physics of the Earth's
   Interior, Cambridge University Press, Cambridge, UK.
"""

from peritheos.eos import (
    EosBase,
    NumericType,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


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

    def __init__(self, V0: float, K0: float) -> None:
        """
        Initialize the 2nd-order Birch-Murnaghan EOS.

        Parameters
        ----------
        V0 : float
            Reference volume (in cubic angstroms or any consistent unit)
        K0 : float
            Bulk modulus at reference volume (in GPa or any consistent unit)
        """
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")

    def pressure(self, V: NumericType) -> NumericType:
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
        V = validate_volume(V)
        f = ((self.V0 / V) ** (2 / 3) - 1) / 2
        return 3 * self.K0 * f * (1 + 2 * f) ** (5 / 2)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Return the bulk modulus at the given volume.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in the same units as V0)

        Returns
        -------
        float or numpy.ndarray
            Bulk modulus (in the same units as K0)
        """
        V = validate_volume(V)
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
    Anderson, O.L., (1995) Equations of State of Solids for Geophysics and Ceramic Science,
    Oxford University Press, Oxford, UK.
    """

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
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
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")

    def pressure(self, V: NumericType) -> NumericType:
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
        V = validate_volume(V)
        f = ((self.V0 / V) ** (2 / 3) - 1) / 2
        return (
            3
            * self.K0
            * f
            * (1 + 2 * f) ** (5 / 2)
            * (1 + (3 / 2) * (self.K0_prime - 4) * f)
        )

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Calculate the bulk modulus at the given volume.

        For the 3rd-order Birch-Murnaghan EOS, the bulk modulus at V0 is:
        K(V) = K0 * (1 + 2f)^(5/2) * [K0 + (3K0*K0'-5K0)f + 27/2(K0*K0'-4K0)f^2]

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in the same units as V0)

        Returns
        -------
        float or numpy.ndarray
            Bulk modulus (in the same units as K0)
        """
        V = validate_volume(V)
        f = ((self.V0 / V) ** (2 / 3) - 1) / 2
        return (1.0 + 2.0 * f) ** (5.0 / 2.0) * (
            self.K0
            + (3.0 * self.K0 * self.K0_prime - 5 * self.K0) * f
            + 27.0 / 2.0 * (self.K0 * self.K0_prime - 4.0 * self.K0) * f * f
        )

class BM4(EosBase):
    """
    4th-order Birch-Murnaghan equation of state.

    The 4th-order Birch-Murnaghan EOS is defined as:
    P(V) = (3/2) * K0 * [(V0/V)^(-7/3) - (V0/V)^(-5/3)] * 
            { 1 + (3/4) * (K0'- 4) * ((V0/V)^(-2/3) - 1) + (3/8) * [K0 * K0'' + (K0'- 4) * (K0' - 3) + (35/9)] * ((V0/V)^(-2/3) - 1)

    or in terms of the Eulerian strain f = [(V0/V)^(2/3) - 1]/2, zeta = (3/4) * (4 - K0'), and xi = (3/8) * [K0 * K0'' + (K0'- 4) * (K0' - 3) + (35/9)]:
    P(f) = 3 * K0 *f * (1 + 2f)^(5/2) * [1 + 2 * zeta * f + 4 * xi * f ** 2] 
    
    Equations from:
    Anderson, O.L., (1995) Equations of State of Solids for Geophysics and Ceramic Science,
    Oxford University Press, Oxford, UK.
    """

    def __init__(self, V0: float, K0: float, K0_prime: float, K0_double_prime: float) -> None:
        """
        Initialize the 4th-order Birch-Murnaghan EOS.

        Parameters
        ----------
        V0 : float
            Reference volume (in cubic angstroms or any consistent unit)
        K0 : float
            Bulk modulus at reference volume (in GPa or any consistent unit)
        K0_prime : float
            Pressure derivative of the bulk modulus at reference volume (dimensionless)
        K0_double_prime : float
            Second pressure derivative of the bulk modulus at reference volume
            (inverse pressure, e.g. GPa^-1 when K0 is in GPa)
        """
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self.K0_double_prime = validate_finite_scalar(
            K0_double_prime, "K0_double_prime"
        )

    def pressure(self, V: NumericType) -> NumericType:
        """
        Calculate pressure using the 4th-order Birch-Murnaghan EOS.

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in the same units as V0)

        Returns
        -------
        float or numpy.ndarray
            Pressure (in the same units as K0)
        """
        V = validate_volume(V)
        f = ((self.V0/V) ** (2 / 3) - 1) / 2
        zeta = (3 / 4) * (4 - self.K0_prime)
        xi = (3 / 8) * (self.K0 * self.K0_double_prime + (self.K0_prime - 4) * (self.K0_prime - 3) + (35 / 9))
        return (
            3 
            * self.K0 
            * f 
            * (1 + 2 * f) ** (5/2) 
            * (1 - 2 * zeta * f + 4 * xi * f ** 2)
        ) 

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Calculate the bulk modulus at the given volume.

        For the 4th-order Birch-Murnaghan EOS, the bulk modulus at V0 with f = f = ((V0/V) ** (2 / 3) - 1) / 2
        zeta = (3 / 4) * (4 - K0_prime)
        xi = (3 / 8) * [K0 * K0_double_prime + (K0_prime - 4) * (K0_prime - 3) + (35 / 9)]
        K(f) = 5 * self.K0 * (1 + 2 * f) ** (5 / 2) * (1 - 2 * zeta * f + 4 * xi * f ** 2) + self.K0 * (1 + 2 * f) ** (7 / 2)
            * (1 - 4 * zeta * f + 12 * xi * f ** 2)

        Parameters
        ----------
        V : float or numpy.ndarray
            Volume (in the same units as V0)

        Returns
        -------
        float or numpy.ndarray
            Bulk modulus (in the same units as K0)
        """
        V = validate_volume(V)
        f = ((self.V0/V) ** (2 / 3) - 1) / 2
        zeta = (3 / 4) * (4 - self.K0_prime)
        xi = (3 / 8) * (self.K0 * self.K0_double_prime + (self.K0_prime - 4) * (self.K0_prime - 3) + (35 / 9))
        return (
            5 
            * f
            * self.K0
            * (1 + 2 * f) ** (5 / 2)
            * (1 - 2 * zeta * f + 4 * xi * f ** 2) 
            + self.K0
            * (1 + 2 * f) ** (7 / 2)
            * (1 - 4 * zeta * f + 12 * xi * f ** 2)
        )
