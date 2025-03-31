import numpy as np
from peritheos.eos import EosBase, NumericType


class Holzapfel(EosBase):
    """
    Holzapfel equation of state. Introduced in Holzapfel et al. 2001. And here implemented as described in
    Sokolova et al. 2016.
    """

    def __init__(self, V0: float, K0: float, K0_prime: float, n: float, Z: int) -> None:
        """
        Initialize the Holzapfel equation of state.

        **Note:**
        Contrary to other EOS implementations, the reference values need to be provided in the correct units.

        Parameters
        ----------
        V0 : float
            Reference volume [JBar^-1] (same as [cm^3/mol]/10)
        K0 : float
            Bulk modulus at reference volume in [kbar]
        K0_prime : float
            Pressure derivative of bulk modulus at reference volume
        n : float
            Number of atoms in a chemical formula
        Z : int
            Atomic number of the forumala unit
        """
        super().__init__()
        self.V0: float = V0
        self.K0: float = K0
        self.K0_prime: float = K0_prime
        self.n: float = n
        self.Z: float = Z

        self._P_FG0 = 1003.6 * (self.Z * self.n / (self.V0 * 10)) ** (5 / 3)
        self._c0 = -np.log(3 * self.K0 / 10 / self._P_FG0)
        self._c2 = 1.5 * (self.K0_prime - 3) - self._c0

    def pressure(self, V: NumericType) -> NumericType:
        """
        Calculate pressure using the Holzapfel equation of state.

        Parameters
        ----------
        V : float or np.ndarray
            Volume in [JBar^-1] (same as [cm^3/mol]/10)

        Returns
        -------
        float or np.ndarray
            Pressure in bar
        """
        x = (V / self.V0) ** (1 / 3)
        P_value = (
            3
            * self.K0
            * 1000
            * np.exp(self._c0 * (1 - x))
            * (1 / x**5 - 1 / x**4)
            * (1 + self._c2 * x - self._c2 * x**2)
        )
        return P_value

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """
        Calculate bulk modulus using the Holzapfel equation of state.

        Parameters
        ----------
        V : float or np.ndarray
            Volume in [JBar^-1] (same as [cm^3/mol]/10)

        Returns
        -------
        float or np.ndarray
            Bulk modulus in [kbar]
        """
        x = (V / self.V0) ** (1 / 3)
        term1 = self.K0 * 1000 * x**-5 * np.exp(self._c0 * (1 - x))

        bracket_1 = (5 - 4 * x) * (1 + self._c2 * x * (1 - x))
        bracket_2 = self._c0 * x * (1 - x) * (1 + self._c2 * x * (1 - x))
        bracket_3 = -(1 - x) * (self._c2 * x - 2 * self._c2 * x**2)
        term2 = bracket_1 + bracket_2 + bracket_3
        return term1 * term2 / 1000

    def bulk_modulus_derivative(self, V: NumericType, eps: float = 1e-6) -> NumericType:
        """
        Compute the second derivative of the bulk modulus with respect to pressure using a
        numerical approximation.

        This method provides a faster alternative to the analytical approach implemented in
        `bulk_modulus_derivative_analytical`, which is more than twice as slow. The numerical
        derivative is computed using a small relative perturbation in volume with a step size of `eps`.

        Parameters
        ----------
        V : float or np.ndarray
            Volume in [JBar^-1] (same as [cm^3/mol]/10)
        eps : float
            Step size for the numerical derivative, default is 1e-6. This value is used as the
            relative difference between the two perturbed volumes. The accuracy of the numerical
            derivative is approximately `eps * 100%`.

        Returns
        -------
        float or np.ndarray
            Second derivative of bulk modulus with respect to pressure (unitless)
        """
        V = np.array([V * (1 + eps), V * (1 - eps)])
        P = self.pressure(V)
        K = self.bulk_modulus(V) * 1000
        return (K[0] - K[1]) / (P[0] - P[1])


def bulk_modulus_derivative_analytical(V0, V, KT, K0, c0, c2):
    """
    Calculate the second derivative of bulk modulus with respect to pressure using the Holzapfel equation of state.
    Equation have been taken from excel sheets provided in Sokolova et al. 2016. The listed equation (5) on page
    163 in this paper seems to be incorrect.
    Values have been verified against the excel sheets provided by Sokolova et al. 2016 and the numerical
    derivative (see test_holzapfel.py).

    Parameters
    ----------
    V : float or np.ndarray
        Volume in [JBar^-1] (same as [cm^3/mol]/10)
    KT : float or np.ndarray
        Bulk modulus in [kbar]
    K0 : float
        Bulk modulus at reference volume in [kbar]
    c0 : float
        Constant in the Holzapfel equation of state
    c2 : float
        Constant in the Holzapfel equation of state
    Returns
    -------
    """
    x = (V / V0) ** (1 / 3)

    # convert moduli from kbar to bar
    K0 = K0 * 1000
    KT = KT * 1000

    term1 = (
        3
        / x**4
        * K0
        * np.exp(c0 * (1 - x))
        * (
            (-5 / x**2 + 4 / x) * (1 + c2 * x - c2 * x**2)
            - (1 / x - 1) * (1 + c2 * x - c2 * x**2) * c0
            + (1 / x - 1) * (c2 - 2 * c2 * x)
        )
    )
    term2 = (
        1
        / x**3
        * K0
        * np.exp(c0 * (1 - x))
        * c0
        * (
            (-5 / x**2 + 4 / x) * (1 + c2 * x - c2 * x**2)
            - (1 / x - 1) * (1 + c2 * x - c2 * x**2) * c0
            + (1 / x - 1) * (c2 - 2 * c2 * x)
        )
    )
    term3 = (
        K0
        * np.exp(c0 * (1 - x))
        * (
            (10 / x**3 - 4 / x**2) * (1 + c2 * x - c2 * x**2)
            + (-5 / x**2 + 4 / x) * (c2 - 2 * c2 * x)
            + c0 / x**2 * (1 + c2 * x - c2 * x**2)
            - (1 / x - 1) * (c2 - 2 * c2 * x) * c0
            - (c2 - 2 * c2 * x) / x**2
            - 2 * c2 * (1 / x - 1)
        )
        / x**3
    )
    return (term1 + term2 - term3) / (-KT / x) / 3
