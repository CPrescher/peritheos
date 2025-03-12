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
        V : float
            Volume in [JBar^-1] (same as [cm^3/mol]/10)

        Returns
        -------
        float
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
        V : float
            Volume in [JBar^-1] (same as [cm^3/mol]/10)

        Returns
        -------
        float
            Bulk modulus in [kbar]
        """
        x = (V / self.V0) ** (1 / 3)
        term1 = self.K0 * x**-5 * np.exp(self._c0 * (1 - x))

        bracket_1 = (5 - 4 * x) * (1 + self._c2 * x * (1 - x))
        bracket_2 = self._c0 * x * (1 - x) * (1 + self._c2 * x * (1 - x))
        bracket_3 = -(1 - x) * (self._c2 * x - 2 * self._c2 * x**2)
        term2 = bracket_1 + bracket_2 + bracket_3
        return term1 * term2

