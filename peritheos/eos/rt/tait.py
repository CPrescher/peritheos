"""Modified Tait equation of state."""

import numpy as np

from peritheos.eos import (
    EosBase,
    NumericType,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class ModifiedTait(EosBase):
    """Modified Tait isothermal equation of state.

    In terms of ``v = V / V0``, the pressure is

    ``P(V) = (((v + a - 1) / a)**(-1 / c) - 1) / b``,

    where ``a``, ``b``, and ``c`` are determined by ``K0``, ``K0_prime``,
    and ``K0_double_prime``. Pressure and both derivatives of the bulk
    modulus therefore have their supplied values at the reference state.

    When ``K0_double_prime`` is zero, this equation reduces exactly to the
    Murnaghan EOS.

    Reference
    ---------
    Holland, T. J. B. & Powell, R. (2011). An improved and extended
    internally consistent thermodynamic dataset for phases of petrological
    interest, involving a new equation of state for solids. Journal of
    Metamorphic Geology, 29, 333-383.
    """

    def __init__(
        self,
        V0: float,
        K0: float,
        K0_prime: float,
        K0_double_prime: float,
    ) -> None:
        """Initialize the modified Tait EOS.

        Parameters
        ----------
        V0 : float
            Reference volume in any consistent volume unit.
        K0 : float
            Bulk modulus at the reference volume.
        K0_prime : float
            First pressure derivative of the reference bulk modulus.
        K0_double_prime : float
            Second pressure derivative of the reference bulk modulus, in
            inverse pressure units.
        """
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self.K0_double_prime = validate_finite_scalar(
            K0_double_prime, "K0_double_prime"
        )

        one_plus_prime = 1.0 + self.K0_prime
        numerator_c = one_plus_prime + self.K0 * self.K0_double_prime
        denominator_c = (
            self.K0_prime**2
            + self.K0_prime
            - self.K0 * self.K0_double_prime
        )
        if one_plus_prime == 0.0 or numerator_c == 0.0 or denominator_c == 0.0:
            raise ValueError(
                "K0_prime and K0_double_prime produce a singular modified Tait EOS"
            )

        self._a = one_plus_prime / numerator_c
        self._b = (
            self.K0_prime / self.K0
            - self.K0_double_prime / one_plus_prime
        )
        self._c = numerator_c / denominator_c

    def _compression_base(self, V: NumericType) -> tuple[NumericType, NumericType]:
        relative_volume = V / self.V0
        base = (relative_volume + self._a - 1.0) / self._a
        if np.any(base <= 0.0):
            raise ValueError("Volume is outside the modified Tait EOS domain")
        return relative_volume, base

    def pressure(self, V: NumericType) -> NumericType:
        """Return pressure at volume *V* in the same units as ``K0``."""
        V = validate_volume(V)
        _, base = self._compression_base(V)
        return np.expm1((-1.0 / self._c) * np.log(base)) / self._b

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """Return the isothermal bulk modulus at volume *V*."""
        V = validate_volume(V)
        relative_volume, base = self._compression_base(V)
        return (
            self.K0
            * relative_volume
            * np.exp((-1.0 / self._c - 1.0) * np.log(base))
        )
