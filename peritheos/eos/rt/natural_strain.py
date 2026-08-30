"""Natural-strain (Poirier-Tarantola) equations of state."""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class _NaturalStrainBase(EosBase):
    def __init__(self, V0: float, K0: float) -> None:
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")

    def _coefficients(self) -> tuple[float, float]:
        raise NotImplementedError

    def pressure(self, V: NumericType) -> NumericType:
        """Return pressure at volume *V* in the same units as ``K0``."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "pressure", V)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """Return isothermal bulk modulus at volume *V*."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "bulk_modulus", V)


class NaturalStrain2(_NaturalStrainBase):
    """Second-order natural-strain EOS, with implied ``K0_prime = 2``."""

    def _coefficients(self) -> tuple[float, float]:
        return 0.0, 0.0

    def __init__(self, V0: float, K0: float) -> None:
        super().__init__(V0, K0)
        self._native = _rust.RtEos.natural_strain2(self.V0, self.K0)


class NaturalStrain3(_NaturalStrainBase):
    """Third-order natural-strain EOS.

    Reference: Poirier, J.-P. & Tarantola, A. (1998), Physics of the Earth
    and Planetary Interiors 109, 1-8, doi:10.1016/S0031-9201(98)00112-5.
    """

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        super().__init__(V0, K0)
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self._native = _rust.RtEos.natural_strain3(self.V0, self.K0, self.K0_prime)

    def _coefficients(self) -> tuple[float, float]:
        return 1.5 * (self.K0_prime - 2.0), 0.0


class NaturalStrain4(NaturalStrain3):
    """Fourth-order natural-strain EOS including ``K0_double_prime``."""

    def __init__(
        self,
        V0: float,
        K0: float,
        K0_prime: float,
        K0_double_prime: float,
    ) -> None:
        super().__init__(V0, K0, K0_prime)
        self.K0_double_prime = validate_finite_scalar(
            K0_double_prime, "K0_double_prime"
        )
        self._native = _rust.RtEos.natural_strain4(
            self.V0, self.K0, self.K0_prime, self.K0_double_prime
        )

    def _coefficients(self) -> tuple[float, float]:
        difference = self.K0_prime - 2.0
        a = 1.5 * difference
        b = 1.5 * (self.K0 * self.K0_double_prime + 1.0 + difference + difference**2)
        return a, b
