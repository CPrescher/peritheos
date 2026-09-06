"""Generalized Rydberg--Stacey equation of state."""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class RydbergStacey(EosBase):
    """Generalized Rydberg--Stacey isothermal equation of state.

    With ``x = (V / V0)**(1/3)`` and
    ``a = 3*K0_prime/2 - 3*K_infinity_prime + 1/2``, the pressure is

    ``P = 3*K0*x**(-3*K_infinity_prime)*(1-x)*exp(a*(1-x))``.

    ``K_infinity_prime`` is the limiting pressure derivative of the bulk
    modulus at infinite pressure. The conventional Rydberg--Vinet equation is
    recovered when ``K_infinity_prime = 2/3``.
    """

    def __init__(
        self,
        V0: float,
        K0: float,
        K0_prime: float,
        K_infinity_prime: float,
    ) -> None:
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self.K_infinity_prime = validate_finite_scalar(
            K_infinity_prime, "K_infinity_prime"
        )
        self._native = _rust.RtEos.rydberg_stacey(
            self.V0,
            self.K0,
            self.K0_prime,
            self.K_infinity_prime,
        )

    def pressure(self, V: NumericType) -> NumericType:
        """Return pressure at volume *V* in the same units as ``K0``."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "pressure", V)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """Return the isothermal bulk modulus at volume *V*."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "bulk_modulus", V)
