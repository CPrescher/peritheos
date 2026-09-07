"""Baonza pseudospinodal equation of state."""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_positive_scalar,
    validate_volume,
)


class Baonza(EosBase):
    r"""Three-parameter Baonza EOS with the published fixed ``beta=0.85``.

    The model is defined through ``B = (P-P_sp)^beta / kappa_star`` and its
    integrated volume-pressure relation.  ``P_sp`` and ``kappa_star`` are
    determined by ``V0``, ``K0``, and ``K0_prime``.
    """

    beta = 0.85

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_positive_scalar(K0_prime, "K0_prime")
        self._native = _rust.RtEos.baonza(self.V0, self.K0, self.K0_prime)

    def pressure(self, V: NumericType) -> NumericType:
        """Return pressure at volume ``V``."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "pressure", V)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """Return the isothermal bulk modulus at volume ``V``."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "bulk_modulus", V)
