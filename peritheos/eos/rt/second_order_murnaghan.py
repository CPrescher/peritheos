"""Second-order Murnaghan with a static pressure offset (Finger 1981 Eq. 6)."""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class SecondOrderMurnaghan(EosBase):
    """Quadratic K(P-P0) equation; distinct from second-order Birch-Murnaghan.

    Supports the real positive-q branch, q^2 = K0_prime^2 - 2*K0*K0_double_prime.
    P0 is the pressure at V0. Volumes use any consistent unit; K0 and P0 use
    GPa, and K0_double_prime uses GPa^-1.
    """

    def __init__(self, V0, K0, K0_prime, K0_double_prime, P0=0.0):
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self.K0_double_prime = validate_finite_scalar(
            K0_double_prime, "K0_double_prime"
        )
        self.P0 = validate_finite_scalar(P0, "P0")
        self._native = _rust.RtEos.second_order_murnaghan(
            self.V0, self.K0, self.K0_prime, self.K0_double_prime, self.P0
        )

    def pressure(self, V: NumericType) -> NumericType:
        return _native_rt_evaluate(self._native, "pressure", validate_volume(V))

    def bulk_modulus(self, V: NumericType) -> NumericType:
        return _native_rt_evaluate(self._native, "bulk_modulus", validate_volume(V))
