"""Murnaghan equation of state."""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class Murnaghan(EosBase):
    """Murnaghan isothermal equation of state.

    The model assumes that the bulk modulus varies linearly with pressure,
    ``K(P) = K0 + K0_prime * P``. Its pressure-volume relation is

    ``P(V) = K0 / K0_prime * ((V0 / V)**K0_prime - 1)``.

    The continuous ``K0_prime = 0`` limit is supported and gives
    ``P(V) = K0 * log(V0 / V)``.

    Reference
    ---------
    Murnaghan, F. D. (1944). The compressibility of media under extreme
    pressures. Proceedings of the National Academy of Sciences, 30, 244-247.
    """

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        """Initialize the Murnaghan EOS.

        Parameters
        ----------
        V0 : float
            Reference volume in any consistent volume unit.
        K0 : float
            Bulk modulus at the reference volume.
        K0_prime : float
            Pressure derivative of the reference bulk modulus.
        """
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")
        self._native = _rust.RtEos.murnaghan(self.V0, self.K0, self.K0_prime)

    def pressure(self, V: NumericType) -> NumericType:
        """Return pressure at volume *V* in the same units as ``K0``."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "pressure", V)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        """Return the isothermal bulk modulus at volume *V*."""
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "bulk_modulus", V)
