"""Morse and Sun Jiu-Xun--Morse isothermal equations of state."""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)


class _ThreeParameterMorseBase(EosBase):
    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        self.V0 = validate_positive_scalar(V0, "V0")
        self.K0 = validate_positive_scalar(K0, "K0")
        self.K0_prime = validate_finite_scalar(K0_prime, "K0_prime")

    def pressure(self, V: NumericType) -> NumericType:
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "pressure", V)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        V = validate_volume(V)
        return _native_rt_evaluate(self._native, "bulk_modulus", V)


class Morse3(_ThreeParameterMorseBase):
    r"""Three-dimensional Morse-potential EOS (Sun et al. 2010, equation 6)."""

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        super().__init__(V0, K0, K0_prime)
        self._native = _rust.RtEos.morse3(self.V0, self.K0, self.K0_prime)


class SunMorse3(_ThreeParameterMorseBase):
    r"""Sun Jiu-Xun--Morse EOS with ``n=3`` (Sun et al. 2010, equation 8)."""

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        super().__init__(V0, K0, K0_prime)
        self._native = _rust.RtEos.sun_morse3(self.V0, self.K0, self.K0_prime)


class SunMorse4(_ThreeParameterMorseBase):
    r"""Sun Jiu-Xun--Morse EOS with ``n=4`` (Sun et al. 2010, equation 8)."""

    def __init__(self, V0: float, K0: float, K0_prime: float) -> None:
        super().__init__(V0, K0, K0_prime)
        self._native = _rust.RtEos.sun_morse4(self.V0, self.K0, self.K0_prime)
