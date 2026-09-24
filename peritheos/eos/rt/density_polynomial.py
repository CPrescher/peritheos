"""Published cubic pressure-density polynomials without an ambient reference."""

from peritheos.eos import (
    EosBase,
    NumericType,
    _native_rt_evaluate,
    _rust,
    validate_positive_scalar,
    validate_volume,
)


class DensityPolynomial3(EosBase):
    """P = c0 + c1*rho + c2*rho**2 + c3*rho**3, rho=rho0*V0/V.

    Pressure is in GPa and density in g/cm^3. V0 is a normalization anchor,
    not a zero-pressure volume. Only the stable positive-compressibility branch
    is supported. c2 must be nonnegative and c3 positive. No thermal reference
    shift or implied adiabatic/isothermal conversion is applied.
    """

    def __init__(
        self, V0: float, rho0: float, c0: float, c1: float, c2: float, c3: float
    ) -> None:
        self.V0 = validate_positive_scalar(V0, "V0")
        self.rho0 = validate_positive_scalar(rho0, "rho0")
        self.c0, self.c1, self.c2, self.c3 = map(float, (c0, c1, c2, c3))
        self._native = _rust.RtEos.density_polynomial_3(
            self.V0, self.rho0, self.c0, self.c1, self.c2, self.c3
        )

    def pressure(self, V: NumericType) -> NumericType:
        return _native_rt_evaluate(self._native, "pressure", validate_volume(V))

    def bulk_modulus(self, V: NumericType) -> NumericType:
        return _native_rt_evaluate(self._native, "bulk_modulus", validate_volume(V))
