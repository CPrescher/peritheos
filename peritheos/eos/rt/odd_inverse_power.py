"""Odd inverse-volume series used by Anderson and Swenson (1975)."""

import numpy as np

from peritheos.eos import (
    EosBase,
    NumericType,
    validate_finite_scalar,
    validate_positive_scalar,
    validate_volume,
)
from peritheos.errors import EosNumericalError


class OddInversePower(EosBase):
    r"""Evaluate ``P = sum(Cn * (V0/V)**n)`` for n=3,5,7,9.

    ``V0`` is a volume normalization, not necessarily the exact zero-pressure
    root. Coefficients have pressure units. For the source's molar coefficients
    An, set Cn=An/V0_molar**n (and convert kbar to GPa). This preserves the
    published polynomial, including rounding; no baseline is subtracted.
    Only the source's finite compression interval is physically warranted.
    Python uses NumPy; the Rust core implements the same family independently.
    """

    def __init__(self, V0: float, C3: float, C5: float, C7: float, C9: float):
        self.V0 = validate_positive_scalar(V0, "V0")
        self.C3 = validate_finite_scalar(C3, "C3")
        self.C5 = validate_finite_scalar(C5, "C5")
        self.C7 = validate_finite_scalar(C7, "C7")
        self.C9 = validate_finite_scalar(C9, "C9")

    def _evaluate(self, V: NumericType, derivative: bool) -> NumericType:
        x = self.V0 / np.asarray(validate_volume(V), dtype=float)
        with np.errstate(over="ignore", invalid="ignore"):
            result = sum(
                c * (n if derivative else 1) * x**n
                for n, c in ((3, self.C3), (5, self.C5), (7, self.C7), (9, self.C9))
            )
        if not np.all(np.isfinite(result)):
            raise EosNumericalError("Odd inverse-volume series overflow")
        return float(result) if result.ndim == 0 else result

    def pressure(self, V: NumericType) -> NumericType:
        return self._evaluate(V, False)

    def bulk_modulus(self, V: NumericType) -> NumericType:
        return self._evaluate(V, True)
