"""Display-only rounding for scientific reports; never change stored numbers."""

from __future__ import annotations

import math
from decimal import ROUND_HALF_UP, Decimal, localcontext
from typing import Any


def number(value: Any, significant: int = 6) -> str:
    if value is None:
        return "Not reported"
    if isinstance(value, bool) or not isinstance(value, (int, float)):
        return str(value)
    if not math.isfinite(value):
        return "Not available"
    return format(value, f".{significant}g")


def value_error(value: Any, error: Any) -> tuple[str, str]:
    """Round error to two significant digits for leading 1/2, otherwise one.

    Match the value's last decimal place to that of the error. Decimal avoids
    binary tie artefacts; retain the selected place when rounding carries.
    Missing/zero errors do not license extra precision or an exact-value claim.
    """
    if (
        isinstance(error, bool)
        or not isinstance(error, (int, float))
        or not math.isfinite(error)
        or error <= 0
    ):
        label = "0 (reported)" if error == 0 else "Not reported"
        return number(value), label
    with localcontext() as context:
        context.prec = 50
        uncertainty = Decimal(str(error))
        exponent = uncertainty.adjusted()
        leading = int(uncertainty.scaleb(-exponent))
        place = exponent - (1 if leading in (1, 2) else 0)
        quantum = Decimal(1).scaleb(place)
        rounded_error = uncertainty.quantize(quantum, rounding=ROUND_HALF_UP)
        if isinstance(value, (int, float)) and math.isfinite(value):
            rounded_value = Decimal(str(value)).quantize(
                quantum, rounding=ROUND_HALF_UP
            )
            if rounded_value == 0:
                rounded_value = abs(rounded_value)
            value_text = format(rounded_value, "f")
        else:
            value_text = number(value)
        # Use a common explicit exponent for very small/large quantities.
        if exponent < -4 or exponent >= 6:
            digits = max(0, exponent - place)
            if isinstance(value, (int, float)) and math.isfinite(value):
                value_text = (
                    f"{rounded_value.scaleb(-exponent):.{digits}f}e{exponent:+d}"
                )
            error_text = f"{rounded_error.scaleb(-exponent):.{digits}f}e{exponent:+d}"
        else:
            error_text = format(rounded_error, "f")
    return value_text, error_text


def rounded_metadata(value: Any) -> Any:
    """Keep protocol text intact, but avoid dumping float tails in diagnostics."""
    if isinstance(value, dict):
        result = {key: rounded_metadata(item) for key, item in value.items()}
        # Auxiliary parameter comparisons (for example acoustic coefficients)
        # live inside protocol metadata rather than the main coefficient table.
        for side in ("published", "refit"):
            if side in value:
                for error_key in (side + "_error", side + "_sigma"):
                    if error_key in value:
                        result[side], result[error_key] = value_error(
                            value[side], value[error_key]
                        )
                        break
        return result
    if isinstance(value, (list, tuple)):
        return [rounded_metadata(item) for item in value]
    if isinstance(value, float):
        return number(value)
    return value
