"""Public exception hierarchy shared by Python and the Rust extension.

Exception classes and their ``code`` values are the machine-readable error
contract.  Exception messages are intended for people and may become more
specific between releases.
"""

from __future__ import annotations

from collections.abc import Mapping
from types import MappingProxyType
from typing import Any


class PeritheosError(Exception):
    """Base class for errors raised deliberately by Peritheos.

    Parameters beyond ``message`` are keyword-only so adding context does not
    change the conventional single-string ``Exception.args`` value.
    """

    default_code = "peritheos.error"

    def __init__(
        self,
        message: str,
        *,
        code: str | None = None,
        operation: str | None = None,
        field: str | None = None,
        context: Mapping[str, Any] | None = None,
    ) -> None:
        super().__init__(message)
        self.code = code or self.default_code
        self.operation = operation
        self.field = field
        self.context = MappingProxyType(dict(context or {}))


class ValidationError(PeritheosError, ValueError):
    """Invalid user input, parameters, shapes, units, or state."""

    default_code = "validation.invalid_input"


class NumericalError(PeritheosError, ArithmeticError):
    """A valid request failed numerically or produced a non-finite result."""

    default_code = "numerical.failure"


class ConfigurationError(PeritheosError, TypeError):
    """Objects or options cannot be combined in the requested configuration."""

    default_code = "configuration.invalid"


class UnsupportedOperationError(PeritheosError, NotImplementedError):
    """The selected model does not implement the requested operation."""

    default_code = "operation.unsupported"


class EosError(PeritheosError):
    """Base class for equation-of-state construction and evaluation errors."""

    default_code = "eos.error"


class EosValidationError(ValidationError, EosError):
    """An EOS parameter or requested thermodynamic state is invalid."""

    default_code = "eos.invalid_input"


class EosNumericalError(NumericalError, EosError):
    """An EOS evaluation or inversion failed numerically."""

    default_code = "eos.numerical_failure"


class FitError(PeritheosError, RuntimeError):
    """Base class for fitting failures."""

    default_code = "fit.error"


class FitValidationError(ValidationError, FitError):
    """Fit observations, bounds, options, or dimensions are invalid."""

    default_code = "fit.invalid_input"


class FitNumericalError(NumericalError, FitError):
    """A fit failed because a numerical system was singular or unstable."""

    default_code = "fit.numerical_failure"


class FitEosValidationError(FitValidationError, EosValidationError):
    """An EOS used by a fit rejected a parameter or state."""

    default_code = "eos.invalid_input"


class FitEosNumericalError(FitNumericalError, EosNumericalError):
    """An EOS used by a fit failed numerically."""

    default_code = "eos.numerical_failure"


class EosmatError(ValidationError):
    """An ``.eosmat`` document is malformed, invalid, or unsupported."""

    default_code = "eosmat.invalid_document"


class MaterialError(ValidationError):
    """A material or EOS record violates the library data contract."""

    default_code = "material.invalid"


class MaterialLookupError(PeritheosError, KeyError):
    """A requested bundled material or EOS record does not exist."""

    default_code = "material.not_found"


__all__ = [
    "ConfigurationError",
    "EosError",
    "EosNumericalError",
    "EosValidationError",
    "EosmatError",
    "FitError",
    "FitEosNumericalError",
    "FitEosValidationError",
    "FitNumericalError",
    "FitValidationError",
    "MaterialError",
    "MaterialLookupError",
    "NumericalError",
    "PeritheosError",
    "UnsupportedOperationError",
    "ValidationError",
]
