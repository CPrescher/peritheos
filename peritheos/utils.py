"""General numerical helpers and compatibility unit-conversion imports."""

from __future__ import annotations

import warnings
from collections.abc import Callable

import peritheos.units as _units
from peritheos.units import Numeric

from .constants import R


def _warn_unit_move(name: str) -> None:
    warnings.warn(
        f"peritheos.utils.{name} is deprecated; import it from peritheos.units",
        DeprecationWarning,
        stacklevel=2,
    )


def convert_pressure(value: Numeric, from_unit: str, to_unit: str) -> Numeric:
    """Deprecated compatibility wrapper for :func:`peritheos.units.convert_pressure`."""
    _warn_unit_move("convert_pressure")
    return _units.convert_pressure(value, from_unit, to_unit)


def convert_temperature(value: Numeric, from_unit: str, to_unit: str) -> Numeric:
    """Deprecated compatibility wrapper for :func:`peritheos.units.convert_temperature`."""
    _warn_unit_move("convert_temperature")
    return _units.convert_temperature(value, from_unit, to_unit)


def compressibility_factor(
    pressure: Numeric,
    volume: Numeric,
    temperature: Numeric,
    moles: Numeric,
) -> Numeric:
    """
    Calculate the compressibility factor Z = PV/nRT

    Parameters
    ----------
    pressure : float
        Pressure in Pascal
    volume : float
        Volume in cubic meters
    temperature : float
        Temperature in Kelvin
    moles : float
        Number of moles

    Returns
    -------
    float
        Compressibility factor (dimensionless)
    """
    return pressure * volume / (moles * R * temperature)


def derivative(
    f: Callable[[Numeric], Numeric],
    x: Numeric,
    dx: float = 1e-6,
) -> Numeric:
    """Compute the derivative of f at x using finite differences"""
    return (f(x + dx) - f(x - dx)) / (2 * dx)
