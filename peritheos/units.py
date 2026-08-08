"""Explicit unit conversions used by Peritheos EOS workflows."""

from typing import Union

import numpy as np
from numpy.typing import NDArray

Numeric = Union[float, NDArray[np.float64]]


def _normalize(unit: str) -> str:
    normalized = (
        unit.lower()
        .replace(" ", "")
        .replace("³", "3")
        .replace("^", "")
        .replace("angstrom", "a")
        .replace("å", "a")
    )
    return normalized.replace("bar-1", "/bar").replace("mol-1", "/mol")


_MOLAR_VOLUME_TO_M3 = {
    "m3/mol": 1.0,
    "cm3/mol": 1.0e-6,
    "l/mol": 1.0e-3,
    "j/bar/mol": 1.0e-5,
}

_DENSITY_TO_KG_M3 = {
    "kg/m3": 1.0,
    "g/cm3": 1000.0,
}

_MOLAR_MASS_TO_KG_MOL = {
    "kg/mol": 1.0,
    "g/mol": 1.0e-3,
}


def _convert(value: Numeric, source: str, target: str, factors: dict[str, float]):
    source_key = _normalize(source)
    target_key = _normalize(target)
    if source_key not in factors:
        raise ValueError(f"Unsupported unit: {source}")
    if target_key not in factors:
        raise ValueError(f"Unsupported unit: {target}")
    result = np.asarray(value, dtype=float) * factors[source_key] / factors[target_key]
    if not np.all(np.isfinite(result)):
        raise ValueError("Values must be finite")
    if result.ndim == 0:
        return float(result)
    return result


def convert_molar_volume(value: Numeric, from_unit: str, to_unit: str) -> Numeric:
    """Convert molar volume among m3/mol, cm3/mol, L/mol, and J/bar/mol."""
    return _convert(value, from_unit, to_unit, _MOLAR_VOLUME_TO_M3)


def convert_density(value: Numeric, from_unit: str, to_unit: str) -> Numeric:
    """Convert density between kg/m3 and g/cm3."""
    return _convert(value, from_unit, to_unit, _DENSITY_TO_KG_M3)


def density_from_molar_volume(
    molar_mass: Numeric,
    molar_volume: Numeric,
    *,
    molar_mass_unit: str = "g/mol",
    volume_unit: str = "cm3/mol",
    density_unit: str = "g/cm3",
) -> Numeric:
    """Calculate density from molar mass and molar volume."""
    mass = _convert(molar_mass, molar_mass_unit, "kg/mol", _MOLAR_MASS_TO_KG_MOL)
    volume = convert_molar_volume(molar_volume, volume_unit, "m3/mol")
    if np.any(np.asarray(mass) <= 0.0):
        raise ValueError("Molar mass must be greater than zero")
    if np.any(np.asarray(volume) <= 0.0):
        raise ValueError("Molar volume must be greater than zero")
    density = np.asarray(mass, dtype=float) / np.asarray(volume, dtype=float)
    return convert_density(density, "kg/m3", density_unit)


def molar_volume_from_density(
    molar_mass: Numeric,
    density: Numeric,
    *,
    molar_mass_unit: str = "g/mol",
    density_unit: str = "g/cm3",
    volume_unit: str = "cm3/mol",
) -> Numeric:
    """Calculate molar volume from molar mass and density."""
    mass = _convert(molar_mass, molar_mass_unit, "kg/mol", _MOLAR_MASS_TO_KG_MOL)
    density_si = convert_density(density, density_unit, "kg/m3")
    if np.any(np.asarray(mass) <= 0.0):
        raise ValueError("Molar mass must be greater than zero")
    if np.any(np.asarray(density_si) <= 0.0):
        raise ValueError("Density must be greater than zero")
    volume = np.asarray(mass, dtype=float) / np.asarray(density_si, dtype=float)
    return convert_molar_volume(volume, "m3/mol", volume_unit)
