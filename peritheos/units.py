"""Explicit unit conversions used by Peritheos EOS workflows.

All helpers accept either scalars or NumPy arrays and preserve that distinction
in their return value.  EOS calculations themselves continue to use the
documented GPa, K, and volume conventions; this module is the single public
location for converting inputs into those conventions.
"""

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

_PRESSURE_TO_PA = {
    "pa": 1.0,
    "mpa": 1.0e6,
    "gpa": 1.0e9,
    "bar": 1.0e5,
    "kbar": 1.0e8,
    "atm": 101325.0,
    "torr": 133.322,
    "psi": 6894.76,
}

_ANGSTROM_CUBED_TO_M3 = 1.0e-30
_AVOGADRO = 6.02214076e23


def _convert(
    value: Numeric, source: str, target: str, factors: dict[str, float]
) -> Numeric:
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


def convert_pressure(value: Numeric, from_unit: str, to_unit: str) -> Numeric:
    """Convert pressure among Pa, MPa, GPa, bar, kbar, atm, torr, and psi."""
    if _normalize(from_unit) not in _PRESSURE_TO_PA:
        raise ValueError(f"Unsupported pressure unit: {from_unit}")
    if _normalize(to_unit) not in _PRESSURE_TO_PA:
        raise ValueError(f"Unsupported pressure unit: {to_unit}")
    return _convert(value, from_unit, to_unit, _PRESSURE_TO_PA)


def convert_temperature(value: Numeric, from_unit: str, to_unit: str) -> Numeric:
    """Convert temperature among kelvin, degrees Celsius, and Fahrenheit."""
    source = from_unit.lower().replace("°", "")
    target = to_unit.lower().replace("°", "")
    values = np.asarray(value, dtype=float)
    if not np.all(np.isfinite(values)):
        raise ValueError("Values must be finite")

    if source == "k":
        kelvin = values
    elif source == "c":
        kelvin = values + 273.15
    elif source == "f":
        kelvin = (values - 32.0) * 5.0 / 9.0 + 273.15
    else:
        raise ValueError(f"Unsupported temperature unit: {from_unit}")

    if target == "k":
        result = kelvin
    elif target == "c":
        result = kelvin - 273.15
    elif target == "f":
        result = (kelvin - 273.15) * 9.0 / 5.0 + 32.0
    else:
        raise ValueError(f"Unsupported temperature unit: {to_unit}")
    if result.ndim == 0:
        return float(result)
    return result


def cell_volume_to_molar_volume(
    cell_volume: Numeric,
    formula_units_per_cell: Numeric,
    *,
    to_unit: str = "J/bar/mol",
) -> Numeric:
    """Convert conventional-cell volume in angstrom^3 to formula molar volume.

    ``formula_units_per_cell`` is the crystallographic number of formula units
    in the conventional cell (often written ``Z``).  The result is per mole of
    formula units, not per mole of atoms.
    """
    cell_values = np.asarray(cell_volume, dtype=float)
    formula_units = np.asarray(formula_units_per_cell, dtype=float)
    try:
        cell_values, formula_units = np.broadcast_arrays(cell_values, formula_units)
    except ValueError as error:
        raise ValueError(
            "cell_volume and formula_units_per_cell must have broadcast-compatible shapes"
        ) from error
    if not np.all(np.isfinite(cell_values)) or np.any(cell_values <= 0.0):
        raise ValueError("Cell volume must be finite and greater than zero")
    if not np.all(np.isfinite(formula_units)) or np.any(formula_units <= 0.0):
        raise ValueError("Formula units per cell must be finite and greater than zero")
    molar_volume_m3 = cell_values * _ANGSTROM_CUBED_TO_M3 * _AVOGADRO / formula_units
    return convert_molar_volume(molar_volume_m3, "m3/mol", to_unit)


def molar_volume_to_cell_volume(
    molar_volume: Numeric,
    formula_units_per_cell: Numeric,
    *,
    from_unit: str = "J/bar/mol",
) -> Numeric:
    """Convert formula molar volume to conventional-cell volume in angstrom^3."""
    molar_volume_m3 = np.asarray(
        convert_molar_volume(molar_volume, from_unit, "m3/mol"), dtype=float
    )
    formula_units = np.asarray(formula_units_per_cell, dtype=float)
    try:
        molar_volume_m3, formula_units = np.broadcast_arrays(
            molar_volume_m3, formula_units
        )
    except ValueError as error:
        raise ValueError(
            "molar_volume and formula_units_per_cell must have broadcast-compatible shapes"
        ) from error
    if np.any(molar_volume_m3 <= 0.0):
        raise ValueError("Molar volume must be greater than zero")
    if not np.all(np.isfinite(formula_units)) or np.any(formula_units <= 0.0):
        raise ValueError("Formula units per cell must be finite and greater than zero")
    result = molar_volume_m3 * formula_units / (_ANGSTROM_CUBED_TO_M3 * _AVOGADRO)
    if result.ndim == 0:
        return float(result)
    return result


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
