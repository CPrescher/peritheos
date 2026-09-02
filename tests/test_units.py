"""Tests for explicit molar-volume and density conversions."""

import numpy as np
import pytest

from peritheos.units import (
    cell_volume_to_molar_volume,
    convert_density,
    convert_molar_volume,
    convert_pressure,
    convert_temperature,
    density_from_molar_volume,
    molar_volume_from_density,
    molar_volume_to_cell_volume,
)


def test_thermal_volume_convention_conversion():
    assert np.isclose(convert_molar_volume(10.0, "cm^3/mol", "J/bar/mol"), 1.0)
    assert np.isclose(convert_molar_volume(1.0, "J bar^-1 mol^-1", "cm3/mol"), 10.0)


def test_array_conversion():
    values = np.array([1.0, 2.0])
    assert np.allclose(convert_density(values, "g/cm3", "kg/m3"), values * 1000.0)


def test_pressure_and_temperature_conversions_are_array_aware():
    assert np.allclose(
        convert_pressure(np.array([1.0, 2.0]), "GPa", "kbar"),
        np.array([10.0, 20.0]),
    )
    assert np.allclose(
        convert_temperature(np.array([273.15, 373.15]), "K", "°C"),
        np.array([0.0, 100.0]),
    )


def test_cell_and_formula_molar_volume_round_trip():
    # MgO has four formula units in the conventional B1 cell.
    cell_volume = 74.698
    molar_volume = cell_volume_to_molar_volume(cell_volume, 4)

    assert molar_volume == pytest.approx(1.1246046762)
    assert molar_volume_to_cell_volume(molar_volume, 4) == pytest.approx(cell_volume)


def test_cell_volume_conversion_broadcasts_formula_unit_counts():
    cell_volumes = np.array([40.0, 80.0])
    formula_units = np.array([2.0, 4.0])

    result = cell_volume_to_molar_volume(cell_volumes, formula_units, to_unit="cm3/mol")

    assert np.allclose(result, result[0])


def test_density_volume_round_trip():
    density = density_from_molar_volume(40.304, 11.25)
    volume = molar_volume_from_density(40.304, density)

    assert np.isclose(volume, 11.25)


@pytest.mark.parametrize("unit", ["cubits/mol", "", "cm/mol"])
def test_unsupported_volume_unit(unit):
    with pytest.raises(ValueError, match="Unsupported unit"):
        convert_molar_volume(1.0, unit, "cm3/mol")


def test_nonpositive_density_and_volume_are_rejected():
    with pytest.raises(ValueError, match="Molar volume"):
        density_from_molar_volume(40.0, 0.0)
    with pytest.raises(ValueError, match="Density"):
        molar_volume_from_density(40.0, 0.0)
    with pytest.raises(ValueError, match="Molar mass"):
        density_from_molar_volume(0.0, 10.0)


@pytest.mark.parametrize("formula_units", [0.0, -1.0, np.nan, np.inf])
def test_cell_volume_conversion_rejects_invalid_formula_unit_count(formula_units):
    with pytest.raises(ValueError, match="Formula units per cell"):
        cell_volume_to_molar_volume(10.0, formula_units)


def test_cell_volume_conversion_rejects_nonpositive_volume():
    with pytest.raises(ValueError, match="Cell volume"):
        cell_volume_to_molar_volume(0.0, 1.0)


def test_cell_volume_conversion_reports_shape_mismatch():
    with pytest.raises(ValueError, match="broadcast-compatible"):
        cell_volume_to_molar_volume(np.ones(2), np.ones(3))
    with pytest.raises(ValueError, match="broadcast-compatible"):
        molar_volume_to_cell_volume(np.ones(2), np.ones(3))
