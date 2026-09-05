#!/usr/bin/env python3
"""Reproduce the defensible volumetric EOS from Criniti et al. (2021)."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "bridgmanite-criniti-2021-table1-density.csv"
)

MOLAR_MASS_MGSIO3_G_MOL = 100.3875
AVOGADRO_MOL_MINUS_1 = 6.02214076e23
FORMULA_UNITS_PER_CELL = 4
PUBLISHED = {"V0": 162.453274, "K0": 254.5, "K0_prime": 3.73}


def cell_volume_from_density(density_g_cm3: np.ndarray) -> np.ndarray:
    """Convert MgSiO3 density to the conventional Pbnm Z=4 cell volume."""
    return (
        FORMULA_UNITS_PER_CELL
        * MOLAR_MASS_MGSIO3_G_MOL
        / density_g_cm3
        / AVOGADRO_MOL_MINUS_1
        * 1.0e24
    )


def bm3_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    """Evaluate the standard Eulerian third-order Birch-Murnaghan equation."""
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def reproduce() -> dict[str, object]:
    """Return deterministic conversion and published-curve diagnostics."""
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    density = np.array([float(row["density_g_cm3"]) for row in rows])
    stored_volume = np.array(
        [float(row["volume_a3_conventional_cell"]) for row in rows]
    )
    converted_volume = cell_volume_from_density(density)
    calculated_pressure = bm3_pressure(
        stored_volume,
        PUBLISHED["V0"],
        PUBLISHED["K0"],
        PUBLISHED["K0_prime"],
    )
    residual = calculated_pressure - pressure
    return {
        "observations": len(rows),
        "published_parameters": PUBLISHED,
        "maximum_volume_conversion_difference_a3": float(
            np.max(np.abs(converted_volume - stored_volume))
        ),
        "published_curve_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_curve_maximum_absolute_pressure_residual_gpa": float(
            np.max(np.abs(residual))
        ),
    }


def main() -> None:
    """Print stable JSON diagnostics."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
