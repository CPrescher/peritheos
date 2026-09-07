"""Reproduce the Matsui et al. (2000) MgO Table 3 EOS checks."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

V0_CELL_A3 = 11.2382 * 4.0e24 / 6.02214076e23
K0_GPA = 160.5
K0_PRIME = 4.10
DATASET = (
    Path(__file__).parents[1] / "peritheos/data/datasets/mgo-matsui-2000-table3-pvt.csv"
)


def bm3_pressure(volume_ratio: np.ndarray) -> np.ndarray:
    """Evaluate the standard third-order Birch-Murnaghan isotherm."""
    x = volume_ratio ** (-1.0 / 3.0)
    return 1.5 * K0_GPA * (x**7 - x**5) * (1.0 + 0.75 * (K0_PRIME - 4.0) * (x**2 - 1.0))


def main() -> None:
    with DATASET.open(encoding="utf-8", newline="") as stream:
        rows = [row for row in csv.DictReader(stream) if row["row_kind"] == "md"]

    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    temperature = np.array([float(row["temperature_k"]) for row in rows])
    volume_ratio = np.array([float(row["volume_ratio"]) for row in rows])
    cold_pressure = bm3_pressure(volume_ratio)

    eta = 1.0 - volume_ratio
    delta_temperature = temperature - 300.0
    design = np.column_stack(
        [
            np.ones_like(eta),
            eta,
            delta_temperature,
            0.5 * eta**2,
            0.5 * delta_temperature**2,
            0.5 * eta * delta_temperature,
        ]
    )
    coefficients, *_ = np.linalg.lstsq(design, pressure - cold_pressure, rcond=None)
    residuals = cold_pressure + design @ coefficients - pressure

    source_selection = (temperature == 300.0) & (pressure <= 20.0)
    source_residuals = cold_pressure[source_selection] - pressure[source_selection]

    print(f"V0 conventional B1 cell (A^3): {V0_CELL_A3:.14g}")
    print("Taylor coefficients:", " ".join(f"{value:.16g}" for value in coefficients))
    print(
        "Source BM3, 300 K and 0-20 GPa: "
        f"RMSE={np.sqrt(np.mean(source_residuals**2)):.12g} GPa, "
        f"max_abs={np.max(np.abs(source_residuals)):.12g} GPa"
    )
    print(
        "Thermal refit, all 72 MD states: "
        f"RMSE={np.sqrt(np.mean(residuals**2)):.12g} GPa, "
        f"max_abs={np.max(np.abs(residuals)):.12g} GPa"
    )


if __name__ == "__main__":
    main()
