#!/usr/bin/env python3
"""Reproduce the Fiquet et al. (2000) 298 K bridgmanite BM3."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "bridgmanite-fiquet-2000-table1-298k-pv.csv"
)
PUBLISHED = {"V0": 162.27, "K0": 253.0, "K0_prime": 3.9}


def birch_murnaghan_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def load_data() -> tuple[np.ndarray, np.ndarray]:
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return (
        np.array([float(row["pressure_gpa"]) for row in rows]),
        np.array([float(row["cell_volume_a3"]) for row in rows]),
    )


def reproduce() -> dict[str, object]:
    pressure, volume = load_data()
    published_pressure = birch_murnaghan_pressure(
        volume, PUBLISHED["V0"], PUBLISHED["K0"], PUBLISHED["K0_prime"]
    )
    published_residual = published_pressure - pressure

    fixed_v0_fit = least_squares(
        lambda parameters: (
            birch_murnaghan_pressure(
                volume, PUBLISHED["V0"], parameters[0], parameters[1]
            )
            - pressure
        ),
        np.array([PUBLISHED["K0"], PUBLISHED["K0_prime"]]),
    )
    free_fit = least_squares(
        lambda parameters: (
            birch_murnaghan_pressure(
                volume, parameters[0], parameters[1], parameters[2]
            )
            - pressure
        ),
        np.array([PUBLISHED["V0"], PUBLISHED["K0"], PUBLISHED["K0_prime"]]),
    )
    return {
        "observations": int(pressure.size),
        "pressure_range_gpa": [float(pressure.min()), float(pressure.max())],
        "published_parameters": PUBLISHED,
        "published_curve_pressure_rmse_gpa": float(
            np.sqrt(np.mean(published_residual**2))
        ),
        "published_curve_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(published_residual))
        ),
        "fixed_v0_refit": {
            "V0": PUBLISHED["V0"],
            "K0": float(fixed_v0_fit.x[0]),
            "K0_prime": float(fixed_v0_fit.x[1]),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(fixed_v0_fit.fun**2))),
        },
        "unweighted_all_free_diagnostic": {
            "V0": float(free_fit.x[0]),
            "K0": float(free_fit.x[1]),
            "K0_prime": float(free_fit.x[2]),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(free_fit.fun**2))),
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
