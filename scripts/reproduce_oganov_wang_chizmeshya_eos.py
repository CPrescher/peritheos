#!/usr/bin/env python3
"""Numerically audit the Oganov (2003), Wang (1996), and Chizmeshya (1996) BM3 records."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OGANOV_DATA = ROOT / "peritheos/data/datasets/mgo-oganov-2003-table3-pv.csv"
WANG_DATA = (
    ROOT / "peritheos/data/datasets/ca-perovskite-wang-1996-table1-room-temperature.csv"
)
MAO_DATA = (
    ROOT / "peritheos/data/datasets/ca-perovskite-mao-1989-table1-compression.csv"
)

OGANOV = {
    "ecp_large_core_static": (77.629, 151.707, 4.212),
    "ecp_small_core_static": (76.595, 150.839, 4.052),
    "paw_large_core_static": (76.049, 154.183, 4.141),
    "paw_small_core_static": (76.947, 150.597, 4.103),
    "pressure_corrected_static": (73.425, 181.240, 3.997),
    "pressure_corrected_0k": (74.439, 173.480, 4.014),
    "pressure_corrected_298k": (74.670, 170.530, 4.036),
    "pressure_corrected_1000k": (76.549, 152.595, 4.130),
    "pressure_corrected_2000k": (79.915, 127.719, 4.244),
    "pressure_corrected_3000k": (83.772, 106.110, 4.331),
    "pressure_corrected_4000k": (88.006, 88.473, 4.385),
}

CHIZMESHYA = {
    "lapw7_static": (45.04, 241.8, 4.15),
    "lapw8_static": (45.02, 241.0, 4.16),
    "lapw9_static": (45.06, 238.2, 4.18),
    "lapw9_300k_kp4": (45.55, 237.9, 4.0),
    "lapw9_300k": (45.62, 227.0, 4.29),
}


def bm3_pressure(volume: np.ndarray | float, v0: float, k0: float, kp: float):
    eta = (v0 / np.asarray(volume)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def _read(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _fit_wang_room_temperature() -> dict[str, float]:
    rows = [row for row in _read(WANG_DATA) if row["fit_included"] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3"]) for row in rows])
    fit = least_squares(
        lambda parameters: (
            bm3_pressure(volume, parameters[0], parameters[1], 4.8) - pressure
        ),
        np.array([45.58, 232.0]),
    )
    fixed_v0 = least_squares(
        lambda parameters: bm3_pressure(volume, 45.58, parameters[0], 4.8) - pressure,
        np.array([232.0]),
    )
    return {
        "observations": len(rows),
        "V0": float(fit.x[0]),
        "K0": float(fit.x[1]),
        "K0_prime": 4.8,
        "pressure_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        "fixed_V0_K0": float(fixed_v0.x[0]),
    }


def _fit_mao_reanalyses() -> dict[str, dict[str, float]]:
    rows = _read(MAO_DATA)
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["unit_cell_volume_a3"]) for row in rows])
    sigma = np.array(
        [
            float(row["pressure_standard_deviation_gpa"])
            if row["pressure_standard_deviation_gpa"]
            else np.nan
            for row in rows
        ]
    )

    def fit(mask: np.ndarray, weights: np.ndarray) -> dict[str, float]:
        result = least_squares(
            lambda parameters: (
                (bm3_pressure(volume[mask], *parameters) - pressure[mask])
                * weights[mask]
            ),
            np.array([45.5, 260.0, 4.2]),
        )
        residual = bm3_pressure(volume[mask], *result.x) - pressure[mask]
        return {
            "observations": int(mask.sum()),
            "V0": float(result.x[0]),
            "K0": float(result.x[1]),
            "K0_prime": float(result.x[2]),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        }

    finite = np.isfinite(sigma)
    safe_weights = np.where(finite, 1.0 / sigma, 1.0)
    return {
        "equal_weight_excluding_below_1_gpa": fit(
            pressure >= 1.0, np.ones_like(pressure)
        ),
        "pressure_weighted_all_with_reported_sigma": fit(finite, safe_weights),
        "pressure_weighted_excluding_below_1_gpa": fit(
            finite & (pressure >= 1.0), safe_weights
        ),
    }


def reproduce() -> dict[str, object]:
    oganov_rows = _read(OGANOV_DATA)
    volume = np.array([float(row["volume_a3"]) for row in oganov_rows])
    uncorrected = np.array([float(row["static_pressure_gpa"]) for row in oganov_rows])
    corrected = np.array([float(row["pressure_corrected_gpa"]) for row in oganov_rows])
    paw_small = bm3_pressure(volume, *OGANOV["paw_small_core_static"])
    corrected_static = bm3_pressure(volume, *OGANOV["pressure_corrected_static"])
    return {
        "oganov": {
            "record_count": len(OGANOV),
            "table3_observations": len(oganov_rows),
            "paw_small_core_table3_pressure_rmse_gpa": float(
                np.sqrt(np.mean((paw_small - uncorrected) ** 2))
            ),
            "pressure_corrected_static_table3_rmse_gpa": float(
                np.sqrt(np.mean((corrected_static - corrected) ** 2))
            ),
            "curve_checkpoints_at_v_over_v0_0_9_gpa": {
                key: float(bm3_pressure(0.9 * values[0], *values))
                for key, values in OGANOV.items()
            },
        },
        "wang": {
            "room_temperature": _fit_wang_room_temperature(),
            "mao_1989_reanalyses": _fit_mao_reanalyses(),
        },
        "chizmeshya": {
            "record_count": len(CHIZMESHYA),
            "curve_checkpoints_at_v_over_v0_0_9_gpa": {
                key: float(bm3_pressure(0.9 * values[0], *values))
                for key, values in CHIZMESHYA.items()
            },
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
