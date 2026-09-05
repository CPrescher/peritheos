"""Reproduce the Katsura et al. (2009) corrected-table BM3 checks."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT / "peritheos/data/datasets/bridgmanite-katsura-2009-corrected-table1-pvt.csv"
)


def bm3_pressure(volume_ratio: np.ndarray, k0: float, k0_prime: float) -> np.ndarray:
    eta = volume_ratio ** (-1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def main() -> None:
    with DATASET.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    ambient = [row for row in rows if float(row["temperature_k"]) <= 308.0]
    pressure = np.asarray([float(row["pressure_gpa"]) for row in ambient])
    pressure_error = np.asarray(
        [float(row["pressure_uncertainty_gpa"]) for row in ambient]
    )
    volume_ratio = np.asarray([float(row["bridgmanite_v_v0"]) for row in ambient])

    preferred_pressure = bm3_pressure(volume_ratio, 256.0, 3.8)
    preferred_conditional = least_squares(
        lambda value: bm3_pressure(volume_ratio, value[0], 3.8) - pressure,
        np.asarray([256.0]),
    )
    preferred_unconstrained = least_squares(
        lambda value: bm3_pressure(volume_ratio, value[0], value[1]) - pressure,
        np.asarray([256.0, 3.8]),
    )
    fixed_k0 = least_squares(
        lambda value: (
            (bm3_pressure(volume_ratio, 253.0, value[0]) - pressure) / pressure_error
        ),
        np.asarray([4.1]),
    )
    fixed_pressure = bm3_pressure(volume_ratio, 253.0, 4.1)

    result = {
        "dataset_rows": len(rows),
        "ambient_rows": len(ambient),
        "preferred": {
            "published_curve_pressure_rmse_gpa": float(
                np.sqrt(np.mean((preferred_pressure - pressure) ** 2))
            ),
            "published_curve_max_absolute_pressure_residual_gpa": float(
                np.max(np.abs(preferred_pressure - pressure))
            ),
            "conditional_K0_gpa_at_published_K0_prime": float(
                preferred_conditional.x[0]
            ),
            "unconstrained_rounded_table_K0_gpa": float(preferred_unconstrained.x[0]),
            "unconstrained_rounded_table_K0_prime": float(preferred_unconstrained.x[1]),
        },
        "K0_fixed_sensitivity": {
            "weighted_refit_K0_prime": float(fixed_k0.x[0]),
            "published_curve_pressure_rmse_gpa": float(
                np.sqrt(np.mean((fixed_pressure - pressure) ** 2))
            ),
        },
    }
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
