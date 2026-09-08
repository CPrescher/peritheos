"""Constrained diagnostics for Irifune et al. (2002), p. 652 and Figure 6.

Only K0 is optimized. The pressure objective has an analytic solution; the
normalized-volume objective uses an independent scalar BM2 inversion. Neither
objective claims to recover unspecified source weights or parameter covariance.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, minimize_scalar

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/mgal2o4-cafe2o4-irifune-2002-text-pv.csv"


def reproduce():
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream) if row["fit_included"] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3"]) for row in rows])
    ratio = volume / 240.6

    def basis(x):
        return 1.5 * (x ** (-7.0 / 3.0) - x ** (-5.0 / 3.0))

    factor = basis(ratio)
    pressure_fit = float(factor @ pressure / (factor @ factor))

    def volume_objective(k0):
        predicted = np.array(
            [brentq(lambda x: k0 * basis(x) - p, 0.5, 1.0) for p in pressure]
        )
        return float(np.sum((predicted - ratio) ** 2))

    volume_fit = minimize_scalar(
        volume_objective, bounds=(150.0, 280.0), method="bounded"
    )
    residual = 213.0 * factor - pressure
    return {
        "compressed_observations": len(rows),
        "free_parameters": ["K0"],
        "fixed_v0_pressure_refit_k0_gpa": pressure_fit,
        "fixed_v0_normalized_volume_refit_k0_gpa": float(volume_fit.x),
        "published_pressures_gpa": (213.0 * factor).tolist(),
        "published_pressure_residuals_gpa": residual.tolist(),
        "published_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "pressure_refit_rmse_gpa": float(
            np.sqrt(np.mean((pressure_fit * factor - pressure) ** 2))
        ),
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, allow_nan=False))
