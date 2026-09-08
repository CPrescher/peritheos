#!/usr/bin/env python3
"""Independent BM2 diagnostics for Ono et al. (2006), Table 3.

No Peritheos evaluator is used. Both coefficients are free; K0-prime=4 is
implicit. Unweighted pressure residuals avoid inventing row pressure errors
or treating zero-rounded volume errors as exact observations. The source does
not specify its objective, weights, or treatment of the recovered cells.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/mgal2o4-cati2o4-ono-2006-table3.csv"


def bm2_pressure(volume, v0, k0):
    """Evaluate the printed p. 203 equation with K0-prime=4."""
    x = (np.asarray(volume) / v0) ** (1.0 / 3.0)
    return 1.5 * k0 * (x**-7 - x**-5)


def reproduce():
    """Return all-row and high-pressure-only sensitivity diagnostics."""
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    volume = np.array([float(row["volume_a3"]) for row in rows])
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    residual = bm2_pressure(volume, 238.9, 219.0) - pressure
    result = {
        "published_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_max_abs_residual_gpa": float(np.max(np.abs(residual))),
        "high_pressure_checkpoint_gpa": float(bm2_pressure(186.4, 238.9, 219.0)),
    }
    for name, mask in (
        ("all_rows", pressure >= 0),
        ("high_pressure_only", pressure > 0),
    ):
        fit = least_squares(
            lambda pars: bm2_pressure(volume[mask], *pars) - pressure[mask],
            [238.9, 219.0],
            bounds=([200.0, 100.0], [270.0, 400.0]),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        if not fit.success:
            raise RuntimeError(fit.message)
        result[name] = {
            "rows": int(mask.sum()),
            "V0": float(fit.x[0]),
            "K0": float(fit.x[1]),
            "rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        }
    return result


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
