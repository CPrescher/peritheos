"""Reproduce Hirose et al. (2005) MgGeO3 post-perovskite fits."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/mggeo3-post-perovskite-hirose-2005-table3-pv.csv"


def reproduce() -> dict[str, dict[str, object]]:
    """Return published-curve residuals and independent diagnostic refits."""
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    pressure_sigma = np.array([float(row["pressure_uncertainty_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3"]) for row in rows])

    cases = {
        "fixed": {
            "published": (183.1, 192.0, 4.0),
            "residual": lambda values: (
                BM3(values[0], values[1], 4.0).pressure(volume) - pressure
            ),
            "start": (183.1, 192.0),
        },
        "free": {
            "published": (182.2, 210.0, 3.5),
            "residual": lambda values: (
                (BM3(*values).pressure(volume) - pressure) / pressure_sigma
            ),
            "start": (182.2, 210.0, 3.5),
        },
    }
    metrics = {}
    for name, case in cases.items():
        published = BM3(*case["published"])
        published_residual = np.asarray(published.pressure(volume)) - pressure
        fit = least_squares(case["residual"], case["start"])
        fit_parameters = list(map(float, fit.x))
        if name == "fixed":
            fit_parameters.append(4.0)
        metrics[name] = {
            "rows": len(rows),
            "published_rmse_gpa": float(np.sqrt(np.mean(published_residual**2))),
            "published_max_abs_gpa": float(np.max(np.abs(published_residual))),
            "refit_parameters": tuple(fit_parameters),
        }
    return metrics


if __name__ == "__main__":
    for name, values in reproduce().items():
        print(
            f"{name}: n={values['rows']}, "
            f"RMSE={values['published_rmse_gpa']:.9f} GPa, "
            f"max={values['published_max_abs_gpa']:.9f} GPa, "
            f"refit={values['refit_parameters']}"
        )
