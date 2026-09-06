"""Reproduce Sakai et al. (2025) Table S9 generalized R-S curves."""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from peritheos.eos.rt import RydbergStacey

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"
PARAMETERS = {
    "copper": (47.218, 133.6, 5.2035, 2.0),
    "rhenium": (29.468, 352.6, 4.411, 5.0 / 3.0),
    "platinum": (60.409, 274.1, 5.128, 5.0 / 3.0),
    "tungsten": (31.724, 296.0, 4.415, 5.0 / 3.0),
    "gold": (67.716, 167.0, 5.780, 5.0 / 3.0),
    "molybdenum": (31.12, 261.0, 4.141, 5.0 / 3.0),
    "mgo": (74.698, 160.6, 4.227, 5.0 / 3.0),
    "nacl_b2": (39.57, 37.12, 4.637, 5.0 / 3.0),
    "iron": (22.42, 161.9, 5.51, 5.0 / 3.0),
}


def reproduce() -> dict[str, dict[str, float]]:
    """Return pressure residual metrics against official derived Table S9."""
    metrics = {}
    for material, parameters in PARAMETERS.items():
        path = DATA / f"{material}-sakai-2025-table-s9-rs-eos-grid.csv"
        with path.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        volumes = np.array([float(row["volume_a3"]) for row in rows])
        observed = np.array([float(row["pressure_gpa"]) for row in rows])
        calculated = np.asarray(RydbergStacey(*parameters).pressure(volumes))
        residual = calculated - observed
        metrics[material] = {
            "rows": float(len(rows)),
            "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "max_abs_gpa": float(np.max(np.abs(residual))),
        }
    return metrics


if __name__ == "__main__":
    for material, values in reproduce().items():
        print(
            f"{material}: n={int(values['rows'])}, "
            f"RMSE={values['rmse_gpa']:.12g} GPa, "
            f"max={values['max_abs_gpa']:.12g} GPa"
        )
