"""Reproduce Martin et al.'s CaIrO3 second-order BM fit."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets/cairo3-martin-2007-table1a-compression.csv"
V0 = 226.632
K0 = 180.2


def bm2_pressure(volume):
    eta = (V0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return 1.5 * K0 * (eta**7 - eta**5)


def reproduce():
    with DATA.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    observed = np.array([float(row["pressure_gpa"]) for row in rows])
    volumes = np.array([float(row["unit_cell_volume_a3"]) for row in rows])
    residuals = bm2_pressure(volumes) - observed
    return {
        "rows": len(rows),
        "rms_pressure_residual_gpa": float(np.sqrt(np.mean(residuals**2))),
        "max_abs_pressure_residual_gpa": float(np.max(np.abs(residuals))),
        "pressure_at_0_9_v0_gpa": float(bm2_pressure(0.9 * V0)),
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
