"""Reproduce Xiao et al. (2013) SrSiO3 EOS diagnostics."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets/srsio3-xiao-2013-table1-pv.csv"

PUBLISHED = {
    "experimental_cubic": (49.18, 211.0, 4.0),
    "gga_cubic": (49.97, 207.7, 4.0),
    "gga_6h": (311.23, 183.8, 4.0),
}


def bm3_pressure(volume, v0, k0, kp):
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def reproduce() -> dict[str, object]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    selected = [row for row in rows if row["fit_included"] == "1"]
    p = np.array([float(row["pressure_gpa"]) for row in selected])
    v = np.array([float(row["volume_a3"]) for row in selected])
    v0, k0, kp = PUBLISHED["experimental_cubic"]
    residual = bm3_pressure(v, v0, k0, kp) - p
    checkpoints = {
        name: float(bm3_pressure(0.9 * pars[0], *pars))
        for name, pars in PUBLISHED.items()
    }
    return {
        "observations": len(rows),
        "selected_observations": len(selected),
        "experimental_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "experimental_max_abs_residual_gpa": float(np.max(np.abs(residual))),
        "pressure_at_0.9_v0_gpa": checkpoints,
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
