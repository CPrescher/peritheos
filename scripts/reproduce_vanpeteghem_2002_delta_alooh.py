"""Evaluate both published delta-AlOOH BM3 fits against Table 1."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets/delta-alooh-vanpeteghem-2002-table1-pv.csv"

PUBLISHED = {
    "fixed_kp": {"V0": 56.54, "K0": 252.0, "K0_prime": 4.0},
    "free_kp": {"V0": 56.54, "K0": 228.0, "K0_prime": 7.0},
}


def bm3_pressure(volume, v0, k0, kp):
    """Evaluate BM3 independently of the library implementation."""
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def reproduce() -> dict[str, object]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3"]) for row in rows])
    results = {}
    for name, pars in PUBLISHED.items():
        calculated = bm3_pressure(volume, pars["V0"], pars["K0"], pars["K0_prime"])
        residual = calculated - pressure
        results[name] = {
            "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "maximum_absolute_pressure_residual_gpa": float(np.max(np.abs(residual))),
            "pressure_at_minimum_volume_gpa": float(calculated[0]),
            "zero_pressure_at_v0_gpa": float(
                bm3_pressure(pars["V0"], pars["V0"], pars["K0"], pars["K0_prime"])
            ),
        }
    return {"observations": len(rows), "results": results}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
