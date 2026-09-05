#!/usr/bin/env python3
"""Reproduce Funamori et al.'s (1998) two MgAl2O4 BM2 fits.

Each source EOS is constrained by a recovered ambient cell and one in-situ
high-pressure cell.  The script multiplies the printed high-pressure lattice
axes, evaluates the published curve, and independently solves for K0 with V0
fixed to the separately printed ambient volume.  K0'=4 is implicit in BM2.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"
CASES = {
    "cafe2o4_type": {
        "path": DATA / "mgal2o4-cafe2o4-funamori-1998-text-pv.csv",
        "V0": 240.3,
        "K0": 211.0,
    },
    "cati2o4_type": {
        "path": DATA / "mgal2o4-cati2o4-funamori-1998-text-pv.csv",
        "V0": 240.3,
        "K0": 206.0,
    },
}


def bm2_pressure(volume: np.ndarray | float, v0: float, k0: float) -> np.ndarray:
    """Return second-order Birch-Murnaghan pressure in GPa."""
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5)


def reproduce() -> dict[str, dict[str, float]]:
    """Return deterministic lattice, curve, and fixed-V0 refit diagnostics."""
    result = {}
    for phase, case in CASES.items():
        with case["path"].open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        high = rows[1]
        lattice_volume = (
            float(high["a_angstrom"])
            * float(high["b_angstrom"])
            * float(high["c_angstrom"])
        )
        pressure = float(high["pressure_gpa"])
        unit_curve_pressure = float(bm2_pressure(lattice_volume, case["V0"], 1.0))
        refit_k0 = pressure / unit_curve_pressure
        published_pressure = float(bm2_pressure(lattice_volume, case["V0"], case["K0"]))
        result[phase] = {
            "lattice_product_volume_a3": lattice_volume,
            "source_pressure_gpa": pressure,
            "published_curve_pressure_gpa": published_pressure,
            "published_curve_pressure_residual_gpa": published_pressure - pressure,
            "fixed_v0_refit_k0_gpa": refit_k0,
        }
    return result


def main() -> None:
    """Print the reproduction as stable JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
