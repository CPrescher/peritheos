"""Reproduce Gong et al.'s shock-derived finite-strain EOS checks."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets/mg092fe008sio3-gong-2004-table1-shock.csv"
MOLAR_MASS_G_MOL = 102.9102
REFERENCE_DENSITY_G_CM3 = 4.19
FORMULA_UNITS = 4.0
AVOGADRO = 6.02214076e23
V0 = MOLAR_MASS_G_MOL / REFERENCE_DENSITY_G_CM3 / AVOGADRO * 1e24 * FORMULA_UNITS
K0 = 260.1
K0_PRIME = 4.18


def bm3_pressure(volume):
    eta = (V0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * K0 * (eta**7 - eta**5) * (1.0 + 0.75 * (K0_PRIME - 4.0) * (eta**2 - 1.0))
    )


def reproduce():
    with DATA.open(newline="") as handle:
        rows = list(csv.DictReader(handle))
    particle_velocity = np.array([float(row["particle_velocity_km_s"]) for row in rows])
    shock_velocity = np.array([float(row["shock_velocity_km_s"]) for row in rows])
    us_residuals = 3.76 + 1.48 * particle_velocity - shock_velocity
    return {
        "rows": len(rows),
        "derived_v0_a3_conventional_z4": V0,
        "us_up_rms_residual_km_s": float(np.sqrt(np.mean(us_residuals**2))),
        "pressure_at_0_9_v0_gpa": float(bm3_pressure(0.9 * V0)),
        "pressure_at_0_8_v0_gpa": float(bm3_pressure(0.8 * V0)),
    }


def main():
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
