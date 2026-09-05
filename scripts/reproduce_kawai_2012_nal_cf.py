#!/usr/bin/env python3
"""Reproduce Kawai and Tsuchiya's (2012) NAL and CF BM3 curves.

The paper reports coefficients in cm^3/mol on a four-oxygen basis, while
Figure 2c-d uses one NaMg2Al5SiO12 (12 oxygen) formula.  Peritheos stores
conventional-cell volumes: one 12-O formula for NAL and four O4 units for CF.
"""

from __future__ import annotations

import csv
import math
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "peritheos" / "data" / "datasets"
MOLAR_A3 = 0.602214076

CASES = {
    "nal": {
        "path": DATA_ROOT
        / "namg2al5sio12-nal-kawai-2012-figure2c-vector-digitized.csv",
        "published": (35.8 * 3.0 / MOLAR_A3, 217.7, 4.08),
    },
    "cf": {
        "path": DATA_ROOT / "namg2al5sio12-cf-kawai-2012-figure2d-vector-digitized.csv",
        "published": (35.2 * 4.0 / MOLAR_A3, 213.2, 4.12),
    },
}


def bm3_pressure(volume, v0, k0, k0_prime):
    """Standard Eulerian third-order Birch-Murnaghan pressure in GPa."""
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def bm3_volumes(pressure, parameters):
    """Invert BM3 independently at each pressure."""
    v0, k0, k0_prime = parameters
    return np.array(
        [
            brentq(
                lambda volume: bm3_pressure(volume, v0, k0, k0_prime) - value,
                0.3 * v0,
                1.2 * v0,
            )
            for value in np.asarray(pressure, dtype=float)
        ]
    )


def load_arrays(phase):
    """Load digitized pressure and conventional-cell volume arrays."""
    with CASES[phase]["path"].open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    return rows, pressure, volume


def reproduce_phase(phase):
    """Return published-curve metrics and a diagnostic volume-residual refit."""
    rows, pressure, volume = load_arrays(phase)
    published = np.array(CASES[phase]["published"], dtype=float)
    published_volume = bm3_volumes(pressure, published)
    published_residual = published_volume - volume
    result = least_squares(
        lambda beta: bm3_volumes(pressure, beta) - volume,
        x0=published,
        bounds=(
            (0.8 * published[0], 50.0, 1.0),
            (1.2 * published[0], 500.0, 10.0),
        ),
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
        max_nfev=5000,
    )
    refit_residual = bm3_volumes(pressure, result.x) - volume
    return {
        "rows": len(rows),
        "published_parameters": published,
        "published_curve_volume_rmse_a3": math.sqrt(
            float(np.mean(published_residual**2))
        ),
        "published_curve_max_abs_volume_residual_a3": float(
            np.max(np.abs(published_residual))
        ),
        "volume_residual_refit": result.x,
        "volume_residual_refit_rmse_a3": math.sqrt(float(np.mean(refit_residual**2))),
    }


def main():
    for phase in ("nal", "cf"):
        result = reproduce_phase(phase)
        print(f"{phase.upper()}: {result['rows']} digitized Figure 2 points")
        print(
            "  published conventional-cell V0, K0, K0': "
            + ", ".join(f"{value:.12g}" for value in result["published_parameters"])
        )
        print(
            "  published-curve volume RMSE / max: "
            f"{result['published_curve_volume_rmse_a3']:.12f} / "
            f"{result['published_curve_max_abs_volume_residual_a3']:.12f} A^3"
        )
        print(
            "  diagnostic unweighted volume-residual refit V0, K0, K0': "
            + ", ".join(f"{value:.12f}" for value in result["volume_residual_refit"])
        )
        print(
            f"  refit volume RMSE: {result['volume_residual_refit_rmse_a3']:.12f} A^3"
        )


if __name__ == "__main__":
    main()
