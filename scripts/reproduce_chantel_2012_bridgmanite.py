#!/usr/bin/env python3
"""Audit Chantel et al. (2012) density and acoustic finite-strain fits."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = (
    ROOT
    / "peritheos/data/datasets/bridgmanite-chantel-2012-table1-density-velocity.csv"
)


def bm3_pressure(
    rho: np.ndarray, rho0: float, k0: float, k0_prime: float
) -> np.ndarray:
    """Evaluate BM3 using density ratio rho/rho0 = V0/V."""
    eta = (rho / rho0) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def reproduce() -> dict[str, object]:
    """Return curve checks and the independently reconstructed acoustic fit."""
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    room = [row for row in rows if row["room_temperature_fit_included"] == "1"]
    rho = np.array([float(row["density_g_cm3"]) for row in room])
    pressure = np.array([float(row["pressure_gpa"]) for row in room])
    predicted = bm3_pressure(rho, 100.387 / 24.42, 252.0, 4.1)

    acoustic = [row for row in room if row["vp_km_s"]]
    arho = np.array([float(row["density_g_cm3"]) for row in acoustic])
    vp = np.array([float(row["vp_km_s"]) for row in acoustic])
    vs = np.array([float(row["vs_km_s"]) for row in acoustic])
    strain = 0.5 * ((arho / 4.110) ** (2.0 / 3.0) - 1.0)
    observed_k = arho * (vp**2 - 4.0 * vs**2 / 3.0)

    def predicted_k(free: np.ndarray) -> np.ndarray:
        k0, k0_prime = free
        return (1.0 + 2.0 * strain) ** 2.5 * (
            k0 + (3.0 * k0 * k0_prime - 5.0 * k0) * strain
        )

    fit = least_squares(lambda free: predicted_k(free) - observed_k, [247.0, 4.5])
    residual = predicted - pressure
    return {
        "observations": len(rows),
        "room_temperature_observations": len(room),
        "accepted_table3_model": {
            "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "maximum_absolute_pressure_residual_gpa": float(np.max(np.abs(residual))),
        },
        "table2_acoustic_reconstruction": {
            "K_S": float(fit.x[0]),
            "K_S_prime": float(fit.x[1]),
            "published_K_S": 247.0,
            "published_K_S_prime": 4.5,
        },
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
