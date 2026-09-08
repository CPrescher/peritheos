#!/usr/bin/env python3
"""Independently reproduce Sueda (2009), Tables 1-2 and Eqs. (1)-(8).

The MGD implementation uses the dimensionally consistent interpretation of the
printed equations documented in the literature audit. No Peritheos evaluator is
used to calculate the reference curves or diagnostic least-squares fits.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.constants import Boltzmann
from scipy.integrate import quad
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = ROOT / "peritheos/data/datasets/mgal2o4-cafe2o4-sueda-2009-table1-pvt.csv"
STATIC = {"V0": 240.1, "K0": 205.0, "K0_prime": 4.1}
THERMAL = {
    "htbm": {"alpha0": 1.96e-5, "alpha1": 1.64e-8, "dK_dT": -0.030},
    "mgd": {"gamma0": 1.73, "q": 2.03, "theta0": 1546.0},
}
ERRORS = {
    "htbm": {"alpha0": 0.13e-5, "alpha1": 0.24e-8, "dK_dT": 0.002},
    "mgd": {"gamma0": 0.07, "q": 0.37, "theta0": 104.0},
}


def load_data():
    """Read every source observation; blank uncertainty fields remain missing."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def bm3(volume, V0=240.1, K0=205.0, K0_prime=4.1):
    """Equation (1), with conventional-cell volumes and pressure in GPa."""
    eta = (V0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * K0 * (eta**7 - eta**5) * (1.0 + 0.75 * (K0_prime - 4.0) * (eta**2 - 1.0))
    )


def htbm_pressure(volume, temperature, alpha0=1.96e-5, alpha1=1.64e-8, dK_dT=-0.030):
    """Equations (1)-(4): alpha0 is the T=0 intercept, c0=0."""
    temperature = np.asarray(temperature, dtype=float)
    v0_t = 240.1 * np.exp(
        alpha0 * (temperature - 300.0) + 0.5 * alpha1 * (temperature**2 - 300.0**2)
    )
    return bm3(volume, V0=v0_t, K0=205.0 + dK_dT * (temperature - 300.0))


def mgd_pressure(volume, temperature, gamma0=1.73, q=2.03, theta0=1546.0):
    """Equations (5)-(8), using 28 atoms/cell and joules per physical cell.

    The native library instead uses seven atoms per formula unit and molar
    formula-unit volume; the independent normalization here checks that adapter.
    """
    volume, temperature = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    gamma = gamma0 * (volume / 240.1) ** q
    theta = theta0 * np.exp((gamma0 - gamma) / q)

    def energy(t, characteristic):
        integral = quad(lambda x: x**3 / np.expm1(x), 0.0, characteristic / t)[0]
        return 9.0 * 28.0 * Boltzmann * t * (t / characteristic) ** 3 * integral

    delta_energy = np.array(
        [
            energy(t, characteristic) - energy(300.0, characteristic)
            for t, characteristic in zip(temperature.flat, theta.flat)
        ]
    ).reshape(volume.shape)
    # A^3 -> m^3 and Pa -> GPa together give 1e-21.
    return bm3(volume) + gamma * delta_energy / (volume * 1.0e-21)


def reproduce():
    """Run the source's staged protocol with explicit diagnostic weighting.

    The paper does not specify weights. These are unweighted pressure-residual
    diagnostics, with the printed static triplet fixed in each thermal stage.
    The 16-row static reconstruction includes all explicitly 300 K observations.
    """
    rows = load_data()
    volume, temperature, pressure = (
        np.array([float(row[key]) for row in rows])
        for key in ("volume_a3", "temperature_k", "pressure_gpa")
    )
    static_mask = temperature == 300.0
    static_fit = least_squares(
        lambda values: bm3(volume[static_mask], *values) - pressure[static_mask],
        list(STATIC.values()),
        x_scale=[240.0, 205.0, 4.0],
    )
    result = {
        "rows": len(rows),
        "static": {
            "observations": int(static_mask.sum()),
            "parameters": dict(zip(STATIC, static_fit.x.tolist())),
            "rmse_gpa": float(np.sqrt(np.mean(static_fit.fun**2))),
            "solver_success": bool(static_fit.success),
        },
    }
    for name, function in (("htbm", htbm_pressure), ("mgd", mgd_pressure)):
        published = THERMAL[name]
        fit = least_squares(
            lambda values: function(volume, temperature, *values) - pressure,
            list(published.values()),
            x_scale=np.abs(list(published.values())),
            bounds=(
                ([-1e-4, -1e-7, -0.09], [1e-4, 1e-7, 0.0])
                if name == "htbm"
                else ([0.1, 0.1, 100.0], [5.0, 10.0, 5000.0])
            ),
        )
        calculated = function(volume, temperature)
        residual = calculated - pressure
        result[name] = {
            "parameters": dict(zip(published, fit.x.tolist())),
            "parameter_differences_in_reported_errors": {
                key: float((value - published[key]) / ERRORS[name][key])
                for key, value in zip(published, fit.x)
            },
            "published_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "published_max_abs_residual_gpa": float(np.max(np.abs(residual))),
            "refit_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
            "solver_success": bool(fit.success),
            "benchmark_pressures_gpa": {
                str(rows[i]["source_row"]): float(calculated[i])
                for i in (0, 25, 30, 35)
            },
        }
    return result


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, allow_nan=False))
