#!/usr/bin/env python3
"""Reproduce the Qin et al. (2023) CF-NaAlSiO4 BM3 checks.

Primary source: American Mineralogist 108, 2331-2337,
doi:10.2138/am-2022-8432. Observations are direct transcriptions of official
deposit AM-23-128432, Tables S3-S4. The publication reports error-weighted
EosFit7c regression but does not publish pressure uncertainties or the exact
weight configuration, so the refits below are explicitly diagnostic.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "peritheos" / "data" / "datasets"

SAMPLES = {
    "fe_free": {
        "dataset": "na093al102si100o4-qin-2023-table-s3-pv.csv",
        "published": np.array([241.6, 220.0, 2.6]),
    },
    "fe_bearing": {
        "dataset": "na088al099fe013si094o4-qin-2023-table-s4-pv.csv",
        "published": np.array([244.2, 211.0, 2.6]),
    },
}


def bm3_pressure(volume: np.ndarray | float, parameters: np.ndarray) -> np.ndarray:
    """Evaluate standard Eulerian third-order Birch-Murnaghan pressure."""
    volume = np.asarray(volume, dtype=float)
    v0, k0, k0_prime = np.asarray(parameters, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def bm3_volume(pressure: float, parameters: np.ndarray) -> float:
    """Invert the compression branch for a non-negative pressure."""
    v0 = float(parameters[0])
    if pressure == 0.0:
        return v0
    return brentq(
        lambda volume: float(bm3_pressure(volume, parameters) - pressure),
        0.5 * v0,
        v0,
    )


def load_data(filename: str) -> dict[str, np.ndarray]:
    """Load one directly transcribed source table."""
    with (DATA_ROOT / filename).open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "volume": np.array([float(row["volume_a3_conventional_cell"]) for row in rows]),
        "volume_sigma": np.array([float(row["volume_uncertainty_a3"]) for row in rows]),
    }


def unweighted_pressure_refit(
    pressure: np.ndarray, volume: np.ndarray, initial: np.ndarray
) -> tuple[np.ndarray, float]:
    """Fit all three BM3 parameters to unweighted pressure residuals."""
    result = least_squares(
        lambda parameters: bm3_pressure(volume, parameters) - pressure,
        initial,
        bounds=([0.8 * initial[0], 1.0, 0.0], [1.2 * initial[0], 1000.0, 10.0]),
        x_scale="jac",
    )
    if not result.success:
        raise RuntimeError(result.message)
    return result.x, float(np.sqrt(np.mean(result.fun**2)))


def volume_uncertainty_weighted_refit(
    pressure: np.ndarray,
    volume: np.ndarray,
    volume_sigma: np.ndarray,
    initial: np.ndarray,
) -> tuple[np.ndarray, float]:
    """Fit volume residuals divided by the printed volume uncertainties."""

    def residual(parameters: np.ndarray) -> np.ndarray:
        calculated = np.array(
            [bm3_volume(float(value), parameters) for value in pressure]
        )
        return (calculated - volume) / volume_sigma

    result = least_squares(
        residual,
        initial,
        bounds=([0.8 * initial[0], 1.0, 0.0], [1.2 * initial[0], 1000.0, 10.0]),
        x_scale="jac",
    )
    if not result.success:
        raise RuntimeError(result.message)
    degrees_of_freedom = len(pressure) - 3
    reduced_chi_square = float(np.sum(result.fun**2) / degrees_of_freedom)
    return result.x, reduced_chi_square


def reproduce_sample(name: str) -> dict[str, object]:
    """Return published-curve and diagnostic-refit results for one sample."""
    configuration = SAMPLES[name]
    data = load_data(str(configuration["dataset"]))
    published = np.asarray(configuration["published"], dtype=float)
    pressure_residual = bm3_pressure(data["volume"], published) - data["pressure"]
    unweighted, unweighted_rmse = unweighted_pressure_refit(
        data["pressure"], data["volume"], published
    )
    volume_weighted, reduced_chi_square = volume_uncertainty_weighted_refit(
        data["pressure"], data["volume"], data["volume_sigma"], published
    )
    return {
        "rows": int(len(data["pressure"])),
        "published_parameters": published.tolist(),
        "published_curve_pressure_rmse_gpa": float(
            np.sqrt(np.mean(pressure_residual**2))
        ),
        "published_curve_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(pressure_residual))
        ),
        "high_pressure_state": {
            "observed_pressure_gpa": float(data["pressure"][-1]),
            "volume_a3": float(data["volume"][-1]),
            "calculated_pressure_gpa": float(
                bm3_pressure(data["volume"][-1], published)
            ),
        },
        "unweighted_pressure_refit": {
            "parameters": unweighted.tolist(),
            "rmse_gpa": unweighted_rmse,
        },
        "volume_uncertainty_weighted_refit": {
            "parameters": volume_weighted.tolist(),
            "reduced_chi_square": reduced_chi_square,
        },
    }


def main() -> None:
    """Print deterministic JSON diagnostics for both source tables."""
    print(
        json.dumps(
            {name: reproduce_sample(name) for name in SAMPLES},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
