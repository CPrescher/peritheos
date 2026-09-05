#!/usr/bin/env python3
"""Reproduce Mao et al. (2011) 300 K ferropericlase branch fits.

Primary source: corrected Figure 1 and Results paragraph 7,
doi:10.1029/2011GL049915; correction doi:10.1029/2011GL050814.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "mg075fe025o-mao-2011-figure1-300k-digitized.csv"
)

PUBLISHED = {
    "high_spin": {"V0": 76.34, "K0": 162.0, "K0_prime": 4.0},
    "low_spin": {"V0": 74.4, "K0": 166.0, "K0_prime": 4.0},
}


def birch_murnaghan_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float = 4.0
) -> np.ndarray:
    """Evaluate the standard Eulerian third-order Birch-Murnaghan EOS."""
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def load_data() -> dict[str, np.ndarray]:
    """Load every digitized 300 K marker, including crossover exclusions."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "volume": np.array([float(row["cell_volume_a3"]) for row in rows]),
        "spin_regime": np.array([row["spin_regime"] for row in rows]),
        "high_spin": np.array([row["used_in_high_spin_fit"] == "1" for row in rows]),
        "low_spin": np.array([row["used_in_low_spin_fit"] == "1" for row in rows]),
    }


def reproduce_branch(branch: str) -> dict[str, object]:
    """Evaluate and independently refit one source-defined 300 K branch."""
    data = load_data()
    selected = data[branch]
    published = PUBLISHED[branch]
    predicted = birch_murnaghan_pressure(
        data["volume"][selected],
        published["V0"],
        published["K0"],
        published["K0_prime"],
    )
    published_residual = predicted - data["pressure"][selected]

    if branch == "high_spin":

        def residual(parameters: np.ndarray) -> np.ndarray:
            return (
                birch_murnaghan_pressure(
                    data["volume"][selected], published["V0"], parameters[0]
                )
                - data["pressure"][selected]
            )

        fit = least_squares(residual, np.array([published["K0"]]))
        refit_parameters = {
            "V0": published["V0"],
            "K0": float(fit.x[0]),
            "K0_prime": 4.0,
        }
        fixed_parameters = ["V0", "K0_prime"]
    else:

        def residual(parameters: np.ndarray) -> np.ndarray:
            return (
                birch_murnaghan_pressure(
                    data["volume"][selected], parameters[0], parameters[1]
                )
                - data["pressure"][selected]
            )

        fit = least_squares(residual, np.array([published["V0"], published["K0"]]))
        refit_parameters = {
            "V0": float(fit.x[0]),
            "K0": float(fit.x[1]),
            "K0_prime": 4.0,
        }
        fixed_parameters = ["K0_prime"]

    return {
        "observations": int(np.count_nonzero(selected)),
        "pressure_range_gpa": [
            float(np.min(data["pressure"][selected])),
            float(np.max(data["pressure"][selected])),
        ],
        "published_parameters": published,
        "published_curve_pressure_rmse_gpa": float(
            np.sqrt(np.mean(published_residual**2))
        ),
        "published_curve_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(published_residual))
        ),
        "fixed_parameters": fixed_parameters,
        "refit_parameters": refit_parameters,
        "refit_curve_pressure_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        "refit_curve_max_abs_pressure_residual_gpa": float(np.max(np.abs(fit.fun))),
    }


def main() -> None:
    """Print deterministic branch diagnostics as JSON."""
    print(
        json.dumps(
            {branch: reproduce_branch(branch) for branch in PUBLISHED},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
