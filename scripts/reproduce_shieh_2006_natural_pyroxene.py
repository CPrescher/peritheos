#!/usr/bin/env python3
"""Reproduce the Shieh et al. (2006) pv and ppv EOS diagnostics.

The primary pressure-volume observations occur only as marker centers in
Figure 2 of doi:10.1073/pnas.0506811103. The bundled CSV files retain the
digitization uncertainty separately from unavailable source fit weights.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_ROOT = ROOT / "peritheos" / "data" / "datasets"
PPV_DATA = DATA_ROOT / "mg090fe009al0005ca0005sio3-shieh-2006-figure2-ppv-digitized.csv"
PV_DATA = DATA_ROOT / "mg090fe009al0005ca0005sio3-shieh-2006-figure2-pv-digitized.csv"

PUBLISHED = {
    "ppv_preferred": {"V0": 164.9, "K0": 219.0, "K0_prime": 4.0},
    "ppv_sensitivity": {"V0": 166.2, "K0": 198.0, "K0_prime": 4.4},
    "pv": {"V0": 163.3, "K0": 255.0, "K0_prime": 3.7},
}


def birch_murnaghan_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    """Evaluate the standard Eulerian BM3 equation; K0'=4 is BM2."""
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def load_data(path: Path) -> dict[str, np.ndarray]:
    """Load all Figure 2 points flagged for the source fit."""
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream) if row["fit_included"] == "1"]
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "volume": np.array([float(row["volume_a3_conventional_cell"]) for row in rows]),
    }


def _curve_diagnostics(
    data: dict[str, np.ndarray], parameters: dict[str, float]
) -> dict[str, float]:
    calculated = birch_murnaghan_pressure(
        data["volume"],
        **{
            "v0": parameters["V0"],
            "k0": parameters["K0"],
            "k0_prime": parameters["K0_prime"],
        },
    )
    residuals = calculated - data["pressure"]
    return {
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "maximum_absolute_pressure_residual_gpa": float(np.max(np.abs(residuals))),
    }


def _refit_ppv(
    data: dict[str, np.ndarray], parameters: dict[str, float]
) -> dict[str, object]:
    k0_prime = parameters["K0_prime"]

    def residuals(free: np.ndarray) -> np.ndarray:
        return (
            birch_murnaghan_pressure(data["volume"], free[0], free[1], k0_prime)
            - data["pressure"]
        )

    result = least_squares(residuals, [parameters["V0"], parameters["K0"]])
    return {
        "objective": "unweighted pressure residuals",
        "fixed_parameters": {"K0_prime": k0_prime},
        "parameters": {
            "V0": float(result.x[0]),
            "K0": float(result.x[1]),
            "K0_prime": k0_prime,
        },
        "pressure_rmse_gpa": float(np.sqrt(np.mean(result.fun**2))),
    }


def _refit_pv(
    data: dict[str, np.ndarray], parameters: dict[str, float]
) -> dict[str, object]:
    v0 = parameters["V0"]
    k0_prime = parameters["K0_prime"]

    def residuals(free: np.ndarray) -> np.ndarray:
        return (
            birch_murnaghan_pressure(data["volume"], v0, free[0], k0_prime)
            - data["pressure"]
        )

    result = least_squares(residuals, [parameters["K0"]])
    return {
        "objective": "unweighted pressure residuals",
        "fixed_parameters": {"V0": v0, "K0_prime": k0_prime},
        "parameters": {
            "V0": v0,
            "K0": float(result.x[0]),
            "K0_prime": k0_prime,
        },
        "pressure_rmse_gpa": float(np.sqrt(np.mean(result.fun**2))),
    }


def reproduce() -> dict[str, object]:
    """Return deterministic published-curve and refit diagnostics."""
    ppv = load_data(PPV_DATA)
    pv = load_data(PV_DATA)
    result: dict[str, object] = {}
    for name in ("ppv_preferred", "ppv_sensitivity"):
        parameters = PUBLISHED[name]
        result[name] = {
            "observations": int(ppv["pressure"].size),
            "published_curve": _curve_diagnostics(ppv, parameters),
            "independent_refit": _refit_ppv(ppv, parameters),
        }
    parameters = PUBLISHED["pv"]
    result["pv"] = {
        "observations": int(pv["pressure"].size),
        "published_curve": _curve_diagnostics(pv, parameters),
        "independent_refit": _refit_pv(pv, parameters),
    }
    return result


def main() -> None:
    """Print the reproduction as stable JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
