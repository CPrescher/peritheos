#!/usr/bin/env python3
"""Reproduce the Kubo et al. (2000) aluminous bridgmanite EOS fits.

Primary source: Proceedings of the Japan Academy, Series B 76, 103-107,
doi:10.2183/pjab.76.103, Table I and the EOS discussion on page 106.
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
    / "mg09al02si09o3-bridgmanite-kubo-2000-table1-compression.csv"
)

PUBLISHED = {
    "fixed": np.array([225.5, 4.0]),
    "free": np.array([215.4, 7.2]),
}


def birch_murnaghan_pressure(
    volume: np.ndarray, reference_volume: np.ndarray, k0: float, k0_prime: float
) -> np.ndarray:
    """Evaluate the standard Eulerian third-order Birch-Murnaghan EOS."""
    eta = (reference_volume / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def load_data() -> dict[str, np.ndarray]:
    """Load the direct transcription of Table I."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_ruby_gpa"]) for row in rows]),
        "volume": np.array([float(row["volume_a3_conventional_cell"]) for row in rows]),
        "axis_product_volume": np.array(
            [
                float(row["lattice_a_angstrom"])
                * float(row["lattice_b_angstrom"])
                * float(row["lattice_c_angstrom"])
                for row in rows
            ]
        ),
        "run_reference_volume": np.array(
            [float(row["run_reference_volume_a3"]) for row in rows]
        ),
        "series": np.array([row["series"] for row in rows]),
    }


def _reference_from_ambient(
    values: np.ndarray, series: np.ndarray, pressure: np.ndarray
) -> np.ndarray:
    """Return the printed ambient value for each run series."""
    return np.array(
        [values[(series == name) & (pressure == 0.0)][0] for name in series]
    )


def reproduce_model(name: str) -> dict[str, object]:
    """Evaluate one published curve against all rounded Table I rows."""
    data = load_data()
    k0, k0_prime = PUBLISHED[name]
    calculated = birch_murnaghan_pressure(
        data["volume"], data["run_reference_volume"], k0, k0_prime
    )
    residual = calculated - data["pressure"]
    canonical = birch_murnaghan_pressure(
        data["volume"], np.full_like(data["volume"], 163.6), k0, k0_prime
    )
    canonical_residual = canonical - data["pressure"]
    return {
        "published_parameters": {"K0": k0, "K0_prime": k0_prime},
        "observations": int(data["pressure"].size),
        "series_normalized_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "series_normalized_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(residual))
        ),
        "single_163_6_a3_reference_pressure_rmse_gpa": float(
            np.sqrt(np.mean(canonical_residual**2))
        ),
        "single_163_6_a3_reference_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(canonical_residual))
        ),
    }


def refit_model(name: str, *, use_axis_products: bool = False) -> dict[str, object]:
    """Refit the eight finite-pressure rows with each series normalized at P=0."""
    data = load_data()
    volume = data["axis_product_volume"] if use_axis_products else data["volume"]
    reference_volume = _reference_from_ambient(volume, data["series"], data["pressure"])
    selected = data["pressure"] > 0.0
    fixed = name == "fixed"

    def residual(parameters: np.ndarray) -> np.ndarray:
        k0 = parameters[0]
        k0_prime = 4.0 if fixed else parameters[1]
        return (
            birch_murnaghan_pressure(
                volume[selected], reference_volume[selected], k0, k0_prime
            )
            - data["pressure"][selected]
        )

    initial = PUBLISHED[name][:1] if fixed else PUBLISHED[name]
    fit = least_squares(
        residual,
        initial,
        xtol=1.0e-14,
        ftol=1.0e-14,
        gtol=1.0e-14,
    )
    degrees_of_freedom = int(np.count_nonzero(selected) - fit.x.size)
    variance = float(np.sum(fit.fun**2) / degrees_of_freedom)
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * variance
    errors = np.sqrt(np.diag(covariance))
    parameters = {"K0": float(fit.x[0])}
    standard_errors = {"K0": float(errors[0])}
    if fixed:
        parameters["K0_prime"] = 4.0
    else:
        parameters["K0_prime"] = float(fit.x[1])
        standard_errors["K0_prime"] = float(errors[1])
    return {
        "input_volumes": "axis products" if use_axis_products else "printed volumes",
        "objective": "unweighted pressure residuals",
        "observations": int(np.count_nonzero(selected)),
        "parameters": parameters,
        "scaled_standard_errors": standard_errors,
        "pressure_rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
    }


def main() -> None:
    """Print deterministic curve and refit diagnostics as JSON."""
    result = {
        name: {
            "published_curve": reproduce_model(name),
            "printed_volume_refit": refit_model(name),
            "axis_product_refit": refit_model(name, use_axis_products=True),
        }
        for name in PUBLISHED
    }
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
