#!/usr/bin/env python3
"""Reproduce the Fu et al. (2024) Fe10-Al14 bridgmanite EOS checks.

Primary source: American Mineralogist 109, 872-881,
doi:10.2138/am-2023-8969. The 22 observations combine direct transcriptions
of official Online Materials Table S1 with conventional-cell volumes from the
peer-reviewed deposited CIF. The article reports weighted EosFit7-GUI fits but
does not disclose whether the numerical inputs were this merged series or the
two separate platelet series plotted in Figure S2. The refits below are
therefore diagnostic and do not claim source-fit parity.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "mg088fe010al014si090o3-bridgmanite-fu-2024-table-s1-cif-pv.csv"
)

PUBLISHED = {
    "BM2": np.array([163.85, 242.0]),
    "BM3": np.array([164.64, 228.0, 4.1]),
}


def birch_murnaghan_pressure(
    volume: np.ndarray | float, parameters: np.ndarray
) -> np.ndarray:
    """Evaluate the standard Eulerian BM2 or BM3 pressure equation."""
    volume = np.asarray(volume, dtype=float)
    parameters = np.asarray(parameters, dtype=float)
    v0, k0 = parameters[:2]
    k0_prime = 4.0 if len(parameters) == 2 else parameters[2]
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def load_data() -> dict[str, np.ndarray]:
    """Load the direct Table S1/CIF transcription."""
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "pressure_sigma": np.array([float(row["pressure_sd_gpa"]) for row in rows]),
        "volume": np.array(
            [float(row["cif_volume_a3_conventional_cell"]) for row in rows]
        ),
    }


def diagnostic_refit(
    pressure: np.ndarray,
    pressure_sigma: np.ndarray,
    volume: np.ndarray,
    initial: np.ndarray,
    *,
    weighted: bool,
) -> tuple[np.ndarray, float, float]:
    """Fit pressure residuals, optionally divided by reported pressure sigma."""

    def residual(parameters: np.ndarray) -> np.ndarray:
        values = birch_murnaghan_pressure(volume, parameters) - pressure
        return values / pressure_sigma if weighted else values

    lower = [0.8 * initial[0], 1.0]
    upper = [1.2 * initial[0], 1000.0]
    if len(initial) == 3:
        lower.append(0.0)
        upper.append(10.0)
    result = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        x_scale="jac",
    )
    if not result.success:
        raise RuntimeError(result.message)
    raw_residual = birch_murnaghan_pressure(volume, result.x) - pressure
    degrees_of_freedom = len(pressure) - len(initial)
    reduced_chi_square = float(np.sum(residual(result.x) ** 2) / degrees_of_freedom)
    return (
        result.x,
        float(np.sqrt(np.mean(raw_residual**2))),
        reduced_chi_square,
    )


def reproduce_model(name: str) -> dict[str, object]:
    """Return published-curve and diagnostic-refit results for one EOS."""
    data = load_data()
    published = PUBLISHED[name]
    calculated = birch_murnaghan_pressure(data["volume"], published)
    residual = calculated - data["pressure"]
    unweighted, unweighted_rmse, _ = diagnostic_refit(
        data["pressure"],
        data["pressure_sigma"],
        data["volume"],
        published,
        weighted=False,
    )
    weighted, weighted_rmse, reduced_chi_square = diagnostic_refit(
        data["pressure"],
        data["pressure_sigma"],
        data["volume"],
        published,
        weighted=True,
    )
    return {
        "rows": int(len(data["pressure"])),
        "published_parameters": published.tolist(),
        "published_curve_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_curve_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(residual))
        ),
        "low_pressure_state": {
            "observed_pressure_gpa": float(data["pressure"][0]),
            "volume_a3": float(data["volume"][0]),
            "calculated_pressure_gpa": float(calculated[0]),
        },
        "high_pressure_state": {
            "observed_pressure_gpa": float(data["pressure"][-1]),
            "volume_a3": float(data["volume"][-1]),
            "calculated_pressure_gpa": float(calculated[-1]),
        },
        "unweighted_pressure_refit": {
            "parameters": unweighted.tolist(),
            "rmse_gpa": unweighted_rmse,
        },
        "pressure_sigma_weighted_refit": {
            "parameters": weighted.tolist(),
            "rmse_gpa": weighted_rmse,
            "reduced_chi_square": reduced_chi_square,
        },
    }


def main() -> None:
    """Print deterministic JSON diagnostics for both published models."""
    print(
        json.dumps(
            {name: reproduce_model(name) for name in PUBLISHED},
            indent=2,
            sort_keys=True,
        )
    )


if __name__ == "__main__":
    main()
