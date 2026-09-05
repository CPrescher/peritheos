#!/usr/bin/env python3
"""Check the Speziale et al. (2007) high-spin ferropericlase EOS.

Primary source: Table 1 and Results section 3.2,
doi:10.1029/2006JB004730.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA_PATH = (
    ROOT / "peritheos" / "data" / "datasets" / "mg080fe020o-speziale-2007-table1-pv.csv"
)

PUBLISHED = {"V0": 76.03, "K0": 158.0, "K0_prime": 4.4}


def birch_murnaghan_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    """Evaluate the standard Eulerian third-order Birch-Murnaghan EOS."""
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def load_data() -> dict[str, np.ndarray]:
    """Load the complete source table and its explicit fit-selection flag."""
    with DATA_PATH.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.array([float(row["pressure_gpa"]) for row in rows]),
        "pressure_sigma": np.array(
            [
                float(row["pressure_sigma_gpa"]) if row["pressure_sigma_gpa"] else 0.0
                for row in rows
            ]
        ),
        "volume": np.array([float(row["cell_volume_a3"]) for row in rows]),
        "volume_sigma": np.array([float(row["cell_volume_sigma_a3"]) for row in rows]),
        "fit_included": np.array([row["fit_included"] == "1" for row in rows]),
    }


def reproduce() -> dict[str, object]:
    """Evaluate the publication curve and expose an unweighted refit audit."""
    data = load_data()
    selected = data["fit_included"]
    pressure = data["pressure"][selected]
    volume = data["volume"][selected]
    pressure_sigma = data["pressure_sigma"][selected]
    volume_sigma = data["volume_sigma"][selected]

    predicted = birch_murnaghan_pressure(
        volume,
        **{
            "v0": PUBLISHED["V0"],
            "k0": PUBLISHED["K0"],
            "k0_prime": PUBLISHED["K0_prime"],
        },
    )
    residual = predicted - pressure

    step = 1.0e-4
    derivative = (
        birch_murnaghan_pressure(
            volume + step,
            PUBLISHED["V0"],
            PUBLISHED["K0"],
            PUBLISHED["K0_prime"],
        )
        - birch_murnaghan_pressure(
            volume - step,
            PUBLISHED["V0"],
            PUBLISHED["K0"],
            PUBLISHED["K0_prime"],
        )
    ) / (2.0 * step)
    effective_sigma = np.sqrt(pressure_sigma**2 + (derivative * volume_sigma) ** 2)
    normalized_residual = residual / effective_sigma

    def unweighted_pressure_residual(parameters: np.ndarray) -> np.ndarray:
        return (
            birch_murnaghan_pressure(
                volume, parameters[0], parameters[1], parameters[2]
            )
            - pressure
        )

    refit = least_squares(
        unweighted_pressure_residual,
        np.array([PUBLISHED["V0"], PUBLISHED["K0"], PUBLISHED["K0_prime"]]),
    )

    return {
        "observations_in_source_table": int(data["pressure"].size),
        "observations_in_high_spin_fit": int(np.count_nonzero(selected)),
        "excluded_spin_transition_observations": int(np.count_nonzero(~selected)),
        "fit_pressure_range_gpa": [float(pressure.min()), float(pressure.max())],
        "published_parameters": PUBLISHED,
        "published_curve_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "published_curve_max_abs_pressure_residual_gpa": float(
            np.max(np.abs(residual))
        ),
        "published_curve_effective_sigma_rms": float(
            np.sqrt(np.mean(normalized_residual**2))
        ),
        "published_curve_max_abs_effective_sigma": float(
            np.max(np.abs(normalized_residual))
        ),
        "unweighted_pressure_space_refit": {
            "V0": float(refit.x[0]),
            "K0": float(refit.x[1]),
            "K0_prime": float(refit.x[2]),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(refit.fun**2))),
        },
    }


def main() -> None:
    """Print deterministic reproduction diagnostics as JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
