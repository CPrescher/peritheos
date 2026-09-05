"""Reproduce the three Mao et al. (2011) Figure 3 BM3 curves.

The article provides no numerical pressure-volume table.  This audit therefore
uses only separable marker centers digitized from the publisher's Figure 3.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"
DOLOMITE = DATA / (
    "ca0988mg0918fe0078mn0016c2o6-mao-2011-figure3-fe-dolomite-digitized.csv"
)
DOLOMITE_III = DATA / (
    "ca0988mg0918fe0078mn0016c2o6-mao-2011-figure3-dolomite-iii-digitized.csv"
)

PUBLISHED = {
    "fe_dolomite": {"V0": 321.77, "K0": 94.1, "K0_prime": 4.0},
    "dolomite_iii_high_spin": {"V0": 239.2, "K0": 164.0, "K0_prime": 4.0},
    "dolomite_iii_low_spin": {"V0": 231.8, "K0": 184.0, "K0_prime": 4.0},
}


def birch_murnaghan_pressure(
    volume: np.ndarray | float, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    """Return third-order Birch-Murnaghan pressure in GPa."""
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def _load(path: Path, flag: str) -> tuple[np.ndarray, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream) if row[flag] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    return pressure, volume


def _volume_at_pressure(pressure: float, parameters: dict[str, float]) -> float:
    return brentq(
        lambda volume: float(
            birch_murnaghan_pressure(
                volume,
                parameters["V0"],
                parameters["K0"],
                parameters["K0_prime"],
            )
            - pressure
        ),
        0.5 * parameters["V0"],
        parameters["V0"],
    )


def _diagnose(
    pressure: np.ndarray, volume: np.ndarray, parameters: dict[str, float]
) -> dict[str, object]:
    calculated_pressure = birch_murnaghan_pressure(
        volume,
        parameters["V0"],
        parameters["K0"],
        parameters["K0_prime"],
    )
    pressure_residual = calculated_pressure - pressure
    calculated_volume = np.array(
        [_volume_at_pressure(value, parameters) for value in pressure]
    )
    volume_residual = calculated_volume - volume

    def pressure_objective(free: np.ndarray) -> np.ndarray:
        return birch_murnaghan_pressure(volume, free[0], free[1], 4.0) - pressure

    pressure_fit = least_squares(
        pressure_objective, [parameters["V0"], parameters["K0"]]
    )

    def volume_objective(free: np.ndarray) -> np.ndarray:
        trial = {"V0": float(free[0]), "K0": float(free[1]), "K0_prime": 4.0}
        return (
            np.array([_volume_at_pressure(value, trial) for value in pressure]) - volume
        )

    volume_fit = least_squares(volume_objective, [parameters["V0"], parameters["K0"]])
    return {
        "observations": int(pressure.size),
        "published_curve": {
            "pressure_rmse_gpa": float(np.sqrt(np.mean(pressure_residual**2))),
            "pressure_maximum_absolute_residual_gpa": float(
                np.max(np.abs(pressure_residual))
            ),
            "volume_rmse_a3": float(np.sqrt(np.mean(volume_residual**2))),
            "volume_maximum_absolute_residual_a3": float(
                np.max(np.abs(volume_residual))
            ),
        },
        "unweighted_pressure_refit": {
            "parameters": {
                "V0": float(pressure_fit.x[0]),
                "K0": float(pressure_fit.x[1]),
                "K0_prime": 4.0,
            },
            "pressure_rmse_gpa": float(np.sqrt(np.mean(pressure_fit.fun**2))),
        },
        "unweighted_volume_refit": {
            "parameters": {
                "V0": float(volume_fit.x[0]),
                "K0": float(volume_fit.x[1]),
                "K0_prime": 4.0,
            },
            "volume_rmse_a3": float(np.sqrt(np.mean(volume_fit.fun**2))),
        },
    }


def reproduce() -> dict[str, object]:
    """Return deterministic curve and unweighted-refit diagnostics."""
    dolomite = _load(DOLOMITE, "fit_included")
    high_spin = _load(DOLOMITE_III, "high_spin_fit_included")
    low_spin = _load(DOLOMITE_III, "low_spin_fit_included")
    return {
        "fe_dolomite": _diagnose(*dolomite, PUBLISHED["fe_dolomite"]),
        "dolomite_iii_high_spin": _diagnose(
            *high_spin, PUBLISHED["dolomite_iii_high_spin"]
        ),
        "dolomite_iii_low_spin": _diagnose(
            *low_spin, PUBLISHED["dolomite_iii_low_spin"]
        ),
    }


def main() -> None:
    """Print the reproduction as stable JSON."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
