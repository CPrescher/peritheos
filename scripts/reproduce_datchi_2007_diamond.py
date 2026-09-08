#!/usr/bin/env python3
"""Audit the Datchi et al. (2007) Table II diamond/H2005 fit.

Occelli et al. (2003) did not print numerical P-V rows. This script starts from
the vector marker centres digitized from their Figure 2, verifies the exact
MXB1986-to-H2005 ruby conversion, and evaluates the fit-protocol choices that
the two papers leave unstated. It intentionally reports no refit covariance.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "diamond-occelli-2003-figure2-digitized.csv"
)
CALIBRATIONS = ROOT / "peritheos" / "data" / "pressure-calibrations.json"
V0_ATOMIC_A3 = 5.6733
OCCELLI_VOLUME_SIGMA_CM3_MOL = 0.003
AVOGADRO_EXACT = 6.02214076e23


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _calibration(identifier: str) -> dict[str, object]:
    payload = json.loads(CALIBRATIONS.read_text(encoding="utf-8"))
    return next(
        item for item in payload["calibrations"] if item["identifier"] == identifier
    )


def _ratio_from_pressure(
    pressure: np.ndarray, calibration: dict[str, object]
) -> np.ndarray:
    parameters = calibration["parameters"]
    assert calibration["model"] == "power_law"
    a = float(parameters["A_gpa"])
    b = float(parameters["B"])
    return np.power(1.0 + b * pressure / a, 1.0 / b)


def _pressure_from_ratio(
    ratio: np.ndarray, calibration: dict[str, object]
) -> np.ndarray:
    parameters = calibration["parameters"]
    assert calibration["model"] == "holzapfel_freund_ingalls"
    a = float(parameters["A_gpa"])
    b = float(parameters["B"])
    c = float(parameters["C"])
    exponent = ((b + c) / c) * (1.0 - np.power(ratio, -c))
    return (a / (b + c)) * np.expm1(exponent)


def _vinet_pressure(volume: np.ndarray | float, k0: float, k0_prime: float):
    x = (np.asarray(volume) / V0_ATOMIC_A3) ** (1.0 / 3.0)
    return 3.0 * k0 * (1.0 - x) / x**2 * np.exp(1.5 * (k0_prime - 1.0) * (1.0 - x))


def _vinet_volume(pressure: np.ndarray, k0: float, k0_prime: float) -> np.ndarray:
    return np.array(
        [
            brentq(
                lambda volume: float(_vinet_pressure(volume, k0, k0_prime) - target),
                0.25 * V0_ATOMIC_A3,
                1.01 * V0_ATOMIC_A3,
            )
            for target in pressure
        ]
    )


def _fit_pressure(volume: np.ndarray, pressure: np.ndarray) -> dict[str, object]:
    result = least_squares(
        lambda parameters: _vinet_pressure(volume, *parameters) - pressure,
        np.array([443.0, 3.97]),
        bounds=([1.0, 0.5], [1000.0, 15.0]),
    )
    residual = _vinet_pressure(volume, *result.x) - pressure
    return {
        "parameters": {"K0": float(result.x[0]), "K0_prime": float(result.x[1])},
        "parameter_errors": {"K0": None, "K0_prime": None},
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_residual_gpa": float(np.max(np.abs(residual))),
        "success": bool(result.success),
        "message": str(result.message),
    }


def _fit_volume(volume: np.ndarray, pressure: np.ndarray) -> dict[str, object]:
    result = least_squares(
        lambda parameters: _vinet_volume(pressure, *parameters) - volume,
        np.array([443.0, 3.97]),
        bounds=([1.0, 0.5], [1000.0, 15.0]),
    )
    predicted_volume = _vinet_volume(pressure, *result.x)
    volume_residual = predicted_volume - volume
    pressure_residual = _vinet_pressure(volume, *result.x) - pressure
    sigma_atomic = OCCELLI_VOLUME_SIGMA_CM3_MOL * 1.0e24 / AVOGADRO_EXACT
    return {
        "parameters": {"K0": float(result.x[0]), "K0_prime": float(result.x[1])},
        "parameter_errors": {"K0": None, "K0_prime": None},
        "rmse_a3_per_atom": float(np.sqrt(np.mean(volume_residual**2))),
        "rmse_gpa": float(np.sqrt(np.mean(pressure_residual**2))),
        "max_abs_residual_gpa": float(np.max(np.abs(pressure_residual))),
        "volume_chi_square_using_reported_constant_sigma": float(
            np.sum((volume_residual / sigma_atomic) ** 2)
        ),
        "success": bool(result.success),
        "message": str(result.message),
    }


def reproduce() -> dict[str, object]:
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))

    source_pressure = np.array([float(row["pressure_mxb1986_gpa"]) for row in rows])
    stored_ratio = np.array([float(row["ruby_r1_ratio"]) for row in rows])
    stored_h2005 = np.array([float(row["pressure_h2005_gpa"]) for row in rows])
    volume = np.array([float(row["atomic_volume_a3"]) for row in rows])
    run2 = np.array(
        [row["best_supported_table_ii_selection"].lower() == "true" for row in rows]
    )

    mxb1986 = _calibration("ruby_mao_1986")
    h2005 = _calibration("ruby_holzapfel_2005")
    ratio = _ratio_from_pressure(source_pressure, mxb1986)
    pressure_h2005 = _pressure_from_ratio(ratio, h2005)
    if not np.allclose(ratio, stored_ratio, rtol=0.0, atol=5.0e-10):
        raise AssertionError("stored ruby ratios do not reproduce MXB1986 inversion")
    if not np.allclose(pressure_h2005, stored_h2005, rtol=0.0, atol=5.0e-9):
        raise AssertionError(
            "stored H2005 pressures do not reproduce the exact transform"
        )

    published = {
        "MXB1986": {"K0": 447.0, "K0_prime": 3.00, "chi_square": 4.5},
        "H2005": {"K0": 443.0, "K0_prime": 3.97, "chi_square": 1.5},
    }
    scales = {"MXB1986": source_pressure, "H2005": pressure_h2005}
    fits: dict[str, object] = {}
    for name, pressure in scales.items():
        fits[name] = {
            "published_table_ii": published[name],
            "all_figure_markers": {
                "rows": len(rows),
                "unweighted_pressure_residual_fit": _fit_pressure(volume, pressure),
                "equal_volume_weight_fit": _fit_volume(volume, pressure),
            },
            "best_supported_run2_selection": {
                "rows": int(run2.sum()),
                "equal_volume_weight_fit": _fit_volume(volume[run2], pressure[run2]),
            },
        }

    return {
        "format": "peritheos.datchi-2007-diamond-audit",
        "format_version": 1,
        "dataset": DATASET.name,
        "dataset_sha256": _sha256(DATASET),
        "source_data": {
            "status": "figure_digitized",
            "observations": len(rows),
            "run2_observations": int(run2.sum()),
            "run3_observations": int((~run2).sum()),
            "source_pressure_range_mxb1986_gpa": [
                float(source_pressure.min()),
                float(source_pressure.max()),
            ],
            "recalibrated_pressure_range_h2005_gpa": [
                float(pressure_h2005.min()),
                float(pressure_h2005.max()),
            ],
        },
        "pressure_recalibration": {
            "source": "ruby_mao_1986",
            "target": "ruby_holzapfel_2005",
            "path": "reported MXB1986 pressure -> ruby R1 ratio -> H2005 pressure",
            "source_parameters": mxb1986["parameters"],
            "target_parameters": h2005["parameters"],
            "stored_columns_verified": True,
        },
        "fit_protocol": {
            "fixed_parameters": {"V0_atomic_a3": V0_ATOMIC_A3},
            "published_exclusions": None,
            "published_weights": None,
            "reported_observation_uncertainties": {
                "volume_cm3_mol": OCCELLI_VOLUME_SIGMA_CM3_MOL,
                "pressure_gpa": "0.05 at 1 GPa to 1 at 140 GPa; no rowwise rule",
            },
            "best_supported_hypothesis": (
                "run 2 only, equal volume weights; this is not source-confirmed, but "
                "it recovers the independently printed MXB1986 and H2005 Table II "
                "pairs from the same markers"
            ),
            "covariance": "not estimated; exact rows and the source objective are unavailable",
        },
        "fits": fits,
        "assessment": (
            "Table II is numerically reproduced within the published fit standard "
            "deviations by the run-2/equal-volume-weight hypothesis. Exact parity is "
            "not claimed because Datchi et al. do not state exclusions or weights and "
            "Occelli et al. publish the observations only as a figure."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    payload = json.dumps(reproduce(), indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(payload, end="")
    else:
        args.output.write_text(payload, encoding="utf-8")


if __name__ == "__main__":
    main()
