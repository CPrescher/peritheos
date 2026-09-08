"""Refit the two Redfern et al. (1993) magnesite parameterizations."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets/magnesite-redfern-1993-table1-pv.csv"

PUBLISHED = {
    "fixed_kp": {"V0": 279.4, "K0": 142.0, "K0_prime": 4.0},
    "free_kp": {"V0": 279.4, "K0": 151.0, "K0_prime": 2.5},
}

PRESSURE_CALIBRATION = {
    "material": "NaCl",
    "reference": "Decker (1971)",
    "reference_doi": "10.1063/1.1660714",
    "pressure_uncertainty_fraction": 0.02,
    "recalculation_status": "missing_calibrant_observations",
}

FIT_PROTOCOL = {
    "fixed_kp": {
        "selection": "all 18 Table 1 rows",
        "objective": "unweighted pressure residuals",
        "fixed_parameters": ["V0", "K0_prime"],
    },
    "free_kp": {
        "selection": "9 Table 1 rows with Eulerian strain greater than 0.01",
        "objective": "volume residuals weighted by printed volume uncertainty",
        "fixed_parameters": ["V0"],
    },
}


def bm3_pressure(volume, v0, k0, kp):
    """Evaluate the standard Eulerian finite-strain BM3 pressure."""
    volume = np.asarray(volume, dtype=float)
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (kp - 4.0) * (eta**2 - 1.0))


def finite_strain_line(k0: float, kp: float) -> dict[str, float]:
    """Return the BM3 normalized-stress intercept and slope versus Eulerian strain."""
    return {
        "intercept_gpa": float(k0),
        "slope_gpa": float(1.5 * k0 * (kp - 4.0)),
    }


def _load_table() -> list[dict[str, str]]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _volume_at_pressure(pressure: float, v0: float, k0: float, kp: float) -> float:
    if pressure == 0.0:
        return v0
    return float(
        brentq(
            lambda volume: bm3_pressure(volume, v0, k0, kp) - pressure,
            0.8 * v0,
            v0,
        )
    )


def _standard_errors(result, *, absolute_sigma: bool) -> np.ndarray:
    degrees_of_freedom = result.fun.size - result.x.size
    covariance = np.linalg.pinv(result.jac.T @ result.jac)
    if not absolute_sigma:
        covariance *= np.sum(result.fun**2) / degrees_of_freedom
    return np.sqrt(np.diag(covariance))


def _fit_fixed_kp(rows: list[dict[str, str]]) -> dict[str, object]:
    selected = [row for row in rows if row["fit_included_bm2"] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in selected])
    volume = np.array([float(row["volume_a3"]) for row in selected])
    v0 = PUBLISHED["fixed_kp"]["V0"]
    kp = PUBLISHED["fixed_kp"]["K0_prime"]
    result = least_squares(
        lambda parameters: bm3_pressure(volume, v0, parameters[0], kp) - pressure,
        [PUBLISHED["fixed_kp"]["K0"]],
        bounds=([20.0], [500.0]),
    )
    k0 = float(result.x[0])
    residual = bm3_pressure(volume, v0, k0, kp) - pressure
    standard_error = float(_standard_errors(result, absolute_sigma=False)[0])
    return {
        "status": "similar",
        "observations": len(selected),
        "selection": FIT_PROTOCOL["fixed_kp"]["selection"],
        "objective": FIT_PROTOCOL["fixed_kp"]["objective"],
        "parameters": {"V0": v0, "K0": k0, "K0_prime": kp},
        "standard_errors": {"V0": None, "K0": standard_error, "K0_prime": None},
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "solver_success": bool(result.success),
    }


def _fit_free_kp(rows: list[dict[str, str]]) -> dict[str, object]:
    selected = [row for row in rows if row["fit_included_bm3"] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in selected])
    volume = np.array([float(row["volume_a3"]) for row in selected])
    sigma = np.array([float(row["volume_uncertainty_a3"]) for row in selected])
    v0 = PUBLISHED["free_kp"]["V0"]

    def residual(parameters):
        k0, kp = parameters
        calculated = np.array(
            [_volume_at_pressure(item, v0, k0, kp) for item in pressure]
        )
        return (calculated - volume) / sigma

    result = least_squares(
        residual,
        [PUBLISHED["free_kp"]["K0"], PUBLISHED["free_kp"]["K0_prime"]],
        bounds=([20.0, 0.0], [500.0, 20.0]),
    )
    k0, kp = map(float, result.x)
    pressure_residual = bm3_pressure(volume, v0, k0, kp) - pressure
    standard_errors = _standard_errors(result, absolute_sigma=True)
    return {
        "status": "similar",
        "observations": len(selected),
        "selection": FIT_PROTOCOL["free_kp"]["selection"],
        "objective": FIT_PROTOCOL["free_kp"]["objective"],
        "parameters": {"V0": v0, "K0": k0, "K0_prime": kp},
        "standard_errors": {
            "V0": None,
            "K0": float(standard_errors[0]),
            "K0_prime": float(standard_errors[1]),
        },
        "pressure_rmse_gpa": float(np.sqrt(np.mean(pressure_residual**2))),
        "weighted_volume_reduced_chi_square": float(
            np.sum(result.fun**2) / (len(selected) - 2)
        ),
        "solver_success": bool(result.success),
    }


def reproduce() -> dict[str, object]:
    """Return direct Table 1 refits and analytical coefficient checks."""
    rows = _load_table()
    fractions = np.array([1.0, 0.98, 0.95, 0.90])
    curves = {}
    for name, pars in PUBLISHED.items():
        pressure = bm3_pressure(
            fractions * pars["V0"], pars["V0"], pars["K0"], pars["K0_prime"]
        )
        curves[name] = dict(zip(map(str, fractions), map(float, pressure)))
    finite_strain = {
        name: finite_strain_line(pars["K0"], pars["K0_prime"])
        for name, pars in PUBLISHED.items()
    }
    return {
        "dataset": {
            "observations": len(rows),
            "pressure_range_gpa": [0.0, 19.7],
            "volume_range_a3": [250.8, 279.4],
        },
        "pressure_calibration": PRESSURE_CALIBRATION,
        "fit_protocol": FIT_PROTOCOL,
        "refits": {
            "fixed_kp": _fit_fixed_kp(rows),
            "free_kp": _fit_free_kp(rows),
        },
        "volume_fractions": fractions.tolist(),
        "curves_gpa": curves,
        "finite_strain_lines": finite_strain,
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
