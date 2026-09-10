"""Reproduce Sun et al. (2019) thermal-EOS diagnostics from Figures S3-S4.

The supporting figures contain rasterized regression lines rather than embedded
source tables.  The bundled digitization therefore constrains the plotted
``Cv(V)`` and ``gamma(V)`` trends, but it cannot recover the authors' unknown
relative weights or make the four-parameter BM4 cold curve well conditioned.
"""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import BM4

ROOT = Path(__file__).parents[1]
DATA = ROOT / "peritheos/data/datasets"
DIGITIZED = DATA / "fesio3-liquid-sun-2019-figures-s3-s4-digitized.csv"
PVT = DATA / "fesio3-liquid-sun-2019-table1-pvt.csv"

VX_CM3_MOL = 40.72
T0_K = 2500.0
CM3_MOL_PER_A3 = 0.602214076

CV_PUBLISHED = np.asarray([204.347, -830.0, 20.291, 29.110])
GAMMA_PUBLISHED = np.asarray([0.311, -0.933, 1.144])
BM4_PUBLISHED = np.asarray([47.25 / CM3_MOL_PER_A3, 2.217, 22.913, -175.628])


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _cv(params: np.ndarray, volume_cm3_mol: np.ndarray) -> np.ndarray:
    cv_prime, area, width, center = params
    return cv_prime + area / (width * math.sqrt(math.pi / 2.0)) * np.exp(
        -2.0 * (volume_cm3_mol - center) ** 2 / width**2
    )


def _gamma(params: np.ndarray, volume_ratio: np.ndarray) -> np.ndarray:
    return (
        params[0]
        + params[1] * (volume_ratio - 1.0)
        + params[2] * (volume_ratio - 1.0) ** 2
    )


def _thermal_pressure(
    cv_params: np.ndarray,
    gamma_params: np.ndarray,
    volume_cm3_mol: np.ndarray,
    temperature_k: np.ndarray,
) -> np.ndarray:
    return (
        _cv(cv_params, volume_cm3_mol)
        * _gamma(gamma_params, volume_cm3_mol / VX_CM3_MOL)
        * (temperature_k - T0_K)
        / volume_cm3_mol
        / 1000.0
    )


def _fit_digitized_thermal_terms() -> tuple[np.ndarray, np.ndarray, dict[str, float]]:
    rows = _rows(DIGITIZED)
    volume = np.asarray([float(row["volume_cm3_mol"]) for row in rows])
    ratio = np.asarray([float(row["volume_ratio_to_vx"]) for row in rows])
    cv = np.asarray([float(row["cv_j_mol_k"]) for row in rows])
    gamma = np.asarray([float(row["gamma"]) for row in rows])

    cv_fit = least_squares(
        lambda params: _cv(params, volume) - cv,
        CV_PUBLISHED,
        bounds=([100.0, -5000.0, 1.0, 0.0], [300.0, 0.0, 100.0, 100.0]),
        x_scale="jac",
        ftol=1.0e-13,
        xtol=1.0e-13,
        gtol=1.0e-13,
    )
    x = ratio - 1.0
    gamma_fit = np.linalg.lstsq(
        np.column_stack((np.ones_like(x), x, x**2)), gamma, rcond=None
    )[0]
    diagnostics = {
        "cv_rmse_j_mol_k": float(np.sqrt(np.mean(cv_fit.fun**2))),
        "gamma_rmse": float(np.sqrt(np.mean((_gamma(gamma_fit, ratio) - gamma) ** 2))),
        "published_cv_rmse_j_mol_k": float(
            np.sqrt(np.mean((_cv(CV_PUBLISHED, volume) - cv) ** 2))
        ),
        "published_gamma_rmse": float(
            np.sqrt(np.mean((_gamma(GAMMA_PUBLISHED, ratio) - gamma) ** 2))
        ),
    }
    return cv_fit.x, gamma_fit, diagnostics


def _fit_bm4(
    target_pressure: np.ndarray,
    volume_a3: np.ndarray,
    *,
    fixed_v0_and_k0_double_prime: bool,
) -> tuple[np.ndarray, np.ndarray, float, float, bool]:
    if fixed_v0_and_k0_double_prime:
        start = BM4_PUBLISHED[[1, 2]]

        def parameters(values: np.ndarray) -> np.ndarray:
            return np.asarray(
                [BM4_PUBLISHED[0], values[0], values[1], BM4_PUBLISHED[3]]
            )

        bounds = ([1.0e-12, -np.inf], [np.inf, np.inf])
    else:
        start = BM4_PUBLISHED

        def parameters(values: np.ndarray) -> np.ndarray:
            return values

        bounds = (
            [
                BM4_PUBLISHED[0] * 0.5,
                BM4_PUBLISHED[1] * 0.05,
                0.0,
                BM4_PUBLISHED[3] * 5.0,
            ],
            [
                BM4_PUBLISHED[0] * 1.5,
                BM4_PUBLISHED[1] * 5.0,
                BM4_PUBLISHED[2] * 5.0,
                0.0,
            ],
        )

    def objective(values: np.ndarray) -> np.ndarray:
        return (
            np.asarray(BM4(*parameters(values)).pressure(volume_a3)) - target_pressure
        )

    fit = least_squares(
        objective,
        start,
        bounds=bounds,
        x_scale="jac",
        max_nfev=10000,
        ftol=1.0e-13,
        xtol=1.0e-13,
        gtol=1.0e-13,
    )
    residuals = objective(fit.x)
    fitted = parameters(fit.x)
    bound_hit = bool(
        not fixed_v0_and_k0_double_prime
        and np.isclose(fitted[3], BM4_PUBLISHED[3] * 5.0, atol=1.0e-5, rtol=0.0)
    )
    return (
        fitted,
        residuals,
        float(np.sqrt(np.mean(residuals**2))),
        float(np.linalg.cond(fit.jac)),
        bound_hit,
    )


def reproduce() -> dict[str, Any]:
    """Return digitized thermal fits and BM4 identifiability diagnostics."""
    cv_fit, gamma_fit, thermal_diagnostics = _fit_digitized_thermal_terms()
    rows = [row for row in _rows(PVT) if row["used_in_published_eos_fit"] == "yes"]
    volume_cm3_mol = np.asarray([float(row["volume_cm3_mol"]) for row in rows])
    volume_a3 = np.asarray([float(row["volume_a3_per_formula_unit"]) for row in rows])
    temperature = np.asarray([float(row["temperature_k"]) for row in rows])
    pressure = np.asarray([float(row["pressure_gpa"]) for row in rows])
    target = pressure - _thermal_pressure(
        cv_fit, gamma_fit, volume_cm3_mol, temperature
    )

    free, _, free_rmse, free_condition, free_bound_hit = _fit_bm4(
        target, volume_a3, fixed_v0_and_k0_double_prime=False
    )
    stable, _, stable_rmse, stable_condition, _ = _fit_bm4(
        target, volume_a3, fixed_v0_and_k0_double_prime=True
    )

    cv_names = ("Cv_prime", "A", "w", "Vc")
    gamma_names = ("gamma_Vx", "gamma_prime", "gamma_double_prime")
    bm4_names = ("V0_a3_per_formula_unit", "K0", "K0_prime", "K0_double_prime")
    return {
        "digitized_dataset": DIGITIZED.name,
        "digitized_volume_states": len(_rows(DIGITIZED)),
        "thermal_parameter_comparison": {
            "Cv": {
                name: {"published": float(source), "digitized_refit": float(refit)}
                for name, source, refit in zip(cv_names, CV_PUBLISHED, cv_fit)
            },
            "gamma": {
                name: {"published": float(source), "digitized_refit": float(refit)}
                for name, source, refit in zip(gamma_names, GAMMA_PUBLISHED, gamma_fit)
            },
            "diagnostics": thermal_diagnostics,
        },
        "bm4_free_sensitivity": {
            "parameters": dict(zip(bm4_names, (float(value) for value in free))),
            "rmse_gpa": free_rmse,
            "jacobian_condition_number": free_condition,
            "k0_double_prime_lower_audit_bound_hit": free_bound_hit,
            "interpretation": "Not a source refit; the four-free-parameter solution remains nearly singular.",
        },
        "bm4_stabilized_check": {
            "fixed_parameters": ["V0_a3_per_formula_unit", "K0_double_prime"],
            "parameters": dict(zip(bm4_names, (float(value) for value in stable))),
            "rmse_gpa": stable_rmse,
            "jacobian_condition_number": stable_condition,
            "interpretation": "Equation-and-unit check only; not the authors' unconstrained regression.",
        },
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
