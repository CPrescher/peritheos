#!/usr/bin/env python3
"""Reconstruct the Akins et al. (2004) candidate MgSiO3 melt Hugoniot.

The paper did not fit its BM3 reference isentrope directly to the tabulated
Hugoniot P-rho pairs.  It offset that isentrope with a Mie-Gruneisen energy
balance.  This reproduction implements that forward model and performs a
deliberately limited inverse check: rho0 and K0S are refined against the three
enstatite-starting melt states, while K0S-prime and every thermal term remain at
the values printed by Akins et al.  It is a validation reconstruction, not a
replacement EOS fit.
"""

from __future__ import annotations

import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATASET = (
    ROOT
    / "peritheos/data/datasets"
    / "mgsio3-liquid-akins-2004-mosenfelder-2009-rereduced-melt.csv"
)

AVOGADRO_MOL = 6.02214076e23
MOLAR_MASS_G_MOL = 100.387

# Final-paper candidate parameters (paragraph 4 and electronic-data appendix).
PUBLISHED_RHO0_G_CM3 = 3.68
PUBLISHED_K0S_GPA = 125.0
PUBLISHED_K0S_PRIME = 4.0
PUBLISHED_GAMMA0 = 2.4
PUBLISHED_Q = 1.0
PUBLISHED_TRANSITION_ENERGY_MJ_KG = 2.4
PUBLISHED_CV_3NR_FACTOR = 0.92

# Broad physical bounds make the two-parameter inverse problem explicit and
# prevent a three-point diagnostic from wandering to an unphysical branch.
RHO0_BOUNDS_G_CM3 = (3.0, 4.5)
K0S_BOUNDS_GPA = (50.0, 250.0)


def formula_volume_a3(rho_g_cm3: float) -> float:
    """Convert MgSiO3 mass density to volume per formula unit."""
    return MOLAR_MASS_G_MOL / (rho_g_cm3 * AVOGADRO_MOL) * 1.0e24


def bm3_isentrope(
    shock_density_g_cm3: np.ndarray,
    rho0_g_cm3: float,
    k0s_gpa: float,
    k0s_prime: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Return reference-isentrope pressure (GPa) and energy (MJ/kg)."""
    density = np.asarray(shock_density_g_cm3, dtype=float)
    v0_m3_kg = 1.0 / (rho0_g_cm3 * 1000.0)
    v_m3_kg = 1.0 / (density * 1000.0)
    f = 0.5 * ((v0_m3_kg / v_m3_kg) ** (2.0 / 3.0) - 1.0)
    a1 = 1.5 * (k0s_prime - 4.0)
    pressure_gpa = 3.0 * k0s_gpa * f * (2.0 * f + 1.0) ** 2.5 * (1.0 + a1 * f)
    energy_j_kg = 9.0 * v0_m3_kg * (k0s_gpa * 1.0e9) * (0.5 * f**2 + (a1 / 3.0) * f**3)
    return pressure_gpa, energy_j_kg / 1.0e6


def candidate_hugoniot_pressure_gpa(
    initial_density_g_cm3: np.ndarray,
    shock_density_g_cm3: np.ndarray,
    rho0_g_cm3: float,
    k0s_gpa: float,
    k0s_prime: float = PUBLISHED_K0S_PRIME,
    gamma0: float = PUBLISHED_GAMMA0,
    q: float = PUBLISHED_Q,
    transition_energy_mj_kg: float = PUBLISHED_TRANSITION_ENERGY_MJ_KG,
) -> np.ndarray:
    """Evaluate the source's energy-balanced theoretical Hugoniot.

    This is the explicit solution of Akins's equations 2.6-2.18 after setting
    the negligible initial pressure to zero.  Volumes are specific volumes.
    """
    initial_density = np.asarray(initial_density_g_cm3, dtype=float)
    shock_density = np.asarray(shock_density_g_cm3, dtype=float)
    v_initial = 1.0 / (initial_density * 1000.0)
    v_shock = 1.0 / (shock_density * 1000.0)
    v0 = 1.0 / (rho0_g_cm3 * 1000.0)
    p_isentrope_gpa, e_isentrope_mj_kg = bm3_isentrope(
        shock_density, rho0_g_cm3, k0s_gpa, k0s_prime
    )
    gamma = gamma0 * (v_shock / v0) ** q
    numerator_gpa = (
        p_isentrope_gpa
        - gamma
        * (transition_energy_mj_kg + e_isentrope_mj_kg)
        * 1.0e6
        / v_shock
        / 1.0e9
    )
    denominator = 1.0 - gamma * (v_initial - v_shock) / (2.0 * v_shock)
    return numerator_gpa / denominator


def _load_fit_rows() -> dict[str, np.ndarray]:
    with DATASET.open(newline="", encoding="utf-8") as stream:
        rows = [row for row in csv.DictReader(stream) if row["fit_included"] == "1"]
    return {
        "shot": np.asarray([int(row["shot"]) for row in rows]),
        "initial_density_g_cm3": np.asarray(
            [float(row["initial_density_g_cm3"]) for row in rows]
        ),
        "shock_density_g_cm3": np.asarray(
            [float(row["shock_density_g_cm3"]) for row in rows]
        ),
        "pressure_gpa": np.asarray([float(row["pressure_gpa"]) for row in rows]),
        "pressure_uncertainty_gpa": np.asarray(
            [float(row["pressure_uncertainty_gpa"]) for row in rows]
        ),
    }


def _fit(rows: dict[str, np.ndarray], weighted: bool) -> Any:
    def residual(parameters: np.ndarray) -> np.ndarray:
        calculated = candidate_hugoniot_pressure_gpa(
            rows["initial_density_g_cm3"],
            rows["shock_density_g_cm3"],
            float(parameters[0]),
            float(parameters[1]),
        )
        values = calculated - rows["pressure_gpa"]
        if weighted:
            values = values / rows["pressure_uncertainty_gpa"]
        return values

    return least_squares(
        residual,
        x0=np.asarray([PUBLISHED_RHO0_G_CM3, PUBLISHED_K0S_GPA]),
        bounds=(
            np.asarray([RHO0_BOUNDS_G_CM3[0], K0S_BOUNDS_GPA[0]]),
            np.asarray([RHO0_BOUNDS_G_CM3[1], K0S_BOUNDS_GPA[1]]),
        ),
        x_scale="jac",
    )


def reproduce() -> dict[str, Any]:
    """Return published-curve residuals and bounded reconstruction diagnostics."""
    rows = _load_fit_rows()
    published_pressure = candidate_hugoniot_pressure_gpa(
        rows["initial_density_g_cm3"],
        rows["shock_density_g_cm3"],
        PUBLISHED_RHO0_G_CM3,
        PUBLISHED_K0S_GPA,
    )
    published_residuals = published_pressure - rows["pressure_gpa"]

    unweighted = _fit(rows, weighted=False)
    rho0, k0s = (float(value) for value in unweighted.x)
    reconstructed_pressure = candidate_hugoniot_pressure_gpa(
        rows["initial_density_g_cm3"], rows["shock_density_g_cm3"], rho0, k0s
    )
    residuals = reconstructed_pressure - rows["pressure_gpa"]
    dof = int(rows["pressure_gpa"].size - unweighted.x.size)
    covariance = np.linalg.inv(unweighted.jac.T @ unweighted.jac)
    covariance *= float(np.sum(unweighted.fun**2)) / dof
    standard_errors = np.sqrt(np.diag(covariance))
    correlation = covariance / np.outer(standard_errors, standard_errors)

    weighted = _fit(rows, weighted=True)
    weighted_pressure = candidate_hugoniot_pressure_gpa(
        rows["initial_density_g_cm3"],
        rows["shock_density_g_cm3"],
        float(weighted.x[0]),
        float(weighted.x[1]),
    )

    return {
        "dataset": DATASET.relative_to(ROOT).as_posix(),
        "observations": int(rows["pressure_gpa"].size),
        "shots": [int(value) for value in rows["shot"]],
        "data_semantics": {
            "directly_measured": [
                "initial density",
                "flyer velocity",
                "shock velocity",
            ],
            "impedance_matched": [
                "particle velocity",
                "pressure",
                "shock density",
                "internal energy",
            ],
            "calculated_model_output": "candidate Hugoniot pressure",
        },
        "published_candidate": {
            "rho0_g_cm3": PUBLISHED_RHO0_G_CM3,
            "V0_a3_per_formula_unit": formula_volume_a3(PUBLISHED_RHO0_G_CM3),
            "K0S_gpa": PUBLISHED_K0S_GPA,
            "K0S_prime": PUBLISHED_K0S_PRIME,
            "gamma0": PUBLISHED_GAMMA0,
            "q": PUBLISHED_Q,
            "transition_energy_mj_kg": PUBLISHED_TRANSITION_ENERGY_MJ_KG,
            "Cv_factor_times_3nR": PUBLISHED_CV_3NR_FACTOR,
            "calculated_pressure_gpa": published_pressure.tolist(),
            "pressure_residuals_gpa": published_residuals.tolist(),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(published_residuals**2))),
        },
        "bounded_unweighted_reconstruction": {
            "refined_parameters": ["rho0", "K0S"],
            "fixed_parameters": [
                "K0S_prime",
                "gamma0",
                "q",
                "transition_energy",
            ],
            "bounds": {
                "rho0_g_cm3": list(RHO0_BOUNDS_G_CM3),
                "K0S_gpa": list(K0S_BOUNDS_GPA),
            },
            "parameters": {
                "rho0_g_cm3": rho0,
                "V0_a3_per_formula_unit": formula_volume_a3(rho0),
                "K0S_gpa": k0s,
                "K0S_prime_fixed": PUBLISHED_K0S_PRIME,
            },
            "standard_errors": {
                "rho0_g_cm3": float(standard_errors[0]),
                "V0_a3_per_formula_unit": (
                    formula_volume_a3(rho0) * float(standard_errors[0]) / rho0
                ),
                "K0S_gpa": float(standard_errors[1]),
            },
            "parameter_correlation": float(correlation[0, 1]),
            "jacobian_condition_number": float(np.linalg.cond(unweighted.jac)),
            "calculated_pressure_gpa": reconstructed_pressure.tolist(),
            "pressure_residuals_gpa": residuals.tolist(),
            "pressure_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
            "degrees_of_freedom": dof,
            "solver_success": bool(unweighted.success),
            "solver_message": str(unweighted.message),
        },
        "weighted_sensitivity_check": {
            "parameters": {
                "rho0_g_cm3": float(weighted.x[0]),
                "K0S_gpa": float(weighted.x[1]),
            },
            "pressure_residuals_gpa": (
                weighted_pressure - rows["pressure_gpa"]
            ).tolist(),
            "weighted_residual_rmse": float(np.sqrt(np.mean(weighted.fun**2))),
            "note": (
                "This is a sensitivity check, not an asserted source objective; "
                "the source does not publish an optimization protocol or covariance."
            ),
        },
        "identifiability": {
            "observations": 3,
            "independently_refined_parameters": 2,
            "degrees_of_freedom": 1,
            "K0S_prime_refittable": False,
            "reason": (
                "Three selected melt states cannot independently determine rho0, "
                "K0S, and K0S-prime, much less the coupled thermal terms. Cv is "
                "used for the source's temperature calculation and does not enter "
                "this pressure-only reconstruction."
            ),
        },
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True, allow_nan=False))
