"""User-requested equal-relative-weight diagnostics on recovered Xiao input data.

This is an alternate objective, not a reconstruction of the published weighted
multiproperty regression. Published coefficients and material records are never
modified. The full nine-parameter fits deliberately expose identifiability and
bound dependence; no covariance-derived parameter error bars are reported.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import NaturalStrain4
from peritheos.eos.thermal import DebyeAnharmonicHelmholtz

ROOT = Path(__file__).resolve().parents[1]
NAMES = ("c1_mpa", "c2_mpa", "c3_mpa", "theta0_k", "gamma0", "q", "b1", "b2", "b3")
PUBLISHED = np.array([2656.5, 7298.0, 10.0, 86.44, 2.68, 0.0024, 0.0128, 0.388, 7.85])
# Diagnostic search bounds, not uncertainty intervals or published constraints.
LOWER = np.array([500.0, 100.0, 0.1, 20.0, 0.1, 1e-6, 1e-6, 1e-4, 0.1])
UPPER = np.array([6000.0, 30000.0, 100000.0, 300.0, 10.0, 5.0, 1.0, 10.0, 30.0])
OUTPUT = ROOT / "docs/data/argon-xiao-2025-equal-weight-refit.json"


def model(parameters):
    c1, c2, c3, theta, gamma, q, b1, b2, b3 = parameters
    k0 = c1 / 1000
    kp = 2 + 2 * c2 / c1
    kpp = (6 * c3 / c1 - kp * kp + 3 * kp - 3) / k0
    return DebyeAnharmonicHelmholtz(
        NaturalStrain4(2.2555, k0, kp, kpp), 70.0, theta, gamma, q, b1, b2, b3
    )


def recovered_data():
    def read(name):
        path = ROOT / "peritheos/data/datasets" / name
        with path.open() as stream:
            rows = list(csv.DictReader(stream))
        return rows, {
            "path": str(path.relative_to(ROOT)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }

    volume, volume_source = read("argon-fcc-xiao-2025-table5.csv")
    bulk, bulk_source = read("argon-anderson-swenson-1975-table1.csv")
    return {
        "volume_t": np.array([float(row["temperature_k"]) for row in volume]),
        "volume": np.array([float(row["molar_volume_cm3_mol"]) for row in volume]),
        "bulk_t": np.array([float(row["temperature_k"]) for row in bulk]),
        "bulk": np.array([float(row["k0_kbar"]) * 0.1 for row in bulk]),
        "volume_source": dict(
            volume_source,
            doi="10.1007/s10765-024-03469-2",
            source_location="Table5",
            kind="measured_neutron_volume",
            pressure_assumption="P=0 approximation; source actual pressures slightly above sublimation and unreported",
        ),
        "bulk_source": dict(
            bulk_source,
            doi="10.1016/0022-3697(75)90004-9",
            source_location="Table1 K0 column",
            kind="published_extrapolated_isothermal_bulk_modulus_estimates",
            note="Five published summary estimates from isotherm fits, not direct independent zero-pressure measurements. Xiao used two in its original regression; all five retained for this user-directed objective. Correlations and uncertainties unavailable.",
        ),
    }


def predictions(parameters, data, include_bulk):
    eos = model(parameters)
    volume = (
        np.asarray(eos.volume(np.zeros_like(data["volume_t"]), data["volume_t"])) * 10
    )
    if not include_bulk:
        return volume
    equilibrium = eos.volume(np.zeros_like(data["bulk_t"]), data["bulk_t"])
    bulk = eos.bulk_modulus(equilibrium, data["bulk_t"])
    return np.concatenate([volume, bulk])


def residuals(log_multipliers, data, include_bulk):
    observed = (
        np.concatenate([data["volume"], data["bulk"]])
        if include_bulk
        else data["volume"]
    )
    try:
        calculated = predictions(
            PUBLISHED * np.exp(log_multipliers), data, include_bulk
        )
    except (ValueError, ArithmeticError):
        return np.full(len(observed), 1000.0)
    if not np.all(np.isfinite(calculated)) or np.any(calculated <= 0):
        return np.full(len(observed), 1000.0)
    return (calculated - observed) / observed


def sensitivity(log_multipliers, data, include_bulk, step=1e-3):
    # Absolute steps in ln(parameter/published) avoid near-zero coordinate
    # finite-difference noise in scipy's relative-step Jacobian.
    basis = np.eye(9) * step
    return np.column_stack(
        [
            (
                residuals(log_multipliers + d, data, include_bulk)
                - residuals(log_multipliers - d, data, include_bulk)
            )
            / (2 * step)
            for d in basis
        ]
    )


def spectrum(jacobian):
    singular = np.linalg.svd(jacobian, compute_uv=False)
    return {
        "parameter_coordinates": "ln(parameter / published_parameter)",
        "singular_values": singular.tolist(),
        "condition_number": float(singular[0] / singular[-1]),
        "rank_at_relative_singular_value_threshold": {
            str(tolerance): int(np.sum(singular > tolerance * singular[0]))
            for tolerance in [1e-4, 1e-6, 1e-8]
        },
        "interpretation": "Sensitivity conditioning only, not statistical confidence or proof of identifiability. Near-null directions and bound hits preclude unique coefficient recovery from these data.",
    }


def metrics(residual, include_bulk):
    def summarize(values):
        return {
            "count": len(values),
            "sum_squared_relative_residuals": float(values @ values),
            "relative_rms_percent": float(np.sqrt(np.mean(values * values)) * 100),
            "relative_aad_percent": float(np.mean(np.abs(values)) * 100),
        }

    result = {"all": summarize(residual), "volume": summarize(residual[:22])}
    if include_bulk:
        result["bulk"] = summarize(residual[22:])
    return result


def thermodynamic_diagnostics(parameters, data, include_bulk):
    temperatures = np.unique(
        np.concatenate([data["volume_t"], data["bulk_t"]])
        if include_bulk
        else data["volume_t"]
    )
    eos = model(parameters)
    volumes = eos.volume(np.zeros_like(temperatures), temperatures)
    cv = np.asarray(eos.molar_heat_capacity_v(volumes, temperatures))
    bulk = np.asarray(eos.bulk_modulus(volumes, temperatures))
    alpha = np.asarray(eos.thermal_expansivity(volumes, temperatures))
    return {
        "scope": "Sampled P=0 states at input temperatures only; no global stability claim.",
        "minimum_cv_j_mol_k": float(cv.min()),
        "minimum_bulk_gpa": float(bulk.min()),
        "minimum_thermal_expansivity_per_k": float(alpha.min()),
        "passes_sampled_positive_cv_and_bulk": bool(
            np.all(cv > 0) and np.all(bulk > 0)
        ),
        "negative_cv_temperatures_k": temperatures[cv < 0].tolist(),
        "note": "No original physical penalties were applied to the alternate equal-weight objective. A negative heat capacity makes a candidate physically unacceptable even when its fitted residuals are small.",
    }


def fit_case(data, include_bulk, free_indices, starts=5, max_nfev=300):
    free_indices = np.array(free_indices)
    lower = np.log(LOWER / PUBLISHED)[free_indices]
    upper = np.log(UPPER / PUBLISHED)[free_indices]
    rng = np.random.default_rng(20250903)
    seeds = [np.zeros(len(free_indices))] + [
        rng.normal(0, 0.12, len(free_indices)) for _ in range(starts - 1)
    ]

    def expand(values):
        full = np.zeros(9)
        full[free_indices] = values
        return full

    solutions = []
    for index, start in enumerate(seeds):
        fit = least_squares(
            lambda values: residuals(expand(values), data, include_bulk),
            np.clip(start, lower + 1e-8, upper - 1e-8),
            bounds=(lower, upper),
            jac=lambda values: sensitivity(
                expand(values), data, include_bulk, step=1e-4
            )[:, free_indices],
            max_nfev=max_nfev,
            xtol=1e-10,
            ftol=1e-10,
            gtol=1e-10,
        )
        logp = expand(fit.x)
        parameters = PUBLISHED * np.exp(logp)
        at_bounds = [
            NAMES[j]
            for j, value in zip(free_indices, fit.x)
            if min(
                value - np.log(LOWER[j] / PUBLISHED[j]),
                np.log(UPPER[j] / PUBLISHED[j]) - value,
            )
            < 1e-4
        ]
        solutions.append(
            {
                "start_index": index,
                "starting_parameters": dict(
                    zip(NAMES, (PUBLISHED * np.exp(expand(start))).tolist())
                ),
                "parameters": dict(zip(NAMES, parameters.tolist())),
                "thermodynamic_diagnostics": thermodynamic_diagnostics(
                    parameters, data, include_bulk
                ),
                "optimizer_success": bool(fit.success),
                "optimizer_status": int(fit.status),
                "optimizer_message": fit.message,
                "function_evaluations": int(fit.nfev),
                "optimality": float(fit.optimality),
                "parameters_at_search_bounds": at_bounds,
                "metrics": metrics(fit.fun, include_bulk),
                "full_nine_parameter_sensitivity": spectrum(
                    sensitivity(logp, data, include_bulk)
                ),
            }
        )
    best = min(
        solutions,
        key=lambda solution: solution["metrics"]["all"][
            "sum_squared_relative_residuals"
        ],
    )
    parameter_array = np.array(
        [[solution["parameters"][name] for name in NAMES] for solution in solutions]
    )
    return {
        "free_parameters": [NAMES[j] for j in free_indices],
        "fixed_parameters": {
            name: float(PUBLISHED[j])
            for j, name in enumerate(NAMES)
            if j not in free_indices
        },
        "starts": starts,
        "seed": 20250903,
        "max_function_evaluations_per_start": max_nfev,
        "selection_note": "Best means lowest attained objective, including iteration-limited candidates; optimizer_success is reported separately. No unique optimum is asserted.",
        "best_optimizer_success": best["optimizer_success"],
        "best_thermodynamic_diagnostics": best["thermodynamic_diagnostics"],
        "converged_starts": sum(
            solution["optimizer_success"] for solution in solutions
        ),
        "best_converged_start_index": min(
            (solution for solution in solutions if solution["optimizer_success"]),
            key=lambda solution: solution["metrics"]["all"][
                "sum_squared_relative_residuals"
            ],
            default={"start_index": None},
        )["start_index"],
        "jacobian_log_parameter_step": 1e-4,
        "sensitivity_log_parameter_step": 1e-3,
        "best_start_index": best["start_index"],
        "best_parameters": best["parameters"],
        "best_metrics": best["metrics"],
        "best_parameters_at_search_bounds": best["parameters_at_search_bounds"],
        "multistart_parameter_ranges": {
            name: [
                float(parameter_array[:, j].min()),
                float(parameter_array[:, j].max()),
            ]
            for j, name in enumerate(NAMES)
        },
        "solutions": solutions,
    }


def reproduce(starts=5, max_nfev=300):
    data = recovered_data()
    result = {
        "status": "alternate_equal_weight_diagnostic",
        "source_objective_reproduced": False,
        "objective": "Minimize sum_i ((calculated_i-observed_i)/observed_i)^2; every recovered input entry has weight one. No uncertainty weights, per-property balancing, original empirical weights, or original penalties.",
        "data_scope": "22 Xiao Table5 neutron volume observations, and separately all five Anderson-Swenson Table1 extrapolated K0 estimates. This is not the complete legacy dataset cited by Xiao.",
        "data_sources": [data["volume_source"], data["bulk_source"]],
        "fixed_v00_cm3_mol": 22.555,
        "excluded_parameter": "Gas-phase triple-point entropy is not constrained by these solid-only properties and is not fitted.",
        "published_parameters": dict(zip(NAMES, PUBLISHED.tolist())),
        "search_bounds": {
            name: [float(LOWER[j]), float(UPPER[j])] for j, name in enumerate(NAMES)
        },
        "bounds_provenance": "Chosen numerical diagnostic limits with positive parameters, not published limits or confidence intervals. Failed/out-of-domain model states return a large residual. Bound hits and alternate minima are reported, never interpreted as parameter recovery.",
        "parameter_uncertainties": None,
        "parameter_uncertainty_note": "No statistical errors asserted: correlated source summary estimates, incomplete property coverage, practical rank deficiency and search-bound dependence.",
        "cases": {},
    }
    for include_bulk, key in [
        (False, "22_neutron_volumes"),
        (True, "22_volumes_plus_5_published_k0_estimates"),
    ]:
        result["cases"][key] = {
            "input_entries": 27 if include_bulk else 22,
            "published_metrics": metrics(
                residuals(np.zeros(9), data, include_bulk), include_bulk
            ),
            "published_full_nine_parameter_sensitivity": spectrum(
                sensitivity(np.zeros(9), data, include_bulk)
            ),
            "full_nine_parameter_fit": fit_case(
                data, include_bulk, range(9), starts, max_nfev
            ),
            "two_parameter_sensitivity_fit": fit_case(
                data, include_bulk, [3, 4], starts, max_nfev
            ),
        }
    return result


if __name__ == "__main__":
    OUTPUT.write_text(json.dumps(reproduce(), indent=2, allow_nan=False) + "\n")
    print(OUTPUT)
