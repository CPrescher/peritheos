#!/usr/bin/env python3
"""Independently reproduce Zhang et al. (2025) hcp-Fe Fits 1, 2, and 5."""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from build_zhang_2025_hcp_iron import DATASET, FIT_SPECS, SOURCE_COVARIANCE
from scipy.optimize import least_squares

SOURCE_GAS_CONSTANT = 8.314
RUNTIME_GAS_CONSTANT = 8.31446261815324
DEBYE_NODES, DEBYE_WEIGHTS = np.polynomial.legendre.leggauss(96)


def load_data(path: Path = DATASET) -> dict[str, np.ndarray]:
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    return {
        "pressure": np.asarray([float(row["pressure_gpa_fit"]) for row in rows]),
        "volume": np.asarray([float(row["molar_volume_cm3_mol"]) for row in rows]),
        "temperature": np.asarray([float(row["temperature_k"]) for row in rows]),
        **{
            f"fit{fit}_residual": np.asarray(
                [float(row[f"fit{fit}_observed_minus_model_gpa"]) for row in rows]
            )
            for fit in ("1", "2", "5")
        },
    }


def _debye_energy(
    volume: np.ndarray,
    temperature: np.ndarray,
    v0: float,
    theta0: float,
    gamma0: float,
    q: float,
    *,
    gas_constant: float,
) -> tuple[np.ndarray, np.ndarray]:
    gamma = gamma0 * (volume / v0) ** q
    theta = theta0 * np.exp((gamma0 - gamma) / q)

    def energy(at_temperature: np.ndarray) -> np.ndarray:
        limit = theta / at_temperature
        points = 0.5 * limit[:, None] * (DEBYE_NODES[None, :] + 1.0)
        integral = (
            0.5
            * limit
            * np.sum(
                DEBYE_WEIGHTS[None, :] * points**3 / np.expm1(points),
                axis=1,
            )
        )
        return (
            9.0
            * gas_constant
            * at_temperature
            * (at_temperature / theta) ** 3
            * integral
        )

    return gamma, energy(temperature) - energy(np.full_like(temperature, 300.0))


def pressure(
    model: str,
    parameters: np.ndarray,
    volume: np.ndarray,
    temperature: np.ndarray,
    *,
    gas_constant: float = SOURCE_GAS_CONSTANT,
) -> np.ndarray:
    v0, k0, k0_prime, theta0, gamma0, q = parameters
    if model == "BM3":
        eta = (v0 / volume) ** (1.0 / 3.0)
        reference = (
            1.5
            * k0
            * (eta**7 - eta**5)
            * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
        )
    elif model == "Vinet":
        x = (volume / v0) ** (1.0 / 3.0)
        reference = (
            3.0 * k0 * (1.0 - x) / x**2 * np.exp(1.5 * (k0_prime - 1.0) * (1.0 - x))
        )
    else:
        raise ValueError(model)
    gamma, delta_energy = _debye_energy(
        volume,
        temperature,
        v0,
        theta0,
        gamma0,
        q,
        gas_constant=gas_constant,
    )
    # E is J/mol and V is cm^3/mol: J/cm^3 = MPa.
    thermal = 1.0e-3 * gamma / volume * delta_energy
    return reference + thermal


def _published(fit: str) -> np.ndarray:
    spec = FIT_SPECS[fit]
    return np.asarray(
        [
            spec["V0"],
            spec["K0"],
            spec["K0p"],
            spec["theta0"],
            spec["gamma0"],
            spec["q"],
        ],
        dtype=float,
    )


def reproduce(path: Path = DATASET) -> dict[str, Any]:
    data = load_data(path)
    observations = data["pressure"].size
    outcomes: dict[str, Any] = {}
    for fit in ("1", "2", "5"):
        spec = FIT_SPECS[fit]
        published = _published(fit)
        free = np.arange(6) if fit == "1" else np.arange(1, 6)

        def expand(values: np.ndarray) -> np.ndarray:
            result = published.copy()
            result[free] = values
            return result

        def residual(values: np.ndarray) -> np.ndarray:
            parameters = expand(values)
            return data["pressure"] - pressure(
                spec["type"], parameters, data["volume"], data["temperature"]
            )

        lower = np.asarray([5.5, 50.0, 2.0, 100.0, 0.2, 0.05])[free]
        upper = np.asarray([8.0, 400.0, 8.0, 3000.0, 6.0, 3.0])[free]
        result = least_squares(
            residual,
            published[free],
            bounds=(lower, upper),
            x_scale="jac",
            xtol=1.0e-13,
            ftol=1.0e-13,
            gtol=1.0e-13,
            max_nfev=3000,
        )
        fitted = expand(result.x)
        fit_residual = residual(result.x)
        degrees_of_freedom = observations - free.size
        covariance_free = (
            np.linalg.inv(result.jac.T @ result.jac)
            * np.sum(fit_residual**2)
            / degrees_of_freedom
        )
        covariance = np.zeros((6, 6))
        covariance[np.ix_(free, free)] = covariance_free
        source_covariance = np.asarray(SOURCE_COVARIANCE[fit], dtype=float)
        covariance_relative_frobenius = float(
            np.linalg.norm(covariance - source_covariance)
            / np.linalg.norm(source_covariance)
        )

        published_model_residual = data["pressure"] - pressure(
            spec["type"], published, data["volume"], data["temperature"]
        )
        deposited_residual = data[f"fit{fit}_residual"]
        residual_delta = published_model_residual - deposited_residual
        runtime_constant_delta = pressure(
            spec["type"],
            published,
            data["volume"],
            data["temperature"],
            gas_constant=RUNTIME_GAS_CONSTANT,
        ) - pressure(spec["type"], published, data["volume"], data["temperature"])
        outcomes[fit] = {
            "model": spec["type"],
            "observations": int(observations),
            "free_parameters": [
                ("V0", "K0", "K0_prime", "theta0", "gamma0", "q")[index]
                for index in free
            ],
            "published_parameters": published.tolist(),
            "refit_parameters": fitted.tolist(),
            "difference": (fitted - published).tolist(),
            "refit_rmse_gpa_sqrt_sse_over_n_minus_1": float(
                np.sqrt(np.sum(fit_residual**2) / (observations - 1))
            ),
            "published_rmse_gpa": spec["rmse"],
            "published_rounded_parameter_rmse_gpa_sqrt_sse_over_n_minus_1": float(
                np.sqrt(np.sum(published_model_residual**2) / (observations - 1))
            ),
            "deposited_residual_reconstruction_rms_gpa": float(
                np.sqrt(np.mean(residual_delta**2))
            ),
            "deposited_residual_reconstruction_max_abs_gpa": float(
                np.max(np.abs(residual_delta))
            ),
            "covariance_relative_frobenius": covariance_relative_frobenius,
            "runtime_gas_constant_max_abs_pressure_delta_gpa": float(
                np.max(np.abs(runtime_constant_delta))
            ),
            "solver_success": bool(result.success),
            "solver_message": result.message,
        }
    return {
        "format": "peritheos.zhang-2025-hcp-iron-reproduction",
        "source_gas_constant_j_mol_k": SOURCE_GAS_CONSTANT,
        "dataset": str(path),
        "fits": outcomes,
    }


def check(outcome: dict[str, Any]) -> None:
    parameter_tolerances = np.asarray([0.001, 0.08, 0.001, 4.0, 0.005, 0.005])
    for result in outcome["fits"].values():
        difference = np.abs(np.asarray(result["difference"]))
        assert np.all(difference <= parameter_tolerances), difference
        assert result["deposited_residual_reconstruction_max_abs_gpa"] < 2.0e-6
        assert result["covariance_relative_frobenius"] < 0.006
        assert result["runtime_gas_constant_max_abs_pressure_delta_gpa"] < 0.1
        assert result["solver_success"]


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset", type=Path, default=DATASET)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    outcome = reproduce(args.dataset)
    if args.check:
        check(outcome)
    print(json.dumps(outcome, indent=2))


if __name__ == "__main__":
    main()
