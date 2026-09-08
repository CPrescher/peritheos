#!/usr/bin/env python3
"""Reproduce the Ye (2017) 300 K Vinet fits adopted by Zhu et al. (2025)."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
DATASET_DIR = ROOT / "peritheos" / "data" / "datasets"

DATASETS = {
    "ye_2017_data_s1_au_mgo_300k": {
        "path": DATASET_DIR / "ye-2017-data-s1-au-mgo-300k.csv",
        "sha256": "e4d30f69fb2a317e3d9aebcbac94e774324ca27919fe63d4b396c0c5305d114c",
        "rows": 28,
    },
    "ye_2017_data_s2_pt_mgo_300k": {
        "path": DATASET_DIR / "ye-2017-data-s2-pt-mgo-300k.csv",
        "sha256": "c8749ee759d6d2d4b01c8e8e964056de6be2703224370a0baf7f6eb5d46c06e4",
        "rows": 39,
    },
}

# Ye fixed K0 and V0 to Dorogokupets-Dewaele (2007). The journal table rounds
# V0 to 67.85, 60.38, and 74.71 A^3. These unrounded reference volumes are the
# defaults in the Ye/Shim group's Apache-2.0 pytheos implementation and are
# required to reproduce the printed K0' values from the corrected supplements.
ANCHORS = {
    "Au": {"V0": 4.07860**3, "K0": 167.0, "K0_prime": 5.90},
    "Pt": {"V0": 3.9231**3, "K0": 277.3, "K0_prime": 5.12},
    "MgO": {"V0": 74.698, "K0": 160.3, "K0_prime": 4.18},
}

PUBLISHED = {
    "Au_from_MgO": {"K0_prime": 5.897, "error": 0.022},
    "Au_from_Pt": {"K0_prime": 5.813, "error": 0.022},
    "MgO_from_Au": {"K0_prime": 4.182, "error": 0.019},
    "MgO_from_Pt": {"K0_prime": 4.109, "error": 0.022},
    "Pt_from_Au": {"K0_prime": 5.230, "error": 0.033},
    "Pt_from_MgO": {"K0_prime": 5.226, "error": 0.033},
}

RECORDS = {
    "gold_ye_2017_vinet_300k": ("gold", "Au_from_MgO"),
    "mgo_ye_2017_vinet_300k": ("mgo", "MgO_from_Au"),
    "platinum_ye_2017_vinet_300k": ("platinum", "Pt_from_Au"),
}

# Coefficients are the rounded outputs carried into the released Mendeley Data
# v3 optimizer/property workflow. beta0 is converted from J g^-1 K^-2 to
# J mol^-1 K^-2 using the molar masses used by those scripts.
THERMAL_RECORDS = {
    "gold_zhu_2025_pvt": {
        "material": "gold",
        "V0": 67.85,
        "K0": 167.0,
        "K0_prime": 5.9,
        "Tr": 300.0,
        "theta0": 180.0,
        "gamma0": 2.93,
        "a": 0.75,
        "b": 2.7,
        "n": 1.0,
        "beta0": 0.0,
        "m": 1.0,
        "source_molar_mass": 196.97,
        "optimizer_beta_mass": 0.0,
        "optimizer_beta_molar_mass": 196.97,
        "maximum_iterations": 20,
        "shock_molar_mass": 196.97,
        "stress_y": 0.12,
        "stress_a": 96.0,
    },
    "platinum_zhu_2025_pvt": {
        "material": "platinum",
        "V0": 60.38,
        "K0": 277.3,
        "K0_prime": 5.23,
        "Tr": 300.0,
        "theta0": 240.0,
        "gamma0": 2.75,
        "a": 0.39,
        "b": 5.1,
        "n": 1.0,
        "beta0": 110.0e-7 * 195.0,
        "m": 0.65,
        "source_molar_mass": 195.0,
        "optimizer_beta_mass": 110.0,
        "optimizer_beta_molar_mass": 195.05,
        "maximum_iterations": 14,
        "shock_molar_mass": 195.08,
        "stress_y": 0.29,
        "stress_a": 72.0,
    },
    "mgo_zhu_2025_pvt": {
        "material": "mgo",
        "V0": 74.71,
        "K0": 160.3,
        "K0_prime": 4.182,
        "Tr": 300.0,
        "theta0": 761.0,
        "gamma0": 1.53,
        "a": 1.0,
        "b": 1.43,
        "n": 2.0,
        "beta0": -200.0e-7 * 40.305,
        "m": 4.8,
        "source_molar_mass": 40.305,
        "optimizer_beta_mass": -200.0,
        "optimizer_beta_molar_mass": 40.305,
        "maximum_iterations": 20,
    },
}

# The standalone v3 calculators diverge from the optimizer/property files in
# three constants. They remain a diagnostic source, not the refit target.
CALCULATOR_OVERRIDES = {
    "gold_zhu_2025_pvt": {"theta0": 170.0},
    "platinum_zhu_2025_pvt": {"theta0": 230.0},
    "mgo_zhu_2025_pvt": {"gamma0": 1.52},
}

THERMAL_DATASETS = {
    "zhu_2025_au_shock": ("zhu-2025-au-shock.csv", "ef2fb3cbfc1636d6a513c40783e6c860d5d3ced3a9e48e0e939ef36f8f8ece01", 12),
    "zhu_2025_au_zero_pressure_thermal_expansion": ("zhu-2025-au-zero-pressure-thermal-expansion.csv", "3baac523dbf58873d01e05499f73934bbebdfc5bbd03de876c9e4f06dfeee611", 10),
    "zhu_2025_pt_shock": ("zhu-2025-pt-shock.csv", "025914a80db03e52fe631c8f18646b722a558026d616079edea3d47c7509b924", 68),
    "zhu_2025_pt_zero_pressure_thermal_expansion": ("zhu-2025-pt-zero-pressure-thermal-expansion.csv", "c338f3db484be4927880b90f7fcf88ce6b29c3c880e726755d3b29cf778cda70", 17),
    "zhu_2025_mgo_pvt": ("zhu-2025-mgo-pvt.csv", "839acb4435c741aaa7266c956d9eba196756e5782ca63d47fb0b1a930fe4285c", 213),
}


def _vinet_pressure(volume: np.ndarray, parameters: dict[str, float]) -> np.ndarray:
    x = np.cbrt(volume / parameters["V0"])
    return (
        3.0
        * parameters["K0"]
        * (1.0 - x)
        / x**2
        * np.exp(1.5 * (parameters["K0_prime"] - 1.0) * (1.0 - x))
    )


def _debye_3(value: float) -> float:
    """Evaluate D3 independently of the Peritheos implementation."""
    integral = quad(lambda x: x**3 / np.expm1(x), 0.0, value)[0]
    return 3.0 * integral / value**3


def _released_v3_pressure(
    volume: float, temperature: float, parameters: dict[str, Any]
) -> float:
    """Translate the released v3 pressure-calculator equations literally."""
    ratio = volume / parameters["V0"]
    static = float(_vinet_pressure(np.asarray(volume), parameters))
    gamma = parameters["gamma0"] * (
        1.0 + parameters["a"] * (ratio ** parameters["b"] - 1.0)
    )
    theta = parameters["theta0"] * np.exp(
        -parameters["gamma0"]
        * (
            (1.0 - parameters["a"]) * np.log(ratio)
            + parameters["a"] * (ratio ** parameters["b"] - 1.0) / parameters["b"]
        )
    )

    # The released calculators use these rounded legacy constants.  A cell
    # contains four formula units for all three source structures.
    gas_constant = 8.314
    # J bar^-1 mol^-1 = cm^3 mol^-1 / 10.  The source calculator uses
    # 0.6022 cm^3 mol^-1 per A^3 atom^-1.
    molar_volume = volume * 0.6022 / 4.0 / 10.0

    def energy(state_temperature: float) -> float:
        return (
            3.0
            * parameters["n"]
            * gas_constant
            * state_temperature
            * _debye_3(theta / state_temperature)
        )

    debye_pressure = gamma * (energy(temperature) - energy(300.0)) / molar_volume
    beta_mass = parameters["beta0"] / parameters["source_molar_mass"] / 1.0e-7
    excess_factor = (
        1.0e-9
        * 0.5
        * ratio ** parameters["m"]
        * beta_mass
        * 1.0e-7
        * (parameters["source_molar_mass"] * 4.0 / volume / 6.023 * 1.0e7)
    )
    excess_pressure = (
        excess_factor
        * (temperature**2 - 300.0**2)
        * parameters["m"]
    )
    return static + debye_pressure / 1.0e4 + excess_pressure


def _calculator_diagnostics() -> dict[str, Any]:
    records: dict[str, Any] = {}
    for identifier, parameters in THERMAL_RECORDS.items():
        calculator_parameters = {
            **parameters,
            **CALCULATOR_OVERRIDES[identifier],
        }
        document = get_material_document(parameters["material"])
        source = next(
            row for row in document["eos_records"] if row["identifier"] == identifier
        )
        model = Material.from_eosmat(
            document, record_identifiers=[identifier]
        ).get_eos_record(identifier)
        checkpoints = []
        for fraction in (0.95, 0.80, 0.65):
            for temperature in (300.0, 1000.0, 2000.0, 3000.0):
                volume = fraction * parameters["V0"]
                expected = _released_v3_pressure(
                    volume, temperature, calculator_parameters
                )
                actual = float(model.pressure(volume, temperature))
                checkpoints.append(
                    {
                        "volume_fraction": fraction,
                        "temperature_k": temperature,
                        "released_calculator_pressure_gpa": expected,
                        "peritheos_pressure_gpa": actual,
                        "difference_gpa": actual - expected,
                    }
                )
        maximum = max(abs(row["difference_gpa"]) for row in checkpoints)
        assert source["equation_kind"] == "thermal"
        records[identifier] = {
            "checkpoints": checkpoints,
            "max_abs_difference_gpa": maximum,
            "optimizer_parameters": {
                key: value
                for key, value in parameters.items()
                if key not in {"material", "source_molar_mass"}
            },
            "calculator_parameters": {
                key: value
                for key, value in calculator_parameters.items()
                if key not in {"material", "source_molar_mass"}
            },
            "fit_datasets": source["fit_datasets"],
        }
    return records


def _load_thermal_csv(identifier: str) -> dict[str, np.ndarray]:
    filename, digest, expected_rows = THERMAL_DATASETS[identifier]
    path = DATASET_DIR / filename
    assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == expected_rows
    return {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in rows[0]
    }


def _debye_energy(
    theta: np.ndarray | float, temperature: np.ndarray | float, n: float
) -> np.ndarray:
    theta_values, temperatures = np.broadcast_arrays(
        np.asarray(theta, dtype=float), np.asarray(temperature, dtype=float)
    )
    output = np.empty_like(theta_values)
    for index in np.ndindex(theta_values.shape):
        ratio = float(theta_values[index] / temperatures[index])
        integral = quad(
            lambda value: 0.0
            if value > 700.0
            else value**3 / np.expm1(value),
            0.0,
            ratio,
        )[0]
        output[index] = (
            9.0
            * n
            * 8.314
            * temperatures[index]
            * ratio**-3
            * integral
        )
    return output


def _robust_bisquare_gamma_fit(
    ratios: np.ndarray,
    targets: np.ndarray,
    *,
    fixed_a: float,
    start: tuple[float, float],
    gamma_bounds: tuple[float, float],
) -> tuple[np.ndarray, np.ndarray]:
    """Reproduce Curve Fitting Toolbox's robust bisquare NLS pathway."""

    def residuals(values: np.ndarray) -> np.ndarray:
        gamma0, exponent = values
        return gamma0 * (1.0 + fixed_a * (ratios**exponent - 1.0)) - targets

    bounds = ([gamma_bounds[0], -np.inf], [gamma_bounds[1], np.inf])
    values = least_squares(residuals, start, bounds=bounds).x
    weights = np.ones_like(targets)
    for _ in range(100):
        raw = residuals(values)
        gamma0, exponent = values
        powers = ratios**exponent
        jacobian = np.column_stack(
            (
                1.0 + fixed_a * (powers - 1.0),
                gamma0 * fixed_a * powers * np.log(ratios),
            )
        )
        leverage = np.einsum(
            "ij,jk,ik->i",
            jacobian,
            np.linalg.pinv(jacobian.T @ jacobian),
            jacobian,
        )
        adjusted = raw / np.sqrt(np.maximum(1.0 - leverage, 1.0e-12))
        scale = np.median(np.abs(adjusted - np.median(adjusted))) / 0.6745
        standardized = adjusted / (4.685 * scale)
        weights = np.where(
            np.abs(standardized) < 1.0,
            (1.0 - standardized**2) ** 2,
            0.0,
        )
        updated = least_squares(
            lambda candidate: np.sqrt(weights) * residuals(candidate),
            values,
            bounds=bounds,
        ).x
        if np.max(np.abs(updated - values)) < 1.0e-10:
            values = updated
            break
        values = updated
    return values, weights


def _source_excess_pressure(
    ratios: np.ndarray,
    temperature: np.ndarray | float,
    parameters: dict[str, Any],
) -> np.ndarray:
    return (
        1.0e-9
        * (
            0.5
            * ratios ** (parameters["m"] - 1.0)
            * parameters["optimizer_beta_mass"]
            * (parameters["optimizer_beta_molar_mass"] / 1.0e7)
            * parameters["m"]
        )
        / (parameters["V0"] * 6.022e-7 / 4.0)
        * np.asarray(temperature, dtype=float) ** 2
    )


def _gamma_theta(
    ratios: np.ndarray, gamma0: float, a: float, exponent: float, theta0: float
) -> tuple[np.ndarray, np.ndarray]:
    gamma = gamma0 * (1.0 + a * (ratios**exponent - 1.0))
    theta = (
        theta0
        * np.exp((gamma0 - gamma) / exponent)
        * ratios ** (gamma0 * (a - 1.0))
    )
    return gamma, theta


def _refit_mgo(parameters: dict[str, Any]) -> dict[str, Any]:
    data = _load_thermal_csv("zhu_2025_mgo_pvt")
    pressure = data["pressure_gpa"]
    temperature = data["temperature_k"]
    volume = data["volume_a3_conventional_cell"]
    ratios = volume / parameters["V0"]
    static = _vinet_pressure(volume, parameters)
    history = [[1.3, 0.1, 1.0]]
    stable = 0
    for _ in range(2, parameters["maximum_iterations"] + 1):
        gamma0, current_a, exponent = history[-1]
        _, theta = _gamma_theta(
            ratios, gamma0, current_a, exponent, parameters["theta0"]
        )
        energy_difference = _debye_energy(
            theta, temperature, parameters["n"]
        ) - _debye_energy(theta, 300.0, parameters["n"])
        excess = _source_excess_pressure(
            ratios, temperature, parameters
        ) - _source_excess_pressure(ratios, 300.0, parameters)
        targets = (
            (pressure - static - excess)
            * 1.0e9
            * (volume * 6.022e-7 / 4.0)
            / energy_difference
        )
        (next_gamma0, next_exponent), weights = _robust_bisquare_gamma_fit(
            ratios,
            targets,
            fixed_a=parameters["a"],
            start=(1.5, 1.0),
            gamma_bounds=(1.0, 1.7),
        )
        current = np.asarray(history[-1])
        following = np.asarray([next_gamma0, parameters["a"], next_exponent])
        history.append(following.tolist())
        stable = stable + 1 if np.all(np.abs(following - current) < 0.2) else 0
        if stable > 5:
            break

    fitted_gamma0, fitted_a, fitted_exponent = np.mean(history[-4:], axis=0)
    gamma, theta = _gamma_theta(
        ratios,
        fitted_gamma0,
        fitted_a,
        fitted_exponent,
        parameters["theta0"],
    )
    energy_difference = _debye_energy(
        theta, temperature, parameters["n"]
    ) - _debye_energy(theta, 300.0, parameters["n"])
    excess = _source_excess_pressure(
        ratios, temperature, parameters
    ) - _source_excess_pressure(ratios, 300.0, parameters)
    predicted = static + excess + (
        1.0e-9 * gamma * energy_difference / (volume * 6.022e-7 / 4.0)
    )
    residuals = pressure - predicted
    return {
        "gamma0": float(fitted_gamma0),
        "a": float(fitted_a),
        "b": float(fitted_exponent),
        "iterations": len(history) - 1,
        "observations": int(volume.size),
        "nonzero_robust_weights": int(np.count_nonzero(weights)),
        "rmse_pressure_gpa": float(np.sqrt(np.mean(residuals**2))),
        "max_abs_pressure_residual_gpa": float(np.max(np.abs(residuals))),
    }


def _refit_shock_and_expansion(
    identifier: str, parameters: dict[str, Any]
) -> dict[str, Any]:
    prefix = "au" if parameters["material"] == "gold" else "pt"
    shock = _load_thermal_csv(f"zhu_2025_{prefix}_shock")
    expansion = _load_thermal_csv(
        f"zhu_2025_{prefix}_zero_pressure_thermal_expansion"
    )
    volume = shock["volume_a3_conventional_cell"]
    pressure = shock["shock_pressure_gpa"]
    particle_velocity = shock["particle_velocity_m_s"]
    high_temperature = expansion["temperature_k"]
    high_volume = expansion["volume_a3_conventional_cell"]
    shock_ratios = volume / parameters["V0"]
    high_ratios = high_volume / parameters["V0"]
    ratios = np.concatenate((shock_ratios, high_ratios))
    static = _vinet_pressure(volume, parameters)
    high_static = _vinet_pressure(high_volume, parameters)
    density = volume * 6.022e-7 / 4.0
    shock_energy = (
        0.5
        * particle_velocity**2
        * parameters["shock_molar_mass"]
        / 1000.0
    )
    compression_energy = np.asarray(
        [
            150.55
            * quad(
                lambda state_volume: float(
                    _vinet_pressure(np.asarray(state_volume), parameters)
                ),
                parameters["V0"],
                state_volume,
            )[0]
            for state_volume in volume
        ]
    )
    stress = (
        2.0
        / 3.0
        * parameters["stress_y"]
        * (
            1.0
            + parameters["stress_a"]
            / 1000.0
            * pressure
            / shock_ratios ** (-1.0 / 3.0)
        )
    )
    start_gamma = 2.8 if identifier.startswith("gold") else 2.5
    history = [[start_gamma, 0.5, 3.0]]
    shock_temperature = np.full_like(volume, 300.0)
    stable = 0
    for _ in range(2, parameters["maximum_iterations"] + 1):
        gamma0, current_a, exponent = history[-1]
        gamma, theta = _gamma_theta(
            shock_ratios,
            gamma0,
            current_a,
            exponent,
            parameters["theta0"],
        )
        energy_300 = _debye_energy(theta, 300.0, parameters["n"])
        excess = _source_excess_pressure(
            shock_ratios, shock_temperature, parameters
        ) - _source_excess_pressure(shock_ratios, 300.0, parameters)
        estimated_energy = energy_300 + (
            (pressure - static - stress - excess)
            * 1.0e9
            * (volume / 4.0 * 0.6022 / 1.0e6)
            / gamma
        )
        shock_temperature = np.asarray(
            [
                brentq(
                    lambda trial: float(
                        _debye_energy(theta_value, trial, parameters["n"])
                        - target_energy
                    ),
                    1.0,
                    1.0e5,
                )
                for theta_value, target_energy in zip(theta, estimated_energy)
            ]
        )
        excess = _source_excess_pressure(
            shock_ratios, shock_temperature, parameters
        ) - _source_excess_pressure(shock_ratios, 300.0, parameters)
        excess_energy = (
            0.5
            * shock_ratios ** parameters["m"]
            * parameters["optimizer_beta_mass"]
            * (parameters["optimizer_beta_molar_mass"] / 1.0e7)
            * (shock_temperature**2 - 300.0**2)
        )
        shock_targets = (
            (pressure - static - stress - excess)
            * 1.0e9
            / (shock_energy + compression_energy - excess_energy)
            * density
        )

        _, high_theta = _gamma_theta(
            high_ratios,
            gamma0,
            current_a,
            exponent,
            parameters["theta0"],
        )
        high_energy_difference = _debye_energy(
            high_theta, high_temperature, parameters["n"]
        ) - _debye_energy(high_theta, 300.0, parameters["n"])
        high_excess = _source_excess_pressure(
            high_ratios, high_temperature, parameters
        ) - _source_excess_pressure(high_ratios, 300.0, parameters)
        high_targets = (
            (-high_excess - high_static)
            * 1.0e9
            * (high_volume * 6.022e-7 / 4.0)
            / high_energy_difference
        )
        (next_gamma0, next_exponent), weights = _robust_bisquare_gamma_fit(
            ratios,
            np.concatenate((shock_targets, high_targets)),
            fixed_a=parameters["a"],
            start=(start_gamma, 3.0),
            gamma_bounds=(2.0, 4.0),
        )
        current = np.asarray(history[-1])
        following = np.asarray([next_gamma0, parameters["a"], next_exponent])
        history.append(following.tolist())
        stable = stable + 1 if np.all(np.abs(following - current) < 0.2) else 0
        if stable > 5:
            break

    fitted_gamma0, fitted_a, fitted_exponent = np.mean(history[-4:], axis=0)
    fitted_gamma, fitted_theta = _gamma_theta(
        shock_ratios,
        fitted_gamma0,
        fitted_a,
        fitted_exponent,
        parameters["theta0"],
    )
    energy_300 = _debye_energy(fitted_theta, 300.0, parameters["n"])
    excess = _source_excess_pressure(
        shock_ratios, shock_temperature, parameters
    ) - _source_excess_pressure(shock_ratios, 300.0, parameters)
    estimated_energy = energy_300 + (
        (pressure - static - stress - excess)
        * 1.0e9
        * (volume / 4.0 * 0.6022 / 1.0e6)
        / fitted_gamma
    )
    shock_temperature = np.asarray(
        [
            brentq(
                lambda trial: float(
                    _debye_energy(theta_value, trial, parameters["n"])
                    - target_energy
                ),
                1.0,
                1.0e5,
            )
            for theta_value, target_energy in zip(fitted_theta, estimated_energy)
        ]
    )
    excess = _source_excess_pressure(
        shock_ratios, shock_temperature, parameters
    ) - _source_excess_pressure(shock_ratios, 300.0, parameters)
    excess_energy = (
        0.5
        * shock_ratios ** parameters["m"]
        * parameters["optimizer_beta_mass"]
        * (parameters["optimizer_beta_molar_mass"] / 1.0e7)
        * (shock_temperature**2 - 300.0**2)
    )
    shock_residuals = (
        pressure
        - static
        - fitted_gamma
        / density
        * (shock_energy + compression_energy - excess_energy)
        / 1.0e9
        - stress
        - excess
    )
    fitted_high_gamma, fitted_high_theta = _gamma_theta(
        high_ratios,
        fitted_gamma0,
        fitted_a,
        fitted_exponent,
        parameters["theta0"],
    )
    high_energy_difference = _debye_energy(
        fitted_high_theta, high_temperature, parameters["n"]
    ) - _debye_energy(fitted_high_theta, 300.0, parameters["n"])
    high_excess = _source_excess_pressure(
        high_ratios, high_temperature, parameters
    ) - _source_excess_pressure(high_ratios, 300.0, parameters)
    high_residuals = (
        -high_static
        - fitted_high_gamma
        * 1.0e-9
        / (high_volume * 6.022e-7 / 4.0)
        * high_energy_difference
        - high_excess
    )
    pressure_residuals = np.concatenate((shock_residuals, high_residuals))
    return {
        "gamma0": float(fitted_gamma0),
        "a": float(fitted_a),
        "b": float(fitted_exponent),
        "iterations": len(history) - 1,
        "observations": int(ratios.size),
        "shock_observations": int(volume.size),
        "thermal_expansion_observations": int(high_volume.size),
        "nonzero_robust_weights": int(np.count_nonzero(weights)),
        "rmse_pressure_gpa": float(np.sqrt(np.mean(pressure_residuals**2))),
        "max_abs_pressure_residual_gpa": float(
            np.max(np.abs(pressure_residuals))
        ),
        "shock_temperature_range_k": [
            float(np.min(shock_temperature)),
            float(np.max(shock_temperature)),
        ],
    }


def _thermal_refits() -> dict[str, Any]:
    results = {}
    for identifier, parameters in THERMAL_RECORDS.items():
        fit = (
            _refit_mgo(parameters)
            if identifier == "mgo_zhu_2025_pvt"
            else _refit_shock_and_expansion(identifier, parameters)
        )
        target = {name: parameters[name] for name in ("gamma0", "a", "b")}
        differences = {name: fit[name] - value for name, value in target.items()}
        if identifier == "gold_zhu_2025_pvt":
            parity = abs(differences["gamma0"]) < 0.01 and abs(
                differences["b"]
            ) < 0.1
        elif identifier == "platinum_zhu_2025_pvt":
            parity = abs(differences["gamma0"]) < 0.01 and abs(
                differences["b"]
            ) < 0.05
        else:
            parity = abs(differences["gamma0"]) < 0.01 and abs(
                differences["b"]
            ) < 0.005
        assert parity
        results[identifier] = {
            **fit,
            "target": target,
            "differences": differences,
            "parity": parity,
            "fit_datasets": (
                ["zhu_2025_mgo_pvt"]
                if identifier == "mgo_zhu_2025_pvt"
                else [
                    f"zhu_2025_{'au' if identifier.startswith('gold') else 'pt'}_shock",
                    f"zhu_2025_{'au' if identifier.startswith('gold') else 'pt'}_zero_pressure_thermal_expansion",
                ]
            ),
        }
    return results


def _load(dataset_id: str) -> dict[str, np.ndarray]:
    metadata = DATASETS[dataset_id]
    path = metadata["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == metadata["sha256"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == metadata["rows"]
    return {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in rows[0]
        if name != "observation_index"
    }


def _fit_k0_prime(
    target_volume: np.ndarray,
    target: dict[str, float],
    standard_volume: np.ndarray,
    standard: dict[str, float],
) -> dict[str, Any]:
    """Fit one K0' to unweighted pressure differences, as in Ye's figures."""

    def residuals(values: np.ndarray) -> np.ndarray:
        candidate = {**target, "K0_prime": float(values[0])}
        return _vinet_pressure(target_volume, candidate) - _vinet_pressure(
            standard_volume, standard
        )

    result = least_squares(residuals, [target["K0_prime"]])
    fitted = float(result.x[0])
    raw_residuals = np.asarray(result.fun, dtype=float)
    # The paper identifies +/-1 GPa as the random fit scatter. With that as an
    # absolute observation sigma, sqrt((J^T J)^-1) reproduces Table 1's errors.
    standard_error = float(np.sqrt(np.linalg.inv(result.jac.T @ result.jac)[0, 0]))
    pressures = _vinet_pressure(standard_volume, standard)
    return {
        "K0_prime": fitted,
        "K0_prime_standard_error": standard_error,
        "observations": int(target_volume.size),
        "pressure_range_gpa": [float(np.min(pressures)), float(np.max(pressures))],
        "rmse_pressure_difference_gpa": float(np.sqrt(np.mean(raw_residuals**2))),
        "max_abs_pressure_difference_gpa": float(np.max(np.abs(raw_residuals))),
        "solver_success": bool(result.success),
        "solver_message": str(result.message),
    }


def reproduce() -> dict[str, Any]:
    au_mgo = _load("ye_2017_data_s1_au_mgo_300k")
    pt_mgo = _load("ye_2017_data_s2_pt_mgo_300k")

    fits: dict[str, dict[str, Any]] = {}
    fits["Au_from_MgO"] = _fit_k0_prime(
        au_mgo["au_volume_a3"],
        ANCHORS["Au"],
        au_mgo["mgo_volume_a3"],
        ANCHORS["MgO"],
    )
    fits["MgO_from_Au"] = _fit_k0_prime(
        au_mgo["mgo_volume_a3"],
        ANCHORS["MgO"],
        au_mgo["au_volume_a3"],
        ANCHORS["Au"],
    )
    fits["MgO_from_Pt"] = _fit_k0_prime(
        pt_mgo["mgo_volume_a3"],
        ANCHORS["MgO"],
        pt_mgo["pt_volume_a3"],
        ANCHORS["Pt"],
    )

    mgo_from_au = {**ANCHORS["MgO"], "K0_prime": fits["MgO_from_Au"]["K0_prime"]}
    mgo_from_pt = {**ANCHORS["MgO"], "K0_prime": fits["MgO_from_Pt"]["K0_prime"]}
    fits["Pt_from_Au"] = _fit_k0_prime(
        pt_mgo["pt_volume_a3"],
        ANCHORS["Pt"],
        pt_mgo["mgo_volume_a3"],
        mgo_from_au,
    )
    fits["Pt_from_MgO"] = _fit_k0_prime(
        pt_mgo["pt_volume_a3"],
        ANCHORS["Pt"],
        pt_mgo["mgo_volume_a3"],
        ANCHORS["MgO"],
    )
    fits["Au_from_Pt"] = _fit_k0_prime(
        au_mgo["au_volume_a3"],
        ANCHORS["Au"],
        au_mgo["mgo_volume_a3"],
        mgo_from_pt,
    )

    selected_fits = {fit_name for _, fit_name in RECORDS.values()}
    for name, fit in fits.items():
        published = PUBLISHED[name]
        fit["difference_from_published"] = fit["K0_prime"] - published["K0_prime"]
        if name in selected_fits:
            assert f"{fit['K0_prime']:.3f}" == f"{published['K0_prime']:.3f}"
            assert f"{fit['K0_prime_standard_error']:.3f}" == (
                f"{published['error']:.3f}"
            )
        fit["published"] = published

    records = {}
    for identifier, (material_identifier, fit_name) in RECORDS.items():
        document = get_material_document(material_identifier)
        source = next(
            row for row in document["eos_records"] if row["identifier"] == identifier
        )
        fit = fits[fit_name]
        assert source["eos"]["parameters"]["K0_prime"] == PUBLISHED[fit_name][
            "K0_prime"
        ]
        assert source["parameter_errors"]["K0_prime"] == PUBLISHED[fit_name]["error"]
        model = Material.from_eosmat(
            document, record_identifiers=[identifier]
        ).get_eos_record(identifier)
        record_v0 = float(source["eos"]["parameters"]["V0"])
        checkpoints = {
            str(fraction): float(model.pressure(fraction * record_v0))
            for fraction in (0.95, 0.80, 0.65)
        }
        records[identifier] = {"fit": fit_name, "pressure_gpa": checkpoints, **fit}

    thermal_datasets = {}
    for identifier, (filename, digest, rows) in THERMAL_DATASETS.items():
        path = DATASET_DIR / filename
        assert hashlib.sha256(path.read_bytes()).hexdigest() == digest
        with path.open(newline="", encoding="utf-8") as stream:
            observed_rows = sum(1 for _ in csv.DictReader(stream))
        assert observed_rows == rows
        thermal_datasets[identifier] = {
            "path": str(path.relative_to(ROOT)),
            "sha256": digest,
            "rows": rows,
        }

    return {
        "source": {
            "fit_publication_doi": "10.1002/2016JB013811",
            "adopting_preprint_doi": "10.22541/essoar.176236186.65259830/v1",
            "mendeley_release_doi": "10.17632/6kxnhc2g73.3",
            "mendeley_v3_archive_sha256": (
                "77fbd3e3d02c2e5d27bf66b1ea69a88aed46266e1207b32b461fed33f48e4a0a"
            ),
            "mendeley_v3_optimizer_sha256": {
                "Au/Au_opt.m": "0ddb1f5d9814ab12373987770112ded3e0cbf745406369ac83c0241501de8ee6",
                "Pt/Pt_opt.m": "4dc5392a470eea53c9cae698a2082b1f698c38783d4dd62218d70568022ea958",
                "MgO/MgO_opt.m": "cd3e5082f9833c901d7895de3f500551db793c61564992355164b892120e0781",
            },
            "publisher_supplement_sha256": {
                "jgrb52101-sup-0002-supinfo.csv": (
                    "fada21e5fc1a5c9bb1329ad08cc396065e158f215653afdf0978b8bfc2ebcbd9"
                ),
                "jgrb52101-sup-0003-supinfo.csv": (
                    "b48ee2264424d4a67c85070aa07c387d4b45bc41877525fa691edc75c37d757c"
                ),
            },
            "pytheos_revision": {
                "gold.py": "a4fe7240979a8bf5f8e225c13adc15d338ab11f1",
                "periclase.py": "d63b07baca1714e445b7cb29d63108a6b37c5324",
                "platinum.py": "5d774329275f8a5ef9cf70a44a46e4227c796f1c",
                "license": "Apache-2.0",
            },
        },
        "method": {
            "objective": "unweighted pressure differences between co-compressed phases",
            "fixed_parameters": ["V0", "K0"],
            "free_parameters": ["K0_prime"],
            "absolute_pressure_scatter_gpa": 1.0,
            "uncertainty": "sqrt(diag((J^T J)^-1)) for 1 GPa absolute scatter",
        },
        "datasets": {
            identifier: {
                "path": str(metadata["path"].relative_to(ROOT)),
                "sha256": metadata["sha256"],
                "rows": metadata["rows"],
            }
            for identifier, metadata in DATASETS.items()
        },
        "thermal_datasets": thermal_datasets,
        "anchors": ANCHORS,
        "fits": fits,
        "records": records,
        "thermal_records": _thermal_refits(),
        "calculator_diagnostics": _calculator_diagnostics(),
        "thermal_reproduction_scope": {
            "status": "thermal_refit_parity",
            "note": "All released v3 fit-input rows are passed through an independent translation of the source's iterative robust-bisquare optimization. The standalone v3 pressure calculators contain inconsistent fixed constants and are retained only as diagnostics.",
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
