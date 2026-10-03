"""Independent BM3/Vinet/Debye checks and conditional fits of Miozzi (2020)."""

from __future__ import annotations

import argparse
import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/miozzi-2020-iron-reproduction.json"
RECORDS = {
    "iron_miozzi_2020_bm3": [22.80, 129.0, 6.2],
    "iron_miozzi_2020_vinet": [22.81, 125.0, 6.5],
    "iron_miozzi_2020_bm3_mgd": [22.81, 129.0, 6.24, 1.11, 0.3],
}
DATASETS = ["iron_miozzi_2020_he_pvt", "iron_miozzi_2020_mgo_pvt"]
QUALIFICATION = (
    "Conditional diagnostic, not the original EosFit-7c objective: the source "
    "does not specify its numerical weights, covariance, iteration settings, "
    "or handling of blank/zero coordinate errors. Missing errors contribute no "
    "term to diagnostic variance; one cold row with no positive effective error "
    "is excluded only from that diagnostic. Printed errors are provisionally "
    "treated as standard deviations solely for this diagnostic. Published "
    "coefficients remain "
    "unchanged. Thermal V0=22.81 A^3 follows Section 3.3; Table 1 instead gives "
    "22.80(2) A^3 and the printed molar value converts to 22.81581 A^3. "
    "The thermal source replay now uses n=2 with volume per mole Fe, as "
    "confirmed in user-relayed personal communication with Miozzi, recorded "
    "2026-10-03. This diagnostic fixes V0 even though the author clarified "
    "that V0 was free in the final stage. "
    "See literature-reproductions/miozzi-2020-iron.md."
)


def observations():
    """Read all source rows; missing errors stay NaN, printed zeros stay zero."""
    rows = []
    for medium in ("he", "mgo"):
        path = ROOT / f"peritheos/data/datasets/iron-miozzi-2020-{medium}-pvt.csv"
        with path.open(encoding="utf-8") as stream:
            rows.extend(csv.DictReader(stream))
    return {
        name: np.array([float(r[name]) if r[name] else np.nan for r in rows])
        for name in (
            "pressure_gpa",
            "pressure_error_gpa",
            "temperature_k",
            "temperature_error_k",
            "volume_a3",
            "volume_error_a3",
        )
    }


_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(64)


def source_pressure(volume, temperature, parameters, *, vinet=False, n=2):
    """Author-reported MGD setup, independent of Peritheos.

    Eq. 5's denominator and comma, and Eq. 7's missing cube/integration limit,
    are documented typographical errors. Cell V contains two Fe atoms.
    Thermal volume is per mole Fe; n=2 reproduces the personal-communication
    setting. Pass n=1 explicitly for a physically normalized per-Fe model.
    """
    v0, k0, kp = parameters[:3]
    volume, temperature = np.broadcast_arrays(volume, temperature)
    x = (v0 / volume) ** (1 / 3)
    if vinet:
        eta = 1 / x
        cold = 3 * k0 * (1 - eta) / eta**2 * np.exp(1.5 * (kp - 1) * (1 - eta))
    else:
        cold = 1.5 * k0 * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))
    if len(parameters) == 3:
        return cold
    gamma0, q = parameters[3:]
    gamma = gamma0 * (volume / v0) ** q
    theta = 420 * np.exp((gamma0 - gamma) / q)

    def energy(temp):
        limit = theta / temp
        z = limit[..., None] * (_NODES + 1) / 2
        integral = limit / 2 * np.sum(_WEIGHTS * z**3 / np.expm1(z), axis=-1)
        return 9 * n * 8.31446261815324 * temp * integral / limit**3

    # J/mol -> GPa using molar cm^3/mol = V_cell*N_A/2e24.
    thermal = gamma * (energy(temperature) - energy(np.full_like(temperature, 300)))
    return cold + thermal / (volume * 6.02214076e23 / 2e24) / 1000


@lru_cache(maxsize=1)
def reproduce():
    data = observations()
    output = {}
    for identifier, published in RECORDS.items():
        thermal = len(published) == 5
        select = (
            np.ones(len(data["pressure_gpa"]), dtype=bool)
            if thermal
            else data["temperature_k"] <= 300
        )
        p, v, t, pe, ve, te = (
            data[n][select]
            for n in (
                "pressure_gpa",
                "volume_a3",
                "temperature_k",
                "pressure_error_gpa",
                "volume_error_a3",
                "temperature_error_k",
            )
        )
        vinet = "vinet" in identifier

        def prediction(values, volume=v, temperature=t):
            return source_pressure(volume, temperature, values, vinet=vinet)

        dv = (prediction(published, v + 1e-5) - prediction(published, v - 1e-5)) / 2e-5
        dt = (
            prediction(published, v, t + 1e-3) - prediction(published, v, t - 1e-3)
        ) / 2e-3
        effective = np.sqrt(
            np.nan_to_num(pe) ** 2
            + (dv * np.nan_to_num(ve)) ** 2
            + (dt * np.nan_to_num(te)) ** 2
        )
        fits = {}
        for mode, errors in [
            ("unweighted", np.ones(len(p))),
            ("effective_errors", effective),
        ]:
            # Zero or absent source errors cannot be infinite-weight constraints.
            usable = errors > 0
            denominator = errors[usable]
            scales = np.array([100, 5, 1, 1] if thermal else [23, 100, 5])
            initial = np.array(published[1:] if thermal else published)

            def full(z):
                return np.r_[published[0], z * scales] if thermal else z * scales

            fit = least_squares(
                lambda z: (prediction(full(z))[usable] - p[usable]) / denominator,
                initial / scales,
                bounds=([0.1] * len(scales), [10] * len(scales)),
                xtol=1e-12,
                ftol=1e-12,
                gtol=1e-12,
                max_nfev=3000,
            )
            fitted = full(fit.x)
            dof = int(np.sum(usable)) - len(scales)
            covariance_scaled = (
                np.linalg.inv(fit.jac.T @ fit.jac) * np.sum(fit.fun**2) / dof
            )
            covariance = covariance_scaled * scales[:, None] * scales[None, :]
            fits[mode] = {
                "row_count": int(np.sum(usable)),
                "parameters": fitted.tolist(),
                "free_parameter_names": ["K0", "K0_prime", "gamma0", "q"]
                if thermal
                else ["V0", "K0", "K0_prime"],
                "standard_errors": np.sqrt(np.diag(covariance)).tolist(),
                "covariance": covariance.tolist(),
                "degrees_of_freedom": dof,
                "rmse_gpa": float(np.sqrt(np.mean((prediction(fitted) - p) ** 2))),
                "solver_success": bool(fit.success),
            }
        residual = prediction(published) - p
        output[identifier] = {
            "row_count": len(p),
            "published_parameters": published,
            "published_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "published_residual_range_gpa": [
                float(min(residual)),
                float(max(residual)),
            ],
            "fits": fits,
            "qualification": QUALIFICATION,
        }
    return output


def ledger_outcome(record):
    result = reproduce()[record["identifier"]]
    thermal = "thermal" in record
    fit = result["fits"]["effective_errors"]
    fitted = fit["parameters"]
    names = ["V0", "K0", "K0_prime"] + (["gamma0", "q"] if thermal else [])
    errors = dict(record["parameter_errors"])
    if thermal:
        errors.update(record["thermal"]["parameter_errors"])
    fit_errors = dict(zip(fit["free_parameter_names"], fit["standard_errors"]))
    comparisons = [
        {
            "parameter": n,
            "published": p,
            "refit": f,
            "difference": f - p,
            "published_error": errors[n],
            "refit_error": fit_errors[n],
            "relative_difference": abs(f - p) / abs(p),
            "within_combined_2sigma": None,
            "within_reported_error": None
            if errors[n] is None
            else abs(f - p) <= errors[n],
            "similar": abs(f - p) <= 0.05 * abs(p),
        }
        for n, p, f in zip(names, result["published_parameters"], fitted)
        if n not in (record.get("fixed_parameters", []) + (["V0"] if thermal else []))
    ]
    return {
        "status": "similar"
        if all(
            c["similar"] and c["within_reported_error"] is not False
            for c in comparisons
        )
        else "parity_not_achieved",
        "fit_kind": "conditional_effective_coordinate_errors",
        "dataset_identifiers": DATASETS,
        "observations": fit["row_count"],
        "selection": "All 131 rows for thermal; 35 of 36 room-temperature rows for the effective-error isothermal diagnostic (one lacks positive effective error).",
        "parameters": comparisons,
        "rmse_pressure_gpa": fit["rmse_gpa"],
        "solver_success": fit["solver_success"],
        "qualification": QUALIFICATION,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = json.dumps(reproduce(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.check:
        if OUTPUT.read_text(encoding="utf-8") != text:
            raise SystemExit("Miozzi reproduction output is stale")
    else:
        OUTPUT.write_text(text, encoding="utf-8")
