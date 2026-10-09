"""Independently reconstruct Fei (2007) Pt and audit its joint fit inputs."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos import get_eos_record

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/fei-2007-platinum-reproduction.json"
RECORD = "platinum_fei_2007_vinet_mgd"
PARAMETERS = {
    "gold": (67.85, 167.0, 6.0, 2.97, 0.6, 170.0),
    "platinum": (60.38, 277.0, 5.08, 2.72, 0.5, 230.0),
}
DATASETS = {
    "platinum_dewaele_2004_table1_compression": "platinum-dewaele-2004-table1-compression.csv",
    "platinum_fei_2004_table2": "platinum-fei-2004-table2.csv",
}
QUALIFICATION = (
    "All 36 Dewaele (2004) cold rows and all 42 paired Fei (2004) Au-Pt rows "
    "are used. Au pressures are independently recalculated with Fei (2007) "
    "Table 1 and its printed variable-exponent Debye law. Dewaele's raw 298 K "
    "provenance remains unchanged; the diagnostic assigns that subset to the "
    "source's 300 K reference branch. Equal pressure weights, exact SI R and "
    "Avogadro constant, fixed V0/K0/theta0, and fitted K0-prime/gamma0/q are "
    "explicit diagnostic choices. The source does not supply its weights, "
    "complete fitted/fixed parameter list, covariance, or unrounded inputs. "
    "Diagnostic coefficients do not replace the published record."
)


def source_pressure(volume, temperature, material="platinum", coefficients=None):
    """Table 1 and Equations 2-3, independent of Peritheos EOS evaluation.

    Volumes are conventional fcc cell A^3 (four atoms). Energy is J/mol of
    atoms. The direct variable-exponent theta law is the printed Fei choice.
    """
    v0, k0, kp, gamma0, q, theta0 = (
        PARAMETERS[material] if coefficients is None else coefficients
    )
    volume, temperature = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    ratio = volume / v0
    gamma = gamma0 * ratio**q
    theta = theta0 * ratio ** (-gamma)
    x = ratio ** (1.0 / 3.0)
    cold = 3 * k0 * (1 - x) * np.exp(1.5 * (kp - 1) * (1 - x)) / x**2

    def energy(t, th):
        y = th / t
        integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11, epsrel=1e-11)[
            0
        ]
        return 9 * R * t * integral / y**3

    increment = np.array(
        [
            energy(t, th) - energy(300.0, th)
            for t, th in zip(temperature.flat, theta.flat)
        ]
    ).reshape(volume.shape)
    molar_volume_si = volume * N_A * 1e-30 / 4
    result = cold + gamma * increment / molar_volume_si / 1e9
    return float(result) if result.ndim == 0 else result


def load_rows(identifier):
    with (DATA / DATASETS[identifier]).open() as stream:
        return list(csv.DictReader(stream))


def observations():
    """Recover pressures from measured Au volumes, never from Pt itself."""
    cold = load_rows("platinum_dewaele_2004_table1_compression")
    paired = load_rows("platinum_fei_2004_table2")
    hot_volume = np.array([float(r["volume_a3"]) for r in paired])
    temperature = np.array([float(r["temperature_k"]) for r in paired])
    gold_volume = np.array([float(r["gold_volume_a3"]) for r in paired])
    reduced_pressure = source_pressure(gold_volume, temperature, "gold")
    volume = np.r_[[4 * float(r["atomic_volume_a3"]) for r in cold], hot_volume]
    temperatures = np.r_[np.full(len(cold), 300.0), temperature]
    pressure = np.r_[
        [float(r["ruby_pressure_revised_gpa"]) for r in cold], reduced_pressure
    ]
    return volume, temperatures, pressure, paired, reduced_pressure


def diagnostic_fit(volume, temperature, pressure, *, free_v0=False):
    """Qualified pressure-residual fit; the author's objective is unknown."""
    base = np.array(PARAMETERS["platinum"])
    free = [0, 2, 3, 4] if free_v0 else [2, 3, 4]
    scales = np.array([60.0, 200.0, 1.0, 1.0, 1.0, 200.0])[free]
    names = ["V0", "K0", "K0_prime", "gamma0", "q", "theta0"]

    def unpack(x):
        c = base.copy()
        c[free] = x * scales
        return c

    fit = least_squares(
        lambda x: (
            source_pressure(volume, temperature, coefficients=unpack(x)) - pressure
        ),
        base[free] / scales,
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=2000,
    )
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * (
        fit.fun @ fit.fun / (len(volume) - len(free))
    )
    errors = np.sqrt(np.diag(covariance)) * scales
    return {
        "observations": len(volume),
        "free_parameters": [names[i] for i in free],
        "fixed_parameters": [names[i] for i in range(6) if i not in free],
        "parameters": dict(zip([names[i] for i in free], unpack(fit.x)[free])),
        "diagnostic_standard_errors": dict(zip([names[i] for i in free], errors)),
        "rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        "published_rmse_gpa": float(
            np.sqrt(np.mean((source_pressure(volume, temperature) - pressure) ** 2))
        ),
        "solver_success": bool(fit.success),
        "qualification": QUALIFICATION,
    }


@lru_cache(maxsize=1)
def reproduce():
    volume, temperature, pressure, paired, reduced = observations()
    record = get_eos_record(RECORD)
    grid_volume = np.array([60.38, 0.9 * 60.38, 0.8 * 60.38])[:, None]
    grid_temperature = np.array([300.0, 1000.0, 1873.0])[None, :]
    reference = source_pressure(grid_volume, grid_temperature)
    native = record.pressure(grid_volume, grid_temperature)
    gold_volume = np.array([float(r["gold_volume_a3"]) for r in paired])
    paired_temperature = np.array([float(r["temperature_k"]) for r in paired])
    gold_native = get_eos_record("gold_fei_2007_vinet_2").pressure(
        gold_volume, paired_temperature
    )
    return {
        "doi": "10.1073/pnas.0609013104",
        "primary_source": "https://pmc.ncbi.nlm.nih.gov/articles/PMC1890468/",
        "record_identifier": RECORD,
        "scope": QUALIFICATION,
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in DATASETS.values()
        },
        "observations": 78,
        "cold_rows": 36,
        "paired_rows": 42,
        "independent_equation_max_difference_gpa": float(
            np.max(abs(native - reference))
        ),
        "independent_au_reduction_max_difference_gpa": float(
            np.max(abs(gold_native - reduced))
        ),
        "equation_grid": {
            "volume_a3": grid_volume[:, 0].tolist(),
            "temperature_k": grid_temperature[0].tolist(),
            "pressure_gpa": reference.tolist(),
        },
        "paired_pressure_reductions": [
            {
                "run": row["run"],
                "temperature_k": float(row["temperature_k"]),
                "platinum_volume_a3": float(row["volume_a3"]),
                "gold_volume_a3": float(row["gold_volume_a3"]),
                "original_fei_2004_pressure_gpa": float(row["pressure_gpa"]),
                "reconstructed_fei_2007_au_pressure_gpa": float(p),
            }
            for row, p in zip(paired, reduced)
        ],
        "fixed_v0_joint_diagnostic": diagnostic_fit(volume, temperature, pressure),
        "free_v0_joint_sensitivity": diagnostic_fit(
            volume, temperature, pressure, free_v0=True
        ),
        "figure2_temperature_selection_sensitivity": diagnostic_fit(
            volume[temperature != 1673],
            temperature[temperature != 1673],
            pressure[temperature != 1673],
        ),
        "fit_reproduction_status": "parity_not_achieved",
    }


def ledger_outcome(record):
    result = reproduce()
    fit = result["fixed_v0_joint_diagnostic"]
    published = {"K0_prime": 5.08, "gamma0": 2.72, "q": 0.5}
    errors = {"K0_prime": 0.02, "gamma0": 0.03, "q": 0.5}
    return {
        "status": "parity_not_achieved",
        "fit_kind": "diagnostic_joint_au_rereduced_pvt",
        "dataset_identifiers": list(DATASETS),
        "observations": 78,
        "free_parameters": fit["free_parameters"],
        "fixed_parameters": fit["fixed_parameters"],
        "parameters": [
            {
                "parameter": name,
                "published": p,
                "published_error": errors[name],
                "refit": fit["parameters"][name],
                "refit_error": fit["diagnostic_standard_errors"][name],
                "relative_difference": abs(fit["parameters"][name] - p) / abs(p),
                "similar": bool(abs(fit["parameters"][name] - p) <= errors[name]),
                "within_combined_2sigma": None,
            }
            for name, p in published.items()
        ],
        "rmse_gpa": fit["rmse_gpa"],
        "published_rmse_gpa": fit["published_rmse_gpa"],
        "solver_success": fit["solver_success"],
        "reason": QUALIFICATION,
        "reproduction": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = json.dumps(reproduce(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text(encoding="utf-8") != payload:
            raise SystemExit("Fei (2007) Pt reproduction is stale")
    else:
        OUTPUT.write_text(payload, encoding="utf-8")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
