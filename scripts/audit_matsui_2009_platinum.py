"""Replay Matsui (2009) with recovered primary electronic pressure nodes.

The oracle uses SI Vinet/Debye quadrature and manual linear interpolation of
Tsuchiya-Kawamura Table I. Interpolation is an explicit implementation choice;
calculated checkpoints are validation output, never author fitting inputs.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad

from peritheos import get_eos_record
from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import MieGruneisenDebye

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/matsui-2009-platinum-audit.json"
RECORD = "platinum_matsui_2009_vinet_300k"
THERMAL_RECORD = "platinum_matsui_2009_vinet_mgd_electronic"
V0, K0, KP, GAMMA0, Q, THETA0 = 60.38, 273.0, 5.20, 2.70, 1.10, 230.0
FILES = {
    "table1": "platinum-matsui-2009-table1-pvt.csv",
    "table3": "platinum-matsui-2009-table3-grid.csv",
    "expansion": "platinum-arblaster-1997-table2-expansion.csv",
    "holmes": "platinum-holmes-1989-table3-shock.csv",
    "electronic": "tsuchiya-kawamura-2002-table1-electronic-pressure.csv",
}
QUALIFICATION = (
    "Published Vinet/integrated-q Debye pressure and all 51 upstream Pt electronic "
    "nodes are recovered. The 2003 erratum corrects Au only. Electronic pressure "
    "is volume independent and referenced by subtracting the raw 300 K value "
    "0.04 GPa once. Linear interpolation is an explicit numerical implementation "
    "choice, not a specified author rule; extrapolation is forbidden. Rounded "
    "grid agreement verifies the pressure parameterization, not the original "
    "shock/expansion optimization."
)


def rows(key):
    with (DATA / FILES[key]).open() as stream:
        return list(csv.DictReader(stream))


def source_components(volume, temperature):
    """Equations 6-10: cell A^3 -> m^3/mol of atoms, n=1, pressure GPa."""
    volume, temperature = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    if np.any(volume <= 0) or np.any(temperature <= 0):
        raise ValueError("Volumes and temperatures must be positive")
    ratio = volume / V0
    gamma = GAMMA0 * ratio**Q
    theta = THETA0 * np.exp((GAMMA0 - gamma) / Q)
    x = ratio ** (1 / 3)
    cold = 3 * K0 * (1 - x) / x**2 * np.exp(1.5 * (KP - 1) * (1 - x))

    def energy(t, th):
        y = th / t
        integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11)[0]
        return 9 * R * t * integral / y**3

    energy_increment = np.array(
        [energy(t, th) - energy(300, th) for t, th in zip(temperature.flat, theta.flat)]
    ).reshape(volume.shape)
    molar_volume_si = volume * N_A * 1e-30 / 4
    phonon = gamma * energy_increment / molar_volume_si / 1e9
    return cold, phonon, theta, gamma


def electronic_increment(temperature):
    """Independent manual interpolation; no production EOS or np.interp."""
    source = rows("electronic")
    temperatures = np.array([float(r["temperature_k"]) for r in source])
    pressures = np.array([float(r["platinum_electronic_pressure_gpa"]) for r in source])
    target = np.asarray(temperature, dtype=float)
    if not np.all(np.isfinite(target)) or np.any((target < 0) | (target > 5000)):
        raise ValueError("Electronic table covers 0-5000 K")
    values = []
    for t in target.flat:
        j = min(int(t // 100), len(temperatures) - 2)
        fraction = (t - temperatures[j]) / (temperatures[j + 1] - temperatures[j])
        values.append(
            pressures[j] + fraction * (pressures[j + 1] - pressures[j]) - pressures[3]
        )
    return np.asarray(values).reshape(target.shape)


def ledger_outcome(record):
    """Keep equation reconstruction separate from original author fit parity."""
    result = audit()
    identifiers = record["scientific_validation"]["primary_data_check"][
        "dataset_identifiers"
    ]
    return {
        "status": "source_reconstruction"
        if record["identifier"] == THERMAL_RECORD
        else "not_refittable",
        "dataset_identifiers": identifiers,
        "observations": 35 if record["identifier"] == THERMAL_RECORD else 0,
        "reconstruction_kind": "published_equation_with_tabulated_electronic_input",
        "coefficient_optimization_performed": False,
        "original_author_fit_status": "not_reproduced",
        "qualification": QUALIFICATION,
        "reason": "Marsh inputs, original Hugoniot/expansion sampling, weights and joint-fit objective remain unavailable. The pressure parameterization is reproduced at source nodes; sub-grid linear interpolation is an implementation choice.",
        "reproduction": result,
    }


def audit():
    cold_record = get_eos_record(RECORD)
    thermal_record = get_eos_record(THERMAL_RECORD)
    # Reuse the production kernel only as the comparator, never as the oracle.
    # An isothermal catalog record keeps cell volumes internally; thermal
    # kernels require J/bar/mol. Convert both the Vinet V0 and evaluated V.
    molar_scale = N_A * 1e-25 / 4
    debye = MieGruneisenDebye(
        Vinet(V0 * molar_scale, K0, KP), 300, THETA0, GAMMA0, Q, 1
    )
    ratio = np.array([float(r["volume_ratio"]) for r in rows("table3")])
    temperature = np.array([300, 500, 1000, 2000, 3000])
    cold, phonon, theta, gamma = source_components(V0 * ratio[:, None], temperature)
    published = np.array(
        [[float(r[f"pressure_{t}k_gpa"]) for t in temperature] for r in rows("table3")]
    )
    full = cold + phonon + electronic_increment(temperature)
    electronic = published - cold - phonon
    # Each printed pressure has a +/-0.005 GPa rounding interval. Intersection
    # tests the volume-independent approximation; it does not fit a function.
    inferred = []
    for j, t in enumerate(temperature):
        lo = float(np.max(electronic[:, j] - 0.005))
        hi = float(np.min(electronic[:, j] + 0.005))
        inferred.append(
            {
                "temperature_k": int(t),
                "electronic_increment_gpa_by_volume": electronic[:, j].tolist(),
                "common_rounding_interval_gpa": [lo, hi] if lo <= hi else None,
                "volume_spread_gpa": float(np.ptp(electronic[:, j])),
            }
        )
    electronic_examples = {1000: 0.21, 3000: 1.60, 5000: 3.78}
    known_grid_checks = {}
    for t in (1000, 3000):
        j = list(temperature).index(t)
        residual = cold[:, j] + phonon[:, j] + electronic_examples[t] - published[:, j]
        known_grid_checks[str(t)] = {
            "reference_subtracted_electronic_gpa": electronic_examples[t],
            "max_abs_pressure_difference_gpa": float(np.max(abs(residual))),
            "pressure_difference_gpa_by_volume": residual.tolist(),
        }
    table1 = []
    for row in rows("table1"):
        t, r = float(row["temperature_k"]), float(row["volume_ratio"])
        c, ph, _, _ = source_components(V0 * r, t)
        obs, calc = (
            float(row["observed_pressure_gpa"]),
            float(row["calculated_pressure_gpa"]),
        )
        table1.append(
            {
                "temperature_k": t,
                "volume_ratio": r,
                "observed_mgo_pressure_gpa": obs,
                "printed_calculated_pt_pressure_gpa": calc,
                "cold_plus_phonon_pressure_gpa": float(c + ph),
                "electronic_increment_inferred_from_output_gpa": float(calc - c - ph),
                "observed_minus_printed_calculated_gpa": obs - calc,
                "reconstructed_pressure_gpa": float(c + ph + electronic_increment(t)),
                "reconstructed_minus_printed_calculated_gpa": float(
                    c + ph + electronic_increment(t) - calc
                ),
                "observed_minus_reconstructed_gpa": float(
                    obs - c - ph - electronic_increment(t)
                ),
            }
        )
    shock = []
    for row in rows("holmes"):
        if float(row["pressure_gpa"]) > 290:
            continue
        up, us = float(row["mass_velocity_km_s"]), float(row["shock_velocity_km_s"])
        rho = float(row["sample_density_g_cm3"])
        shock.append(
            {
                "shot": row["shot"],
                "published_pressure_gpa": float(row["pressure_gpa"]),
                "momentum_pressure_gpa": rho * us * up,
                "mass_conservation_volume_ratio": 1 - up / us,
                "published_density_volume_ratio": rho / float(row["density_g_cm3"]),
                "us_minus_matsui_linear_fit_km_s": us - (3.604 + 1.543 * up),
            }
        )
    expansion = []
    a300 = next(
        float(r["lattice_parameter_nm"])
        for r in rows("expansion")
        if float(r["temperature_k"]) == 300
    )
    for row in rows("expansion"):
        r = (float(row["lattice_parameter_nm"]) / a300) ** 3
        c, ph, _, _ = source_components(V0 * r, float(row["temperature_k"]))
        expansion.append(
            {
                "temperature_k": float(row["temperature_k"]),
                "volume_ratio_from_rounded_lattice_parameter": r,
                "cold_plus_phonon_pressure_gpa": float(c + ph),
                "electronic_increment_required_for_zero_pressure_gpa": float(-c - ph),
                "zero_pressure_residual_with_linear_electronic_gpa": float(
                    c + ph + electronic_increment(float(row["temperature_k"]))
                ),
            }
        )
    kernel_phonon = debye.thermal_pressure(
        V0 * ratio[:, None] * molar_scale, temperature
    )
    kernel_cold = cold_record.pressure(V0 * ratio[:, None])
    residuals = np.array([r["observed_minus_printed_calculated_gpa"] for r in table1])
    return {
        "doi": "10.1063/1.3054331",
        "audit_date": "2026-10-09",
        "record_identifier": RECORD,
        "scope": QUALIFICATION,
        "continuous_thermal_record_status": "registered_with_explicit_linear_interpolation",
        "thermal_record_identifier": THERMAL_RECORD,
        "electronic_source": {
            "doi": "10.1103/PhysRevB.66.094115",
            "erratum_doi": "10.1103/PhysRevB.67.019902",
            "erratum_changes_platinum": False,
            "raw_pressure_reference": "0 K",
            "raw_300k_pressure_gpa": 0.04,
            "nodes": 51,
            "temperature_range_k": [0, 5000],
            "interpolation": "linear, implementation choice not specified by authors",
            "extrapolation": "forbidden",
            "original_analytic_function_or_DOS": "not published numerically",
            "uncertainties": None,
        },
        "author_fit_reproduction_status": "not_reproduced",
        "input_sha256": {
            f: hashlib.sha256((DATA / f).read_bytes()).hexdigest()
            for f in FILES.values()
        },
        "parameters": {
            "V0_cell_a3": V0,
            "K0_gpa": K0,
            "K0_prime": KP,
            "gamma0": GAMMA0,
            "q": Q,
            "theta0_k": THETA0,
            "n": 1,
            "Tr_k": 300,
        },
        "published_parameter_errors": None,
        "published_parameter_covariance": None,
        "independent_cold_max_difference_gpa": float(np.max(abs(kernel_cold - cold))),
        "independent_phonon_max_difference_gpa": float(
            np.max(abs(kernel_phonon - phonon))
        ),
        "independent_full_max_difference_gpa": float(
            np.max(
                abs(thermal_record.pressure(V0 * ratio[:, None], temperature) - full)
            )
        ),
        "table3": {
            "volume_ratio": ratio.tolist(),
            "temperature_k": temperature.tolist(),
            "published_pressure_gpa": published.tolist(),
            "reconstructed_pressure_gpa": full.tolist(),
            "reconstructed_minus_published_gpa": (full - published).tolist(),
            "max_abs_pressure_difference_gpa": float(np.max(abs(full - published))),
            "rmse_gpa": float(np.sqrt(np.mean((full - published) ** 2))),
            "reference_subtracted_electronic_gpa": electronic_increment(
                temperature
            ).tolist(),
            "phonon_increment_gpa": phonon.tolist(),
            "theta_k": theta[:, 0].tolist(),
            "gamma": gamma[:, 0].tolist(),
            "inferred_electronic": inferred,
        },
        "electronic_examples_reference_subtracted_gpa": {
            str(t): p for t, p in electronic_examples.items()
        },
        "quadratic_coefficient_implied_by_each_example_gpa_per_k2": {
            str(t): p / (t * t - 300**2) for t, p in electronic_examples.items()
        },
        "known_electronic_table3_checks": known_grid_checks,
        "table1_validation": table1,
        "table1_max_abs_calculated_reconstruction_difference_gpa": float(
            max(abs(r["reconstructed_minus_printed_calculated_gpa"]) for r in table1)
        ),
        "table1_printed_residual_rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
        "table1_printed_max_abs_residual_gpa": float(np.max(abs(residuals))),
        "holmes_subset": shock,
        "shock_scope": "Three of the seven primary Holmes Table III shots are below 290 GPa. The other four are excluded, not relabeled as static or fitted to 300 K. Marsh rows and Matsui's exact Hugoniot sampling/weights remain unavailable.",
        "arblaster_expansion": expansion,
        "expansion_scope": "32 original 1997 Table II assessed crystallographic states at 100-2000 K, normalized via rounded a(T)/a(300). These are literature-assessed expansion values, not Matsui measurements. The 2006 methodology/erratum is separate; Matsui cites the 1997 assessment. The recovered electronic table is linearly interpolated for expansion diagnostics; this does not recover the original author selection or objective.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = json.dumps(audit(), indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.check:
        if OUTPUT.read_text() != payload:
            raise SystemExit(
                "Matsui audit is stale; run python -m scripts.audit_matsui_2009_platinum"
            )
    else:
        OUTPUT.write_text(payload)


if __name__ == "__main__":
    main()
