"""Audit Yokoo (2009) Au fits separately from full thermal EOS checkpoints."""

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

from scripts.reproduce_yokoo_2009_platinum import check_saved

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/yokoo-2009-gold-reproduction.json"
SOURCE = "gold-yokoo-2009-source.json"
TABLE = "gold-yokoo-2009-table3-isochores.csv"
ELECTRONIC = "tsuchiya-kawamura-2002-table1-electronic-pressure.csv"
V0 = 67.72  # four-atom fcc cell, converted from rho0=19.32 g/cm^3
RECORDS = {
    "gold_yokoo_2009_bm3_300k": ("BM3", 5.79),
    "gold_yokoo_2009_vinet_300k": ("Vinet", 5.94),
}


def source_bm3(ratio, reference_ratio=1.0, k0=167.5, kp=5.79):
    """Equation (2), with explicitly chosen reference state and coefficients."""
    z = np.asarray(ratio, dtype=float) / reference_ratio
    return (
        1.5
        * k0
        * (z ** (-7 / 3) - z ** (-5 / 3))
        * (1 + 0.75 * (kp - 4) * (z ** (-2 / 3) - 1))
    )


def source_vinet(ratio):
    """Table II footnote a: the separate fitted 300 K Vinet branch."""
    x = np.asarray(ratio, dtype=float) ** (1 / 3)
    return 3 * 167.5 * (1 - x) / x**2 * np.exp(1.5 * (5.94 - 1) * (1 - x))


def gamma_theta(ratio):
    """Equations (10) and (12), using Au coefficients and ambient V0."""
    ratio = np.asarray(ratio, dtype=float)
    gamma = 2.96 * (1 + 0.45 * (ratio**4.2 - 1))
    theta = 170 * ratio ** (-(1 - 0.45) * 2.96) * np.exp(-(gamma - 2.96) / 4.2)
    return gamma, theta


def phonon_pressure(ratio, temperature):
    """Absolute phonon pressure; diagnostic only, not the full thermal EOS."""
    if temperature == 0:
        return 0.0
    gamma, theta = gamma_theta(ratio)
    x = float(theta) / temperature
    integral = quad(lambda z: z**3 / np.expm1(z), 0, x, epsabs=1e-11)[0]
    energy = 9 * R * temperature * integral / x**3
    molar_volume_si = V0 * ratio * N_A * 1e-30 / 4
    return float(gamma * energy / molar_volume_si / 1e9)


def table_rows():
    with (DATA / TABLE).open() as stream:
        return list(csv.DictReader(stream))


def shock_relation_checkpoints():
    """Derived Eq. (13)-(15) outputs, never represented as measured shots."""
    checkpoints = []
    for up in [0.0, 1.0, 2.0, 3.0, 3.5]:
        us = 2.995 + 1.653 * up - 0.013 * up**2
        checkpoints.append(
            {
                "particle_velocity_km_s": up,
                "shock_velocity_km_s": us,
                "pressure_increment_gpa": 19.32 * us * up,
                "volume_ratio": 1 - up / us,
                "kind": "derived_shock_regression_output_not_observation",
            }
        )
    return checkpoints


@lru_cache(maxsize=1)
def reproduce():
    from peritheos import get_eos_record
    from scripts.reproduce_yokoo_2009_platinum import reproduce as platinum_audit

    rows = table_rows()
    cold = [row for row in rows if float(row["temperature_k"]) == 0]
    ratios = np.array([float(row["volume_ratio"]) for row in cold])
    p0 = np.array([float(row["pressure_gpa"]) for row in cold])
    p300 = np.array(
        [
            float(row["pressure_gpa"])
            for row in rows
            if float(row["temperature_k"]) == 300
        ]
    )
    fit = least_squares(
        lambda vc: source_bm3(ratios, vc[0], 180.0, 5.61) - p0,
        [0.990],
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
    )
    grid = np.array([1.0, 0.98, 0.9, 0.8, 0.7, 0.6])
    branches = {}
    for identifier, (kind, _) in RECORDS.items():
        equation = source_bm3 if kind == "BM3" else source_vinet
        expected = equation(grid)
        record = get_eos_record(identifier)
        native = record.pressure(grid * V0)
        differences = equation(ratios) - p300
        branches[identifier] = {
            "volume_ratio": grid.tolist(),
            "pressure_gpa": expected.tolist(),
            "independent_equation_max_difference_gpa": float(
                np.max(abs(native - expected))
            ),
            "inversion_max_difference_a3": float(
                np.max(abs(record.volume(expected) - grid * V0))
            ),
            "fit_minus_full_model_table3_300k_gpa": differences.tolist(),
            "max_abs_fit_minus_full_model_table3_300k_gpa": float(
                np.max(abs(differences))
            ),
            "pressure_at_volume_ratio_0p6_gpa": float(equation(0.6)),
        }
    with (DATA / ELECTRONIC).open() as stream:
        electronic = list(csv.DictReader(stream))
    electronic_by_t = {
        float(row["temperature_k"]): float(
            row["gold_corrected_electronic_pressure_gpa"]
        )
        for row in electronic
    }
    thermal = []
    for row in rows:
        ratio = float(row["volume_ratio"])
        t = float(row["temperature_k"])
        phonon = phonon_pressure(ratio, t)
        thermal.append(
            {
                **row,
                "phonon_pressure_gpa": phonon,
                "table_total_minus_0k_minus_phonon_gpa": float(row["pressure_gpa"])
                - p0[np.argmin(abs(ratios - ratio))]
                - phonon,
                "upstream_electronic_pressure_at_source_node_gpa": electronic_by_t[t],
                "table_minus_0k_phonon_electronic_gpa": float(row["pressure_gpa"])
                - p0[np.argmin(abs(ratios - ratio))]
                - phonon
                - electronic_by_t[t],
            }
        )
    pairs = []
    for row in platinum_audit()["dewaele_simultaneous_volume_pairs"]:
        au_ratio = 4 * row["gold_atomic_volume_a3"] / V0
        pairs.append(
            {
                **row,
                "derived_yokoo_vinet_au_pressure_gpa": float(source_vinet(au_ratio)),
                "derived_yokoo_bm3_au_pressure_gpa": float(source_bm3(au_ratio)),
                "fitted_vinet_pt_minus_fitted_vinet_au_gpa": row[
                    "derived_yokoo_vinet_pt_pressure_gpa"
                ]
                - float(source_vinet(au_ratio)),
                "scope": "Comparison of fitted isotherms, not exact full-model Figure 6 replay.",
            }
        )
    files = [
        TABLE,
        SOURCE,
        "gold-dewaele-2004-table1-compression.csv",
        "platinum-dewaele-2004-table1-compression.csv",
        ELECTRONIC,
    ]
    return {
        "doi": "10.1103/PhysRevB.80.104114",
        "record_identifiers": list(RECORDS),
        "author_fit_status": "not_refittable",
        "qualification": "Published 300 K BM3/Vinet fits independently implemented. Table III contains derived full-model outputs, not measurements. The inferred cold volume is a circular checkpoint diagnostic, not an original author refit. Corrected upstream electronic pressure nodes are recovered, but exact Yokoo cold/phonon normalization, adopted electronic implementation and numerical caloric/raw optimization inputs remain unresolved.",
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in files
        },
        "density_derived_V0_cell_a3": 4 * 196.96657 / 19.32 / N_A * 1e24,
        "catalog_V0_cell_a3": V0,
        "branches": branches,
        "cold_checkpoint_diagnostic": {
            "Vc_over_V0": float(fit.x[0]),
            "K0_fixed_gpa": 180.0,
            "K0_prime_fixed": 5.61,
            "rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
            "max_abs_residual_gpa": float(np.max(abs(fit.fun))),
            "kind": "derived_output_only_not_independent_fit",
        },
        "gamma_theta_checkpoints": [
            {
                "volume_ratio": float(r),
                "gamma": float(gamma_theta(r)[0]),
                "theta_k": float(gamma_theta(r)[1]),
            }
            for r in [1, 0.8, 0.6]
        ],
        "table3_thermal_term_diagnostic": thermal,
        "max_abs_table_minus_0k_phonon_electronic_gpa": max(
            abs(row["table_minus_0k_phonon_electronic_gpa"]) for row in thermal
        ),
        "electronic_source_nodes": electronic,
        "shock_relation_checkpoints": shock_relation_checkpoints(),
        "dewaele_simultaneous_volume_pairs": pairs,
    }


def ledger_outcome(record):
    return {
        "record_identifier": record["identifier"],
        "status": "not_refittable",
        "fit_kind": "published_reference_isotherm_from_joint_shock_thermodynamic_model",
        "observations": 0,
        "dataset_identifiers": [
            "gold_yokoo_2009_table3_isochores",
            "tsuchiya_kawamura_2002_table1_electronic_pressure",
        ],
        "reason": "Table II publishes separate 300 K BM3/Vinet fits. Table III is derived full thermal model output, not independent observations. Corrected upstream electronic pressure nodes are recovered. Complete Au shock/ambient selections, numerical electronic energy/DOS, original objective/weights and exact cold/phonon normalization remain unavailable. See [Au audit](literature-reproductions/yokoo-2009-gold.md).",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    assert len(result["table3_thermal_term_diagnostic"]) == 162
    assert len(result["dewaele_simultaneous_volume_pairs"]) == 36
    for branch in result["branches"].values():
        assert branch["independent_equation_max_difference_gpa"] < 1e-9
        assert branch["inversion_max_difference_a3"] < 1e-7
    if args.check:
        check_saved(json.loads(OUTPUT.read_text()), result)
    else:
        OUTPUT.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    print(
        json.dumps(
            {
                identifier: branch["max_abs_fit_minus_full_model_table3_300k_gpa"]
                for identifier, branch in result["branches"].items()
            }
        )
    )


if __name__ == "__main__":
    main()
