"""Audit Yokoo (2009) Pt branches and recovered electronic pressure nodes."""

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

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/yokoo-2009-platinum-reproduction.json"
RECORD = "platinum_yokoo_2009_vinet_300k"
V0 = 60.55  # density-derived fcc cell volume, four atoms
TABLE = "platinum-yokoo-2009-table5-isochores.csv"
SOURCE = "platinum-yokoo-2009-source.json"
ELECTRONIC = "tsuchiya-kawamura-2002-table1-electronic-pressure.csv"


def source_vinet(volume):
    """Table IV 300 K Vinet fit, independent of the library evaluator."""
    x = (np.asarray(volume, dtype=float) / V0) ** (1 / 3)
    return 3 * 276.4 * (1 - x) / x**2 * np.exp(1.5 * (5.48 - 1) * (1 - x))


def gamma_theta(ratio):
    """Equations 10 and 12, normalized at ambient V0, not cold Vc."""
    ratio = np.asarray(ratio, dtype=float)
    gamma = 2.63 * (1 + 0.39 * (ratio**5.2 - 1))
    theta = 230 * ratio ** (-(1 - 0.39) * 2.63) * np.exp(-(gamma - 2.63) / 5.2)
    return gamma, theta


def phonon_pressure(ratio, temperature):
    """Absolute pressure, Eq. 4/5/8; 0 K energy excludes zero point."""
    gamma, theta = gamma_theta(ratio)
    if temperature == 0:
        return 0.0
    x = float(theta) / temperature
    integral = quad(lambda z: z**3 / np.expm1(z), 0, x, epsabs=1e-11)[0]
    energy = 9 * R * temperature * integral / x**3
    molar_volume_si = V0 * ratio * N_A * 1e-30 / 4
    return float(gamma * energy / molar_volume_si / 1e9)


def bm3(ratio, cold_volume_ratio):
    """Eq. 2, Vc distinct from the 300 K normalization V0."""
    z = np.asarray(ratio) / cold_volume_ratio
    return (
        1.5
        * 288.4
        * (z ** (-7 / 3) - z ** (-5 / 3))
        * (1 + 0.75 * (5.05 - 4) * (z ** (-2 / 3) - 1))
    )


def table_rows():
    with (DATA / TABLE).open() as stream:
        return list(csv.DictReader(stream))


@lru_cache(maxsize=1)
def reproduce():
    from peritheos import get_eos_record

    rows = table_rows()
    cold = [r for r in rows if float(r["temperature_k"]) == 0]
    r = np.array([float(row["volume_ratio"]) for row in cold])
    p0 = np.array([float(row["pressure_gpa"]) for row in cold])
    p300 = np.array(
        [
            float(row["pressure_gpa"])
            for row in rows
            if float(row["temperature_k"]) == 300
        ]
    )
    # Effective Vc inferred ONLY from model-output checkpoints. Not an author refit.
    fit = least_squares(
        lambda vc: bm3(r, vc[0]) - p0, [0.994], xtol=1e-13, ftol=1e-13, gtol=1e-13
    )
    grid = V0 * np.array([1.0, 0.98, 0.9, 0.8, 0.7, 0.6])
    expected = source_vinet(grid)
    native = get_eos_record(RECORD).pressure(grid)
    thermal = []
    with (DATA / ELECTRONIC).open() as stream:
        electronic_rows = list(csv.DictReader(stream))
    electronic_by_temperature = {
        float(r["temperature_k"]): float(r["platinum_electronic_pressure_gpa"])
        for r in electronic_rows
    }
    for row in rows:
        ratio = float(row["volume_ratio"])
        t = float(row["temperature_k"])
        base = p0[np.argmin(abs(r - ratio))]
        phonon = phonon_pressure(ratio, t)
        thermal.append(
            {
                **row,
                "phonon_pressure_gpa": phonon,
                "table_total_minus_0k_minus_phonon_gpa": float(row["pressure_gpa"])
                - base
                - phonon,
                "source_electronic_pressure_gpa": electronic_by_temperature[t],
                "printed_components_minus_table_gpa": base
                + phonon
                + electronic_by_temperature[t]
                - float(row["pressure_gpa"]),
            }
        )
    # Holmes rows remain shock states; independently check their momentum reduction.
    with (DATA / "platinum-holmes-1989-table3-shock.csv").open() as stream:
        shock = list(csv.DictReader(stream))
    shock_checks = [
        {
            "shot": row["shot"],
            "source_pressure_gpa": float(row["pressure_gpa"]),
            "momentum_pressure_gpa": float(row["sample_density_g_cm3"])
            * float(row["shock_velocity_km_s"])
            * float(row["mass_velocity_km_s"]),
            "yokoo_linear_Us_residual_km_s": float(row["shock_velocity_km_s"])
            - (3.635 + 1.543 * float(row["mass_velocity_km_s"])),
        }
        for row in shock
    ]
    for row, check in zip(shock, shock_checks):
        rho = float(row["sample_density_g_cm3"])
        us = float(row["shock_velocity_km_s"])
        up = float(row["mass_velocity_km_s"])
        # Half-last-digit intervals, separate from the experimental 2-sigma bounds.
        lower = (rho - 0.005) * (us - 0.0005) * (up - 0.0005)
        upper = (rho + 0.005) * (us + 0.0005) * (up + 0.0005)
        check["momentum_last_digit_interval_gpa"] = [lower, upper]
        check["rounding_intervals_overlap"] = (
            lower <= check["source_pressure_gpa"] + 0.05
            and upper >= check["source_pressure_gpa"] - 0.05
        )
    paired = []
    with (DATA / "platinum-dewaele-2004-table1-compression.csv").open() as stream:
        static_pt = list(csv.DictReader(stream))
    with (DATA / "gold-dewaele-2004-table1-compression.csv").open() as stream:
        static_au = list(csv.DictReader(stream))
    au_by_pressure = {
        (row["ruby_pressure_classical_gpa"], row["ruby_pressure_revised_gpa"]): (
            index,
            row,
        )
        for index, row in enumerate(static_au, 1)
    }
    assert len(au_by_pressure) == len(static_au) == 37
    assert len(static_pt) == 36
    for index, pt in enumerate(static_pt, 1):
        au_index, au = au_by_pressure[
            (pt["ruby_pressure_classical_gpa"], pt["ruby_pressure_revised_gpa"])
        ]
        assert pt["ruby_pressure_classical_gpa"] == au["ruby_pressure_classical_gpa"]
        assert pt["ruby_pressure_revised_gpa"] == au["ruby_pressure_revised_gpa"]
        paired.append(
            {
                "platinum_csv_row_one_based": index,
                "gold_csv_row_one_based": au_index,
                "platinum_atomic_volume_a3": float(pt["atomic_volume_a3"]),
                "platinum_atomic_volume_error_a3": float(
                    pt["atomic_volume_uncertainty_a3"]
                ),
                "gold_atomic_volume_a3": float(au["atomic_volume_a3"]),
                "gold_atomic_volume_error_a3": float(
                    au["atomic_volume_uncertainty_a3"]
                ),
                "original_classical_ruby_pressure_gpa": float(
                    pt["ruby_pressure_classical_gpa"]
                ),
                "original_revised_ruby_pressure_gpa": float(
                    pt["ruby_pressure_revised_gpa"]
                ),
                "derived_yokoo_vinet_pt_pressure_gpa": float(
                    source_vinet(4 * float(pt["atomic_volume_a3"]))
                ),
                "role": "comparison_only_not_yokoo_fit_target",
            }
        )
    files = [
        TABLE,
        SOURCE,
        ELECTRONIC,
        "platinum-holmes-1989-table3-shock.csv",
        "platinum-dewaele-2004-table1-compression.csv",
        "gold-dewaele-2004-table1-compression.csv",
    ]
    return {
        "doi": "10.1103/PhysRevB.80.104114",
        "record_identifier": RECORD,
        "author_fit_status": "not_refittable",
        "qualification": "Independent implementation of the published 300 K Vinet fit. Table V is derived full-model output, not Vinet output or observations. Cold-volume diagnostic is circular checkpoint reconstruction; no author-fit parity. Electronic pressure nodes and the Au-only erratum are recovered. Exact author phonon/electronic normalization and original joint-fit inputs remain unresolved; a separate explicitly derived PVT reconstruction is registered.",
        "input_sha256": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in files
        },
        "density_derived_V0_cell_a3": 4 * 195.084 / 21.4 / N_A * 1e24,
        "catalog_V0_cell_a3": V0,
        "equation_grid": {
            "volume_a3": grid.tolist(),
            "pressure_gpa": expected.tolist(),
        },
        "independent_equation_max_difference_gpa": float(
            np.max(abs(native - expected))
        ),
        "vinet_minus_full_model_300k_table_gpa": (source_vinet(r * V0) - p300).tolist(),
        "max_abs_vinet_minus_full_model_300k_table_gpa": float(
            np.max(abs(source_vinet(r * V0) - p300))
        ),
        "cold_checkpoint_diagnostic": {
            "Vc_over_V0": float(fit.x[0]),
            "K0_fixed_gpa": 288.4,
            "K0_prime_fixed": 5.05,
            "rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
            "max_abs_residual_gpa": float(np.max(abs(fit.fun))),
            "kind": "derived_output_only_not_independent_fit",
        },
        "gamma_theta_checkpoints": [
            {
                "volume_ratio": float(v),
                "gamma": float(gamma_theta(v)[0]),
                "theta_k": float(gamma_theta(v)[1]),
            }
            for v in [1, 0.8, 0.6]
        ],
        "table5_thermal_term_diagnostic": thermal,
        "recovered_electronic_pressure_check": {
            "rows": len(electronic_rows),
            "erratum_changes_platinum": False,
            "reference": "Absolute relative to 0 K; no 300 K subtraction",
            "printed_components_max_abs_difference_gpa": max(
                abs(r["printed_components_minus_table_gpa"]) for r in thermal
            ),
            "pressure_reconstruction_report": "docs/data/yokoo-2009-platinum-thermal-reconstruction.json",
        },
        "upstream_shock_momentum_checks": shock_checks,
        "dewaele_simultaneous_volume_pairs": paired,
        "max_abs_shock_momentum_rounding_difference_gpa": max(
            abs(s["source_pressure_gpa"] - s["momentum_pressure_gpa"])
            for s in shock_checks
        ),
    }


def ledger_outcome(record):
    return {
        "record_identifier": record["identifier"],
        "status": "not_refittable",
        "fit_kind": "published_reference_isotherm_from_joint_shock_thermodynamic_model",
        "observations": 0,
        "dataset_identifiers": [
            "platinum_yokoo_2009_table5_isochores",
            "platinum_holmes_1989_table3_shock",
        ],
        "reason": "Table IV reconstructs the 300 K Vinet fit, but Table V is full-model output and the seven Holmes rows are shock states. Electronic-pressure nodes are recovered and support the separately registered derived PVT reconstruction. Original Pt ambient thermodynamic inputs, complete shock selection, electronic energy, exact author normalization and row-level objective/weights remain unavailable. These cannot be fitted as static 300 K observations. See [Yokoo audit](literature-reproductions/yokoo-2009-platinum.md).",
    }


def check_saved(saved, current):
    """Check every archived result, allowing only floating-point roundoff."""
    if isinstance(current, dict):
        assert saved.keys() == current.keys()
        for key in current:
            check_saved(saved[key], current[key])
    elif isinstance(current, list):
        assert len(saved) == len(current)
        for left, right in zip(saved, current):
            check_saved(left, right)
    elif isinstance(current, float):
        np.testing.assert_allclose(saved, current, atol=1e-9, rtol=1e-10)
    else:
        assert saved == current


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    assert result["independent_equation_max_difference_gpa"] < 1e-9
    assert len(result["table5_thermal_term_diagnostic"]) == 168
    assert len(result["dewaele_simultaneous_volume_pairs"]) == 36
    assert all(
        row["rounding_intervals_overlap"]
        for row in result["upstream_shock_momentum_checks"]
    )
    if args.check:
        saved = json.loads(OUTPUT.read_text(encoding="utf-8"))
        check_saved(saved, result)
    else:
        OUTPUT.write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "author_fit_status",
                    "independent_equation_max_difference_gpa",
                    "max_abs_vinet_minus_full_model_300k_table_gpa",
                ]
            }
        )
    )


if __name__ == "__main__":
    main()
