#!/usr/bin/env python3
"""Reconstruct the staged Anderson et al. (1989) gold EOS derivation.

This deliberately does not fit Equation (29) to Table V.  Table V is output
from that equation, while Tables I-IV combine heterogeneous literature inputs
and derived thermodynamic quantities.  The reproducible claims are the source's
separate one-dimensional regressions, its four thermodynamic-consistency
trials, and evaluation of the published equation against the published grid.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"
DEFAULT_REPORT = (
    ROOT / "docs" / "data" / "anderson-1989-gold-thermal-eos-reproduction.json"
)
RECORD_IDENTIFIER = "gold_anderson_1989_bm3_1"
SOURCE_PDF_SHA256 = "38a4a6276270c3777d20c0437e22103612c32f5c0b622ebba32f3532499f9024"
DATASETS = {
    "table1": "gold-anderson-1989-table1-thermodynamic-inputs.csv",
    "table2": "gold-anderson-1989-table2-temperature-derivatives.csv",
    "table3": "gold-anderson-1989-table3-bulk-moduli.csv",
    "table4": "gold-anderson-1989-table4-anharmonic-diagnostics.csv",
    "table5": "gold-anderson-1989-table5-pressure-grid.csv",
}
PUBLISHED_DKDT = {
    "6.39": -0.0052,
    "6.12": -0.0071,
    "5.5": -0.0115,
    "5.21": -0.0135,
}
CONSISTENCY_TRIALS = (
    # K0', Eq. (18) left side, Eq. (18) right side, unrounded dKT/dT|V
    (6.39, 0.66, 0.72, -0.00515),
    (6.12, 0.93, 0.99, None),
    (5.50, 1.55, 1.60, None),
    (5.21, 1.84, 1.89, None),
)


def _rows(name: str) -> list[dict[str, str]]:
    with (DATA / DATASETS[name]).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _linear_fit(
    rows: list[dict[str, str]], x_name: str, y_name: str, *, minimum_x: float
) -> tuple[float, float, int]:
    selected = [row for row in rows if row[y_name] and float(row[x_name]) >= minimum_x]
    x = np.asarray([float(row[x_name]) for row in selected])
    y = np.asarray([float(row[y_name]) for row in selected])
    slope, intercept = np.polyfit(x, y, 1)
    return float(slope), float(intercept), len(selected)


def reproduce() -> dict[str, Any]:
    table1 = _rows("table1")
    table2 = _rows("table2")
    table3 = _rows("table3")
    table4 = _rows("table4")
    table5 = _rows("table5")

    low_300 = next(
        row
        for row in table1
        if row["temperature_k"] == "300" and row["elastic_series"] == "neighbours_alers"
    )
    table3_300 = next(
        row
        for row in table3
        if row["temperature_k"] == "300" and row["elastic_series"] == "neighbours_alers"
    )
    table4_300 = next(row for row in table4 if row["temperature_k"] == "300")

    temperature = 300.0
    molar_mass_g_mol = 196.967
    density_g_cm3 = 19.30
    alpha = float(low_300["volume_expansivity_1e6_k_inverse"]) * 1e-6
    ks_gpa = float(table3_300["ks_ambient_pressure_gpa"])
    gamma = float(table4_300["gamma"])
    cp_printed_scaled = float(low_300["cp_printed_0p1_j_g_k"])
    cp_reconstructed_scaled = float(low_300["cp_reconstructed_0p1_j_g_k"])
    cp_equation4_scaled = (
        (5.66 + 1.222e-3 * temperature + 3.4e3 * temperature**-2)
        * 4.184
        / molar_mass_g_mol
        * 10.0
    )
    cp_from_gamma_scaled = alpha * ks_gpa * 1000.0 / (density_g_cm3 * gamma) * 10.0

    def gamma_from_cp(cp_scaled: float) -> float:
        return alpha * ks_gpa * 1000.0 / (density_g_cm3 * cp_scaled * 0.1)

    table2_dkt = next(row for row in table2 if row["property"] == "dKt_dT_at_P")
    ambient_low = [
        row
        for row in table3
        if row["elastic_series"] == "neighbours_alers"
        and float(row["temperature_k"]) >= 150.0
    ]
    ambient_high = [row for row in table3 if row["elastic_series"] == "chang_himmel"]

    def ambient_kt_fit(rows: list[dict[str, str]]) -> tuple[float, float, int]:
        x = np.asarray([float(row["temperature_k"]) for row in rows])
        y = np.asarray([float(row["kt_ambient_pressure_gpa"]) for row in rows])
        slope, intercept = np.polyfit(x, y, 1)
        return float(slope), float(intercept), len(rows)

    low_slope, low_intercept, low_observations = ambient_kt_fit(ambient_low)
    high_slope, high_intercept, high_observations = ambient_kt_fit(ambient_high)

    trial_columns = {
        "6.39": "kt_constant_v_kprime_6p39_gpa",
        "6.12": "kt_constant_v_kprime_6p12_gpa",
        "5.5": "kt_constant_v_kprime_5p5_gpa",
        "5.21": "kt_constant_v_kprime_5p21_gpa",
    }
    constant_volume_fits: dict[str, dict[str, float | int]] = {}
    for label, column in trial_columns.items():
        slope, intercept, observations = _linear_fit(
            table3, "temperature_k", column, minimum_x=200.0
        )
        constant_volume_fits[label] = {
            "observations": observations,
            "slope_gpa_k": slope,
            "intercept_gpa": intercept,
            "published_slope_gpa_k": PUBLISHED_DKDT[label],
        }

    thermal_slope, thermal_intercept, thermal_observations = _linear_fit(
        table4,
        "temperature_k",
        "thermal_pressure_gpa",
        minimum_x=225.0,
    )

    consistency = []
    for kprime, left, right, unrounded_slope in CONSISTENCY_TRIALS:
        consistency.append(
            {
                "k0_prime_trial": kprime,
                "equation18_left": left,
                "equation18_right": right,
                "relative_mismatch": abs(left - right) / right,
                "unrounded_dkt_dt_v_gpa_k_if_printed": unrounded_slope,
            }
        )

    document = get_material_document("gold")
    record = Material.from_eosmat(
        document, record_identifiers=[RECORD_IDENTIFIER]
    ).eos_records[0]
    observed = np.asarray([float(row["pressure_this_study_gpa"]) for row in table5])
    volumes = np.asarray(
        [record.reference_volume * (1.0 - float(row["compression"])) for row in table5]
    )
    temperatures = np.asarray([float(row["temperature_k"]) for row in table5])
    calculated = np.asarray(record.pressure(volumes, temperatures), dtype=float)
    residuals = calculated - observed

    avogadro = 6.02214076e23
    volume_a3 = molar_mass_g_mol / density_g_cm3 * 4.0e24 / avogadro
    resources = {}
    for name, filename in DATASETS.items():
        path = DATA / filename
        resources[name] = {
            "path": str(path.relative_to(ROOT)),
            "rows": len(_rows(name)),
            "sha256": _sha256(path),
        }

    return {
        "format": "peritheos.anderson-1989-gold-thermal-eos-reproduction",
        "format_version": 1,
        "record_identifier": RECORD_IDENTIFIER,
        "primary_source": {
            "doi": "10.1063/1.342969",
            "source_pdf_sha256": SOURCE_PDF_SHA256,
        },
        "conclusion": {
            "global_fit_status": "not_defined_by_source",
            "classification": "not_refittable",
            "reason": (
                "The source performs a staged thermodynamic synthesis with "
                "heterogeneous literature properties, graphical smoothing, "
                "separate linear fits, numerical integration, and qualitative "
                "thermodynamic-consistency trials. It publishes neither a common "
                "observation matrix nor a global objective, weights, exact "
                "integration protocol, parameter covariance, or fitted residuals."
            ),
        },
        "resources": resources,
        "ambient_volume_reconstruction": {
            "atomic_mass_g_mol": molar_mass_g_mol,
            "density_g_cm3": density_g_cm3,
            "formula_units_per_fcc_cell": 4,
            "calculated_volume_a3": volume_a3,
            "stored_volume_a3": record.reference_volume,
        },
        "table1_300k_cp_audit": {
            "printed_0p1_j_g_k": cp_printed_scaled,
            "reconstructed_0p1_j_g_k": cp_reconstructed_scaled,
            "equation4_0p1_j_g_k": cp_equation4_scaled,
            "equation12_inverted_0p1_j_g_k": cp_from_gamma_scaled,
            "gamma_from_printed_cp": gamma_from_cp(cp_printed_scaled),
            "gamma_from_reconstructed_cp": gamma_from_cp(cp_reconstructed_scaled),
            "table4_gamma": gamma,
            "finding": (
                "The literal Table I value is internally inconsistent; about 1.288 "
                "is required by both Equation (4) and Table IV."
            ),
        },
        "table2_ambient_kt_temperature_derivatives": {
            "low_temperature": {
                "selection": "Neighbours-Alers Table III rows from 150 through 300 K",
                "observations": low_observations,
                "slope_gpa_k": low_slope,
                "intercept_gpa": low_intercept,
                "published_slope_gpa_k": float(
                    table2_dkt["low_temperature_coefficient_1e3"]
                )
                * 1e-3,
            },
            "high_temperature": {
                "selection": "all Chang-Himmel Table III rows from 300 through 550 K",
                "observations": high_observations,
                "slope_gpa_k": high_slope,
                "intercept_gpa": high_intercept,
                "published_slope_gpa_k": float(
                    table2_dkt["high_temperature_coefficient_1e3"]
                )
                * 1e-3,
            },
        },
        "table3_constant_volume_linear_fits": {
            "selection": "all printed nonblank rows from 200 through 550 K",
            "trials": constant_volume_fits,
        },
        "table4_thermal_pressure_linear_fit": {
            "selection": (
                "all printed rows at or above 225 K, the first tabulated point "
                "above 1.3*theta_D for theta_D=170 K"
            ),
            "observations": thermal_observations,
            "slope_gpa_k": thermal_slope,
            "intercept_gpa": thermal_intercept,
            "published_slope_gpa_k": 0.00714,
            "published_intercept_gpa": -0.40,
        },
        "equation18_thermodynamic_consistency_trials": {
            "selection_result": "K0' near 5.2-5.5 favored; no scalar optimization reported",
            "trials": consistency,
            "final_equation29_k0_prime": 5.4823,
            "final_value_provenance": (
                "adopted from the revised Heinz-Jeanloz 300 K equation of state, "
                "not fitted by Anderson et al."
            ),
        },
        "table5_equation29_reproduction": {
            "states": len(table5),
            "rmse_gpa": float(np.sqrt(np.mean(residuals**2))),
            "max_abs_residual_gpa": float(np.max(np.abs(residuals))),
            "rounding_tolerance_gpa": 0.005,
        },
        "nonreproducible_source_steps": [
            "graphical interpolation and hand smoothing of expansivity below 250 K",
            "exact numerical-integration implementation for Equations (8), (13), and (19)",
            "a global objective function and cross-property weighting",
            "a joint parameter covariance matrix",
            "row-level upstream observations not reprinted in Table I",
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--write", action="store_true", help="write the JSON report")
    parser.add_argument(
        "--check", action="store_true", help="check the committed report"
    )
    args = parser.parse_args()
    report = reproduce()
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.write:
        DEFAULT_REPORT.write_text(rendered, encoding="utf-8")
    elif args.check:
        if DEFAULT_REPORT.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale report: run {Path(__file__).name} --write")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
