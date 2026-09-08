#!/usr/bin/env python3
"""Reproduce the defensible platinum subset of Dorogokupets--Oganov (2007).

The paper's numerical global objective cannot be recreated because its fitted
point inventory and weights are not published.  This audit therefore keeps
three evidence classes separate: exact EOS outputs, recoverable source rows,
and explicitly unavailable observations.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import DorogokupetsOganov2007

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "peritheos" / "data" / "datasets"
DEFAULT_OUTPUT = ROOT / "docs" / "data" / "dorogokupets-oganov-2007-platinum-audit.json"
AVOGADRO_ANGSTROM3_TO_INTERNAL_MOLAR = 0.060_221_407_6


def platinum_eos() -> DorogokupetsOganov2007:
    """Return the rounded Table-I platinum parameterization."""
    return DorogokupetsOganov2007(
        Vinet(0.9091, 276.07, 5.30),
        Tr=298.15,
        theta_B1=95.2,
        d_B1=8.199,
        m_B1=0.329,
        theta_B2=148.4,
        d_B2=4.005,
        m_B2=0.383,
        theta_E1=214.6,
        m_E1=1.211,
        theta_E2=140.8,
        m_E2=1.077,
        gamma0=2.802,
        gamma_inf=1.538,
        beta=5.550,
        anharmonic_a=160.9,
        anharmonic_m=4.06,
        electronic_e=260.0,
        electronic_g=2.4,
        defect_H=32572.0,
        defect_S=0.631,
        n=1.0,
    )


def read_rows(filename: str) -> list[dict[str, str]]:
    with (DATASETS / filename).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def vinet_pressure(volume: np.ndarray, parameters: np.ndarray) -> np.ndarray:
    volume0, bulk_modulus0, derivative = parameters
    x = np.cbrt(volume / volume0)
    return (
        3.0
        * bulk_modulus0
        * (1.0 - x)
        / x**2
        * np.exp(1.5 * (derivative - 1.0) * (1.0 - x))
    )


def final_ruby_pressure(classical_pressure: np.ndarray) -> np.ndarray:
    """Convert the classical Mao ruby coordinate to final Equation (16)."""
    ruby_shift = np.expm1(np.log1p(classical_pressure * 7.665 / 1904.0) / 7.665)
    return 1884.0 * ruby_shift * (1.0 + 5.5 * ruby_shift)


def table_vi_audit(eos: DorogokupetsOganov2007) -> dict[str, Any]:
    temperatures = (298.15, 1000.0, 2000.0, 3000.0)
    rows = {
        1.00: (0.0, 5.309, 12.864, 20.349),
        0.95: (16.207, 21.219, 28.485, 35.822),
        0.90: (38.307, 43.111, 50.199, 57.476),
        0.85: (68.389, 73.069, 80.079, 87.376),
        0.80: (109.388, 114.016, 121.043, 128.434),
        0.75: (165.473, 170.119, 177.253, 184.808),
        0.70: (242.676, 247.403, 254.730, 262.523),
    }
    nodes = []
    residuals = []
    for ratio, published_pressures in rows.items():
        for temperature, published in zip(temperatures, published_pressures):
            calculated = float(eos.pressure(ratio * eos.rt_eos.V0, temperature))
            residual = calculated - published
            residuals.append(residual)
            nodes.append(
                {
                    "volume_ratio": ratio,
                    "temperature_k": temperature,
                    "published_pressure_gpa": published,
                    "calculated_pressure_gpa": calculated,
                    "residual_gpa": residual,
                }
            )
    residuals_array = np.asarray(residuals)
    return {
        "classification": "exact_reconstruction_from_published_rounded_coefficients",
        "source_location": "Table VI, platinum isochors",
        "rows": len(nodes),
        "rmse_gpa": float(np.sqrt(np.mean(residuals_array**2))),
        "max_abs_gpa": float(np.max(np.abs(residuals_array))),
        "nodes": nodes,
    }


def static_diffraction_audit(eos: DorogokupetsOganov2007) -> dict[str, Any]:
    rows = read_rows("platinum-dewaele-2004-table1-compression.csv")
    classical_pressure = np.asarray(
        [float(row["ruby_pressure_classical_gpa"]) for row in rows]
    )
    volume = np.asarray([float(row["atomic_volume_a3"]) for row in rows])
    pressure = final_ruby_pressure(classical_pressure)
    published_curve = np.asarray(
        eos.rt_eos.pressure(volume * AVOGADRO_ANGSTROM3_TO_INTERNAL_MOLAR)
    )
    curve_residuals = published_curve - pressure

    fit = least_squares(
        lambda parameters: vinet_pressure(volume, parameters) - pressure,
        x0=np.asarray((15.099, 270.8, 5.5)),
    )
    fit_residuals = vinet_pressure(volume, fit.x) - pressure
    return {
        "classification": "partial_validation_recoverable_stage_2_rows",
        "dataset_identifier": "platinum_dewaele_2004_table1_compression",
        "source_location": "Dewaele et al. (2004), Table I, all platinum rows",
        "rows": len(rows),
        "pressure_coordinate": (
            "classical Mao pressure inverted to ruby shift, then evaluated with "
            "Dorogokupets--Oganov (2007) Equation (16)"
        ),
        "pressure_range_gpa": [float(np.min(pressure)), float(np.max(pressure))],
        "published_global_curve": {
            "parameters": {
                "V0_a3_per_atom": 15.096,
                "K0_gpa": 276.07,
                "K0_prime": 5.30,
            },
            "rmse_gpa": float(np.sqrt(np.mean(curve_residuals**2))),
            "mean_residual_gpa": float(np.mean(curve_residuals)),
            "max_abs_gpa": float(np.max(np.abs(curve_residuals))),
        },
        "static_only_unweighted_fit": {
            "parameters": {
                "V0_a3_per_atom": float(fit.x[0]),
                "K0_gpa": float(fit.x[1]),
                "K0_prime": float(fit.x[2]),
            },
            "rmse_gpa": float(np.sqrt(np.mean(fit_residuals**2))),
            "max_abs_gpa": float(np.max(np.abs(fit_residuals))),
            "interpretation": (
                "This is a diagnostic one-material pressure-residual fit, not the "
                "paper's coupled multi-observable, six-material objective."
            ),
        },
    }


def cold_energy(eos: DorogokupetsOganov2007, ratio: float) -> float:
    volume0 = eos.rt_eos.V0
    volume = volume0 * ratio
    integral, _ = quad(lambda value: float(eos.rt_eos.pressure(value)), volume0, volume)
    return -integral * 1.0e4


def shock_diagnostic(eos: DorogokupetsOganov2007) -> dict[str, Any]:
    rows = read_rows("platinum-holmes-1989-table3-shock.csv")
    diagnostics = []
    residuals = []
    for row in rows:
        ratio = float(row["sample_density_g_cm3"]) / float(row["density_g_cm3"])
        volume = eos.rt_eos.V0 * ratio
        gamma = eos.gamma_inf + (eos.gamma0 - eos.gamma_inf) * ratio**eos.beta
        reference_pressure = float(eos.rt_eos.pressure(volume))
        energy = cold_energy(eos, ratio)
        calculated = (reference_pressure - gamma / volume * energy / 1.0e4) / (
            1.0 - gamma * (1.0 - ratio) / (2.0 * ratio)
        )
        observed = float(row["pressure_gpa"])
        residual = calculated - observed
        residuals.append(residual)
        diagnostics.append(
            {
                "shot": row["shot"],
                "volume_ratio": ratio,
                "observed_pressure_gpa": observed,
                "equation_15_pressure_gpa": calculated,
                "residual_gpa": residual,
            }
        )
    residuals_array = np.asarray(residuals)
    return {
        "classification": "independent_partial_validation_not_fit_reconstruction",
        "dataset_identifier": "platinum_holmes_1989_table3_shock",
        "source_location": "Holmes et al. (1989), Table III, all seven new shots",
        "rows": len(rows),
        "rmse_gpa": float(np.sqrt(np.mean(residuals_array**2))),
        "mean_residual_gpa": float(np.mean(residuals_array)),
        "max_abs_gpa": float(np.max(np.abs(residuals_array))),
        "diagnostics": diagnostics,
        "qualification": (
            "Dorogokupets--Oganov cite Holmes in the room-temperature comparison "
            "figure, but identify the Shock Wave Database as the stage-1 shock "
            "input. These seven recoverable rows therefore test Equation (15); "
            "they are not asserted to be the unpublished fitted shock point set."
        ),
    }


def build_audit() -> dict[str, Any]:
    eos = platinum_eos()
    entropy = float(eos.thermal_entropy(eos.rt_eos.V0, eos.Tr))
    return {
        "format": "peritheos.dorogokupets-oganov-2007-platinum-audit",
        "format_version": 1,
        "audit_date": "2026-09-08",
        "record_identifier": "platinum_dorogokupets_oganov_2007_vinet_4",
        "source": {
            "doi": "10.1103/PhysRevB.75.024115",
            "author_pdf": "https://uspex-team.org/static/file/PressureScales-PRB-2007.pdf",
        },
        "global_objective": {
            "status": "not_exactly_reconstructible",
            "published_protocol": [
                "Stage 1 weighted least squares to pressure-calibration-independent ambient thermochemical, dilatometric, ultrasonic, and shock-Hugoniot constraints.",
                "Stage 2 addition of static compression data, refitting all coefficients under trial pressure scales.",
                "Final joint treatment of Al, Cu, Ta, W, Au, and Pt while varying the ruby-scale B coefficient.",
            ],
            "residual_families": [
                "zero-pressure heat capacity",
                "relative enthalpy",
                "zero-pressure volume",
                "thermal expansivity",
                "adiabatic bulk modulus",
                "shock Hugoniot pressure from Equation (15)",
                "static X-ray compression",
            ],
            "unavailable_numerical_inputs": [
                "complete material-by-material point inventory",
                "row weights and uncertainty normalization",
                "Shock Wave Database snapshot and selected rows",
                "Collard--McLellan Figure 1 platinum K_S(T) coordinates used in the fit",
                "material-specific thermochemical row mapping",
                "parameter covariance and optimizer details",
            ],
            "conclusion": (
                "The symbolic residual families are recoverable, but no unique "
                "numerical global objective can be reconstructed without inventing "
                "observations or weights. The catalog status remains not_refittable."
            ),
        },
        "exact_reconstructions": {
            "table_vi": table_vi_audit(eos),
            "standard_entropy": {
                "classification": "exact_scalar_reconstruction",
                "calculated_j_mol_k": entropy,
                "paper_calculated_j_mol_k": 41.45,
                "paper_cited_observed_j_mol_k": 41.63,
                "paper_cited_observed_uncertainty_j_mol_k": 0.21,
                "source_location": "text below Table I and reference 97",
            },
        },
        "partial_validations": {
            "static_diffraction": static_diffraction_audit(eos),
            "shock": shock_diagnostic(eos),
        },
        "source_trace": [
            {
                "role": "stage-2 platinum static X-ray compression",
                "reference": "Dewaele, Loubeyre, and Mezouar (2004)",
                "doi": "10.1103/PhysRevB.70.094112",
                "availability": "exact 36-row platinum Table I subset bundled",
            },
            {
                "role": "platinum shock comparison",
                "reference": "Holmes et al. (1989)",
                "doi": "10.1063/1.344177",
                "availability": "exact seven-row Table III new-shot subset bundled",
            },
            {
                "role": "stage-1 platinum adiabatic bulk modulus",
                "reference": "Collard and McLellan (1992)",
                "doi": "10.1016/0956-7151(92)90011-3",
                "availability": (
                    "used Figure 1 coordinates are not tabulated; the article's "
                    "Table 1 polynomial corresponds to the explicitly excluded "
                    "open squares in Dorogokupets--Oganov Figure 7"
                ),
            },
            {
                "role": "stage-1 shock compilation",
                "reference": "Shock Wave Database (2001), paper reference 30",
                "availability": "cited legacy database snapshot and selected rows unavailable",
            },
            {
                "role": "ambient thermochemistry and dilatometry",
                "reference": "paper references 46--52 and 97",
                "availability": (
                    "standard-entropy comparison is printed, but platinum-specific "
                    "row selection and fitting weights are not identified"
                ),
            },
        ],
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="write the canonical JSON")
    parser.add_argument("--check", action="store_true", help="check the canonical JSON")
    args = parser.parse_args()
    rendered = json.dumps(build_audit(), indent=2, sort_keys=True) + "\n"
    if args.check:
        if DEFAULT_OUTPUT.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale audit artifact: {DEFAULT_OUTPUT}")
        return
    if args.write:
        DEFAULT_OUTPUT.write_text(rendered, encoding="utf-8")
        return
    print(rendered, end="")


if __name__ == "__main__":
    main()
