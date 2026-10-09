#!/usr/bin/env python3
"""Audit Morard (2026) FeS against unchanged Table S1 observations.

The quadrature below is independent of Peritheos's EOS implementations.
Conditional gamma-only fits are diagnostics, not replacement catalog records.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro, R
from scipy.optimize import least_squares
from scipy.special import roots_legendre

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/fes-morard-2026-table-s1.csv"
OUTPUT = ROOT / "docs/data/morard-2026-fes-audit.json"
RECORD_ID = "fes_vi_morard_2026_bm3_mgd"
DATASET_ID = "fes_morard_2026_table_s1_pvt"
CELL_PER_MOLAR = 4e24 / Avogadro
NODES, WEIGHTS = roots_legendre(48)
NODES = (NODES + 1) / 2
WEIGHTS = WEIGHTS / 2
QUALIFICATION = (
    "The 146 supplied Table S1 rows are available, but the complete Sata (2008) "
    "and Ohfuji (2007) cold-fit inputs, final row selection, EosFit file, weights "
    "and explicit Debye-temperature law are unavailable. Gamma-only fits hold "
    "the printed cold coefficients fixed and are conditional diagnostics. "
    "Neither thermal convention reproduces the stated +/-3 GPa bound. "
    "The source record is deferred; no diagnostic replaces its coefficients."
)


def observations():
    """Return every workbook row, in original order and units."""
    with DATA.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def pressure(volume, temperature, gamma0=2.42, law="integrated_gruneisen"):
    """BM3 plus Debye energy difference; V in cm3/mol of FeS, T in K."""
    volume, temperature = np.broadcast_arrays(volume, temperature)
    x = volume / 15.40
    gamma = gamma0 * x
    theta = (
        417 * np.exp(-gamma0 * (x - 1))
        if law == "integrated_gruneisen"
        else 417 * x ** (-gamma)
    )

    def energy(t):
        z = theta[..., None] / t[..., None] * NODES
        d3 = 3 * np.sum(WEIGHTS * NODES**2 * z / np.expm1(z), axis=-1)
        return 6 * R * t * d3

    cold = (
        1.5
        * 115.5
        * (x ** (-7 / 3) - x ** (-5 / 3))
        * (1 + 0.75 * (4.99 - 4) * (x ** (-2 / 3) - 1))
    )
    return cold + gamma * (energy(temperature) - energy(np.full_like(x, 300))) / (
        volume * 1000
    )


def residual_summary(residual, rows):
    """Retain the location of discrepancies, not only their average."""
    return {
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "mean_gpa": float(np.mean(residual)),
        "min_gpa": float(np.min(residual)),
        "max_gpa": float(np.max(residual)),
        "max_abs_gpa": float(np.max(np.abs(residual))),
        "rows_exceeding_3_gpa": [int(x) for x in rows[np.abs(residual) > 3]],
    }


def calibration_replay(d):
    """Match GSAS FeS lattice triples to Table S1 and independently replay KCl."""
    raw = DATA.parent / "morard_2026_sources/FeS-GSAS.csv"
    with raw.open(newline="", encoding="utf-8-sig") as handle:
        rows = list(csv.reader(handle, delimiter=";"))
    blocks = []
    current = None
    for line, row in enumerate(rows, 1):
        if row and row[0].startswith("#"):
            current = {"source_scan": row[0], "source_line": line}
            blocks.append(current)
        elif current is not None and row[0] == "FeSVI":
            current["fes_lattice"] = [float(x) for x in row[3:6]]
            current["temperature_k"] = float(row[13])
        elif current is not None and row[0] == "KCl B2":
            current["marker_a_angstrom"] = float(row[3])
            current["marker_a_error_angstrom"] = float(rows[line][3])
    results = []
    for i in range(len(d["source_row"])):
        lattice = [d[f"{axis}_angstrom"][i] for axis in "abc"]
        matches = [
            b
            for b in blocks
            if b.get("fes_lattice") == lattice
            and b.get("temperature_k") == d["temperature_k"][i]
        ]
        if (
            not matches
            or len(
                {
                    (b.get("marker_a_angstrom"), b.get("marker_a_error_angstrom"))
                    for b in matches
                }
            )
            != 1
        ):
            raise ValueError(f"Table S1 row {i + 3}: {len(matches)} GSAS matches")
        source = matches[0]
        volume = source["marker_a_angstrom"] ** 3
        # Dewaele (2012), Tables III/V and Equation (2), cell volume in A3.
        x = (volume / 54.5) ** (1 / 3)
        cold = 3 * 17.2 * (1 - x) / x**2 * np.exp(1.5 * (5.89 - 1) * (1 - x))
        temperature = 0.75 * d["temperature_k"][i] + 0.25 * 295
        p = cold + 0.00224 * (temperature - 300)
        results.append(
            {
                **source,
                "matching_source_lines": [b["source_line"] for b in matches],
                "source_row": int(d["source_row"][i]),
                "marker_volume_angstrom3": volume,
                "derived_marker_temperature_k": float(temperature),
                "replayed_pressure_gpa": float(p),
                "pressure_difference_gpa": float(p - d["pressure_gpa"][i]),
            }
        )
    return {
        "matched_rows": len(results),
        "max_abs_pressure_difference_gpa": max(
            abs(r["pressure_difference_gpa"]) for r in results
        ),
        "reference_eos_record": "kcl_b2_dewaele_2012_vinet_3",
        "marker_temperature_law": "0.75*T_FeS + 0.25*295 K",
        "qualification": "This is an exact numerical reconstruction with the public GSAS markers. The final fit file is unavailable. The main article describes an average marker temperature but does not print these weights. Earlier PT-workbook marker refinements are not substituted.",
        "rows": results,
    }


def reproduce():
    """Return the source replay, source arithmetic checks and partial refits."""
    d = observations()
    v = d["volume_angstrom3"] / CELL_PER_MOLAR
    t, p = d["temperature_k"], d["pressure_gpa"]
    a, b, c = (d[f"{axis}_angstrom"] for axis in "abc")
    errors = np.array([d[f"sigma_{axis}_angstrom"] for axis in "abc"])
    relative = errors / np.array([a, b, c])
    replays, fits = {}, {}
    for law in ("integrated_gruneisen", "variable_exponent"):
        residual = pressure(v, t, law=law) - p
        replays[law] = residual_summary(residual, d["source_row"])
        replays[law]["residuals_gpa"] = residual.tolist()
        replays[law]["over_3_gpa_at_pressure_ge_40"] = int(
            np.sum((np.abs(residual) > 3) & (p >= 40))
        )
        for objective in (
            "unweighted_pressure",
            "sigma_pressure",
            "effective_variance",
        ):

            def residuals(value):
                gamma = value[0]
                scale = np.ones_like(p)
                if objective != "unweighted_pressure":
                    scale = d["sigma_pressure_gpa"].copy()
                if objective == "effective_variance":
                    dv = v * 1e-5
                    dpdv = (
                        pressure(v + dv, t, gamma, law)
                        - pressure(v - dv, t, gamma, law)
                    ) / (2 * dv)
                    dt = 0.01
                    dpdt = (
                        pressure(v, t + dt, gamma, law)
                        - pressure(v, t - dt, gamma, law)
                    ) / (2 * dt)
                    scale = np.sqrt(
                        scale**2
                        + (dpdv * d["sigma_volume_angstrom3"] / CELL_PER_MOLAR) ** 2
                        + (dpdt * d["sigma_temperature_k"]) ** 2
                    )
                return (pressure(v, t, gamma, law) - p) / scale

            result = least_squares(
                residuals, [2.42], bounds=(0.1, 8), xtol=1e-12, ftol=1e-12, gtol=1e-12
            )
            variance = float(np.sum(result.fun**2) / (len(p) - 1))
            standard_error = float(np.sqrt(variance / np.sum(result.jac**2)))
            fitted = float(result.x[0])
            fits[f"{law}:{objective}"] = {
                "gamma0": fitted,
                "conditional_standard_error": standard_error,
                "conditional_variance": standard_error**2,
                "error_convention": "Local Jacobian covariance scaled by residual variance; printed cold coefficients fixed. Source confidence and cross-row covariance unknown.",
                "converged": bool(result.success),
                "observations": len(p),
                "residuals": residual_summary(
                    pressure(v, t, fitted, law) - p, d["source_row"]
                ),
            }
    checkpoints = [
        {
            "volume_cm3_mol": volume,
            "temperature_k": temp,
            "pressure_gpa": float(pressure(np.array(volume), np.array(temp))),
        }
        for volume, temp in [(15.4, 300), (12, 300), (12, 1500), (10.1, 2500)]
    ]
    return {
        "format": "peritheos.morard-2026-fes-audit",
        "audit_date": "2026-10-02",
        "record_identifier": RECORD_ID,
        "dataset_identifier": DATASET_ID,
        "dataset_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
        "observations": len(p),
        "calibration_replay": calibration_replay(d),
        "source_rows": [int(x) for x in d["source_row"]],
        "volume0_angstrom3": 15.4 * CELL_PER_MOLAR,
        "volume0_error_angstrom3": 0.45 * CELL_PER_MOLAR,
        "table_arithmetic": {
            "max_abs_volume_product_difference_angstrom3": float(
                np.max(abs(a * b * c - d["volume_angstrom3"]))
            ),
            "max_abs_linear_volume_error_difference_angstrom3": float(
                np.max(
                    abs(
                        d["volume_angstrom3"] * relative.sum(axis=0)
                        - d["sigma_volume_angstrom3"]
                    )
                )
            ),
            "max_abs_quadrature_volume_error_difference_angstrom3": float(
                np.max(
                    abs(
                        d["volume_angstrom3"] * np.sqrt((relative**2).sum(axis=0))
                        - d["sigma_volume_angstrom3"]
                    )
                )
            ),
            "max_abs_paper_C_definition_difference": float(
                np.max(abs(np.sqrt(3) * b / c - d["order_parameter_c"]))
            ),
            "max_abs_reciprocal_C_definition_difference": float(
                np.max(abs(c / (np.sqrt(3) * b) - d["order_parameter_c"]))
            ),
        },
        "ranges": {
            key: [float(d[key].min()), float(d[key].max())]
            for key in ["pressure_gpa", "temperature_k", "volume_angstrom3"]
        },
        "published_coefficient_replays": replays,
        "conditional_gamma_only_fits": fits,
        "independent_checkpoints": checkpoints,
        "qualification": QUALIFICATION,
    }


def ledger_outcome(record):
    """Classify staged coefficient reproduction separately from residual claims."""
    report = reproduce()
    path = ROOT / "docs/data/fes-author-staged/manifest.json"
    staged = json.loads(path.read_text(encoding="utf-8"))
    fit = staged["cases"]["reproduced-errors-tr300"]
    cold = fit["cold_stage"]["refined_parameters_final_cycle"]
    thermal = fit["thermal_stage"]["refined_parameters_final_cycle"]
    comparisons = []
    for name, published, error, fitted, factor in (
        ("eos.V0", 15.4, 0.45, cold["V0"], CELL_PER_MOLAR),
        ("eos.K0", 115.5, 27.91, cold["K0"], 1),
        ("eos.K0_prime", 4.99, 0.51, cold["Kp"], 1),
        ("thermal.gamma0", 2.42, 0.03, thermal["Gamm0"], 1),
    ):
        value, esd = fitted["value"], fitted["esd"]
        relative = abs(value - published) / abs(published)
        limit = (
            0.05
            if name.endswith("V0")
            else 0.15
            if name.endswith("K0")
            else 0.20
            if name.endswith("K0_prime")
            else 0.25
        )
        comparisons.append(
            {
                "parameter": name,
                "published": published * factor,
                "published_error": error * factor,
                "refit": value * factor,
                "refit_error": esd * factor,
                "relative_difference": relative,
                "similar": relative <= limit,
                "within_combined_2sigma": None,
                "reported_intervals_overlap": abs(value - published) <= error + esd,
                "uncertainty_qualification": "Source confidence and cross-fit covariance are unknown; thermal esd is conditional on the fixed cold curve. No formal independent combined-sigma parity is asserted.",
            }
        )
    assert all(p["similar"] and p["reported_intervals_overlap"] for p in comparisons)
    return {
        "status": "similar",
        "reproduction_status": "similar",
        "residual_claim_reproduction_status": "not_reproduced",
        "dataset_identifiers": ["fes_morard_2026_author_without_additional_cold_pvt"],
        "observations": 167,
        "stage_observations": {"cold": 21, "thermal": 146, "quenched": 0},
        "parameters": comparisons,
        "free_parameters": ["eos.V0", "eos.K0", "eos.K0_prime", "thermal.gamma0"],
        "fixed_parameters": ["theta0", "q", "n", "Tr=300 K"],
        "fit_kind": "cold_then_fixed_cold_thermal",
        "fitting_software": "Official EosFit7c 7.60; independent Peritheos pressure replay",
        "objective": "EosFit pressure residuals with supplied P/V/T uncertainty weights",
        "published_rmse_gpa": staged["cases"]["published-errors-tr300"][
            "thermal_stage"
        ]["initial_source_replay"]["rmse_gpa"],
        "rmse_gpa": fit["thermal_residual_summary"]["rmse_gpa"],
        "rmse_scope": "146 thermal observations; published replay uses printed cold coefficients and gamma0 at Tr=300 K, refit uses reproduced weighted cold coefficients",
        "reduced_chi_square": fit["thermal_stage"]["eosfit_weighted_chi2_per_dof"],
        "qualification": "Classified similar: staged reproduced-cold coefficients are compatible within reported uncertainty intervals, and the equal-weight control puts all four central estimates inside published error widths. Exact GUI output and the complete +/-3 GPa residual envelope are separate unresolved findings. Published coefficients are preserved; this is not an automatic replacement EOS or a formal independent uncertainty-parity claim.",
        "reason": record.get("scientific_validation", {}).get("note", QUALIFICATION),
        "author_input_followup": record.get("author_input_followup", {}),
        "staged_manifest": path.relative_to(ROOT).as_posix(),
        "staged_manifest_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        "partial_validation": report["conditional_gamma_only_fits"],
        "published_curve_validation": report["published_coefficient_replays"],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    encoded = json.dumps(report, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.check:
        if json.loads(OUTPUT.read_text(encoding="utf-8")) != report:
            raise SystemExit("Morard FeS audit is stale")
    else:
        OUTPUT.write_text(encoded, encoding="utf-8")
    print(
        json.dumps(
            {
                "observations": report["observations"],
                "replays": {
                    law: {k: v for k, v in values.items() if k != "residuals_gpa"}
                    for law, values in report["published_coefficient_replays"].items()
                },
            }
        )
    )


if __name__ == "__main__":
    main()
