"""Reproduce Dewaele 2008 Table IV fits and check Table VI against Figure 4.

Repository-only source audit. Unweighted pressure least squares is an explicit
diagnostic choice: the publication does not identify weights or the dependent
residual coordinate. Source 95% coefficient widths are not treated as 1 sigma.
"""

from __future__ import annotations

import csv
import itertools
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/dewaele-2008-metals-reproduction.json"
NAMES = ("V0", "K0", "K0_prime")
PARITY_POLICY = (
    "Dewaele 2008 Table IV parity compares free coefficients with printed 95% "
    "parameter widths under the documented fixed constraints. Table VI parity "
    "has source-equation scope only: identified Figure 4 checkpoints pass "
    "declared reading and printed-coefficient bounds. It is not PAW refit or "
    "author-workflow reproduction; the four unresolved Fe LDA sets are excluded."
)
Z = {
    "fe": 2,
    "iron": 2,
    "cobalt_hcp": 2,
    "nickel": 4,
    "zinc_hcp": 2,
    "molybdenum": 2,
    "silver": 4,
    "platinum": 4,
    "gold": 4,
}
LIMITS = {
    "fe": 15,
    "iron": 197,
    "cobalt_hcp": 65.2,
    "nickel": 156.6,
    "zinc_hcp": 156.6,
    "molybdenum": 100,
    "silver": 100,
    "platinum": 94,
    "gold": 140,
}


def vinet(volume, parameters):
    v0, k0, kp = parameters
    x = np.cbrt(np.asarray(volume) / v0)
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def inverse(pressure, parameters):
    return brentq(
        lambda v: float(vinet(v, parameters)) - pressure,
        0.3 * parameters[0],
        parameters[0] * 1.02,
    )


def source():
    return json.loads((DATA / "dewaele-2008-source.json").read_text())


def table4_specs():
    for material, row in source()["tables"]["IV"].items():
        v0, k0, kp, ek0, ekp, constrained = row
        yield (
            material,
            f"{material}_dewaele_2008_table4_vinet",
            [v0, k0, kp],
            [None, ek0, ekp],
            ["V0"],
        )
        for k0, kp, ekp in constrained:
            yield (
                material,
                f"{material}_dewaele_2008_table4_us_k0_{k0:g}_vinet",
                [v0, k0, kp],
                [None, None, ekp],
                ["V0", "K0"],
            )


def table6_specs():
    for row in source()["tables"]["VI"]:
        for state, key in (
            ("static", "static_lattice_tokens"),
            ("298k", "at_298k_tokens"),
        ):
            rid = (
                f"{row['material']}_dewaele_2008_paw_{row['functional']}_{state}_vinet"
            )
            yield (
                row["material"],
                rid,
                list(map(float, row[key])),
                row["functional"],
                state,
            )


def experiment_rows(material):
    name = (
        "dewaele-2008-table3-fe.csv"
        if material == "fe"
        else "dewaele-2008-table2-complete.csv"
    )
    with (DATA / name).open() as stream:
        rows = list(csv.DictReader(stream))
    if material == "fe":
        return [r for r in rows if r["phase"] == "bcc"]
    return [
        r for r in rows if r["material"] == material and r["fit_included"] == "True"
    ]


def table4_outcome(material, rid, parameters, errors, fixed):
    rows = experiment_rows(material)
    pcol = "pressure_revised_gpa" if material == "fe" else "ruby_pressure_dewaele_gpa"
    p = np.array([float(r[pcol]) for r in rows])
    v = np.array([float(r["atomic_volume_a3"]) for r in rows])
    free = [i for i, name in enumerate(NAMES) if name not in fixed]

    def full(x):
        values = np.array(parameters, dtype=float)
        values[free] = x
        return values

    fit = least_squares(
        lambda x: vinet(v, full(x)) - p,
        np.array(parameters)[free],
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    fitted = full(fit.x)
    comparisons = []
    for i in free:
        delta = float(fitted[i] - parameters[i])
        comparisons.append(
            {
                "parameter": NAMES[i],
                "published": parameters[i],
                "published_error": errors[i],
                "refit": float(fitted[i]),
                "refit_error": None,
                "difference": delta,
                "relative_difference": abs(delta / parameters[i]),
                "within_reported_95pct": abs(delta) <= errors[i],
                "within_combined_2sigma": None,
                "similar": abs(delta) <= errors[i],
            }
        )
    within = all(c["within_reported_95pct"] for c in comparisons)
    published_residual = vinet(v, parameters) - p
    return {
        "record_identifier": rid,
        "material_identifier": material,
        "status": "parity" if within else "parity_not_achieved",
        "parity_basis": "within_reported_parameter_errors",
        "source_error_confidence": 0.95,
        "original_author_workflow_parity": False,
        "dataset_identifiers": [f"{material}_dewaele_2008_source_table"],
        "observations": len(rows),
        "selection": "bcc rows only"
        if material == "fe"
        else "all rows except run 2 revised pressure >100 GPa",
        "fit_kind": "isothermal_pv",
        "objective": "unweighted pressure residuals",
        "absolute_sigma": None,
        "weighting": "unweighted diagnostic; source weights not specified",
        "columns": {
            "pressure": pcol,
            "volume": "atomic_volume_a3",
            "temperature": None,
        },
        "observed_pressure_range_gpa": [float(p.min()), float(p.max())],
        "observed_volume_range": [
            float(v.min() * Z[material]),
            float(v.max() * Z[material]),
        ],
        "free_parameters": [NAMES[i] for i in free],
        "fixed_parameters": fixed,
        "fixed_values": {
            name: parameters[NAMES.index(name)] * (Z[material] if name == "V0" else 1)
            for name in fixed
        },
        "parameters": comparisons,
        "rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
        "published_rmse_gpa": float(np.sqrt(np.mean(published_residual**2))),
        "degrees_of_freedom": len(rows) - len(free),
        "solver_success": bool(fit.success),
        "solver_message": str(fit.message),
        "qualification": "Source-fixed V0 and, for the constrained branch, source-adopted ultrasonic K0 are replayed. Every free coefficient is compared directly with its printed 95% width; no refit covariance or source weighting is inferred. This does not reproduce the original low-pressure V0 determination or author workflow.",
    }


def table6_outcome(material, rid, parameters, functional, state):
    with (DATA / "dewaele-2008-figure4-vector.csv").open() as stream:
        rows = [
            r
            for r in csv.DictReader(stream)
            if r["material"] == material
            and r["functional"] == functional
            and (17 if material == "iron" else 5)
            <= float(r["pressure_gpa"])
            <= LIMITS[material]
        ]
    shift = source()["thermal_pressure_298k_gpa"][material] if state == "298k" else 0
    checks = []
    # Printed precision intervals, not statistical parameter uncertainties.
    corners = [
        np.array(parameters) + signs * np.array([0.0005, 0.5, 0.005])
        for signs in map(np.array, itertools.product((-1, 1), repeat=3))
    ]
    for row in rows:
        p = float(row["pressure_gpa"]) + shift
        v = float(row["atomic_volume_a3"])
        pwidth = float(row["pressure_reading_halfwidth_gpa"])
        vwidth = float(row["volume_reading_halfwidth_a3"])
        lo = min(inverse(p + pwidth, q) for q in corners) - vwidth
        hi = max(inverse(max(0, p - pwidth), q) for q in corners) + vwidth
        predicted = inverse(p, parameters)
        checks.append(
            {
                "pressure_gpa": p,
                "source_plot_atomic_volume_a3": v,
                "model_atomic_volume_a3": predicted,
                "volume_difference_a3": predicted - v,
                "source_volume_acceptance_interval_a3": [lo, hi],
                "within_reading_and_print_precision": lo <= v <= hi,
                "source_curve_index": int(row["curve_index"]),
                "source_vertex": int(row["vertex"]),
            }
        )
    passed = bool(checks) and all(
        c["within_reading_and_print_precision"] for c in checks
    )
    if not checks:
        reason = "No separately identified with-semicore alpha-Fe LDA finite-compression curve in Figure 4; raw PAW points are not tabulated. Coefficients are preserved for interchange, pending an independent source benchmark."
    elif not passed:
        reason = f"The hcp-Fe {functional.upper()} Table VI parameterization disagrees with the thick {functional.upper()} Figure 4 path assigned to with-semicore PAW beyond digitization and printed-precision bounds. The follow-up finds agreement with the thin path instead (static maximum difference 0.00244 A3/atom; 298 K 0.00138 A3/atom). The source's table/semicore identity remains unresolved; diagnostic only, not recommended for predictions or pressure calibration. See docs/literature-reproductions/dewaele-2008-metals.md#follow-up-hcp-fe-lda-curve-identity."
    else:
        reason = "Source-equation parity: printed Table VI coefficients pass finite-compression Figure 4 vector checks within the declared reading and printed-coefficient bounds. The 298 K branch is checked against that static curve plus the source's stated constant thermal pressure. Raw PAW fit points, full regression grid and covariance are unavailable; this is graphical equation parity, not parameter-refit or author-workflow reproduction."
    return {
        "record_identifier": rid,
        "material_identifier": material,
        "status": "parity"
        if passed
        else "not_refittable"
        if not checks
        else "parity_not_achieved",
        "comparison_scope": "source_equation",
        **(
            {"parity_basis": "source_curve_reading_and_print_precision"}
            if passed
            else {}
        ),
        "original_author_workflow_parity": False,
        "dataset_identifiers": [f"{material}_dewaele_2008_figure4"] if rows else [],
        "observations": len(checks),
        "reason": reason,
        "fit_kind": "source_curve_validation_only",
        "free_parameters": [],
        "fixed_parameters": [],
        "curve_validation": {
            "passed": passed,
            "checks": checks,
            "source_thermal_pressure_shift_gpa": shift,
            "method": "1 GPa pressure reading bound, 0.04 A3/atom volume reading bound, and independent half-last-digit coefficient bounds; no statistical confidence inferred",
        },
    }


@lru_cache(maxsize=1)
def reproduce():
    records = {
        rid: table4_outcome(m, rid, p, e, fixed)
        for m, rid, p, e, fixed in table4_specs()
    }
    records.update(
        {
            rid: table6_outcome(m, rid, p, f, state)
            for m, rid, p, f, state in table6_specs()
        }
    )
    return {
        "format": "peritheos.dewaele-2008-reproduction",
        "records": records,
        "source_pdf_sha256": source()["source_pdf_sha256"],
    }


def ledger_outcome(record):
    rid = record["identifier"]
    return dict(reproduce()["records"][rid])


if __name__ == "__main__":
    report = reproduce()
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    from collections import Counter

    print(dict(Counter(r["status"] for r in report["records"].values())))
    for rid, row in report["records"].items():
        if row["status"] == "parity_not_achieved":
            print(rid, row.get("reason", "coefficient interval not met"))
