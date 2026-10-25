"""Reconstruct conditional uncertainties for the bounded 2026-10-25 ledger audit.

Repository maintenance only. Source coefficients and observations are never
replaced. Every calculation reuses the recorded selection, objective and fixed
parameters. Unknown source confidence and calibration covariance stay unknown.
"""

from __future__ import annotations

import argparse
import copy
import csv
import hashlib
import json
import math
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.stats import t as student_t

from peritheos import convert_eosmat, load_eosmat
from scripts.reconstruction_plot_data import executable, observations

ROOT = Path(__file__).resolve().parents[1]
INPUT = ROOT / "docs/data/uncertainty-reconstruction-inputs.json"
OUTPUT = ROOT / "docs/data/uncertainty-reconstruction.json"
REPORT = "docs/literature-reproductions/uncertainty-reconstruction.md"
DATA = ROOT / "peritheos/data/datasets"
PARITY_POLICY = "The bounded 2026-10-25 uncertainty follow-up accepts conditional numerical parameter parity within printed source error widths, with explicit reconstructed or recovered standard errors. Published 95% widths remain 95% widths; unknown source confidence, original weighting, calibration covariance and author-workflow scope remain qualified. At parameter_agreement scope, digitization/run-selection hypotheses and diagnostic fixed masks are explicit conditions of the numerical comparison; unresolved source-data, calibration and benchmark findings remain separate. Parameter parity does not validate those hypotheses, recover author workflow or change use recommendations."
COMBINED_ERROR_CASES = {
    "ca_perovskite_tetragonal_chen_2018_vinet": 2.0,
    "superhydrous_phase_b_lt_litasov_2007_htbm_t03": 1.0,
}
COMBINED_ERROR_POLICY = "Chen (2018) and Litasov (2007) HTBM/T03 instead use the ledger's combined-error and numerical-similarity criterion: Chen's printed two-sigma source widths are divided by two; Litasov's unspecified-confidence widths remain diagnostic comparison quantities. Neither assessment establishes a formal confidence level or author-workflow reproduction."
PARAMETER_AGREEMENT_QUALIFICATIONS = {
    "diamond_datchi_2007_vinet_1": "Agreement is conditional on the source-unconfirmed run-2/equal-volume hypothesis: 14 of 24 digitized markers, fixed source V0 and H2005 recalibration. Selection/weighting was supported by coefficient agreement, so this is hypothesis-dependent parameter agreement, not independent validation of that selection. Residual-scaled errors exclude graphical and calibration systematics; raw numerical observations and author workflow remain unavailable.",
    "platinum_fei_2004_bm3_300k": "Agreement concerns K0 and K0-prime only, conditional on diagnostic V0=60.38 A^3 fixed and equal pressure-residual weights for seven Table 2 300 K rows. The source's fixed/free mask is unrecovered; V0 is not independently recovered. This preliminary narrow-range fit has a poorly constrained derivative and is not the preferred Pt pressure scale. No thermal or author-workflow reproduction is claimed.",
    "siderite_fe095mn005_litasov_2013_bm3": "Agreement concerns K0 and K0-prime on all 27 printed 300 K rows, conditional on source-fixed V0=293.4 A^3, original printed Au pressures and equal pressure-residual weights. The published-curve pressure RMS remains 0.22958 GPa, outside the manuscript's reported 0.11-0.17 GPa range; that source-data agreement finding is unresolved separately. Exact pressure-scale realization, original objective/weights and author workflow remain unrecovered. No thermal-coefficient conflict or use restriction is resolved by this parameter comparison.",
}


def table(filename):
    with (DATA / filename).open(newline="") as stream:
        return list(csv.DictReader(stream))


def col(rows, name):
    return np.array([float(r[name] or 0) for r in rows])


def semantic_digest(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def verify_preserved_sources():
    inputs = json.loads(INPUT.read_text())
    for material, digest in inputs["source_dataset_declarations_sha256"].items():
        document = load_eosmat(
            ROOT / "peritheos/data/materials" / (material + ".eosmat")
        )
        if semantic_digest(document.get("datasets", [])) != digest:
            raise ValueError("Source dataset declarations changed: " + material)
        for record in document["eos_records"]:
            expected = inputs["source_record_sha256"].get(record["identifier"])
            if expected is not None:
                source = copy.deepcopy(record)
                source["evidence"]["reproduction"].pop("assessment", None)
                if semantic_digest(source) != expected:
                    raise ValueError(
                        "Source record changed outside its reconstruction assessment: "
                        + record["identifier"]
                    )


def jacobian(function, values, step=1e-4):
    """Five-point central derivative in dimensionless parameter coordinates."""
    values = np.asarray(values, dtype=float)
    scale = np.maximum(np.abs(values), 1e-20)
    columns = []
    for index in range(len(values)):
        delta = np.zeros(len(values))
        delta[index] = step * scale[index]
        columns.append(
            (
                -function(values + 2 * delta)
                + 8 * function(values + delta)
                - 8 * function(values - delta)
                + function(values - 2 * delta)
            )
            / (12 * step)
        )
    return np.column_stack(columns), scale


def covariance_from_residual(function, values, names):
    """Residual-scaled local covariance; no source uncertainty is manufactured."""
    residual = np.asarray(function(np.asarray(values)), dtype=float)
    j, scale = jacobian(function, values)
    check, _ = jacobian(function, values, step=5e-5)
    derivative_change = float(np.linalg.norm(check - j) / np.linalg.norm(j))
    if derivative_change > 2e-4:
        raise ValueError(f"Jacobian fails step-halving check: {derivative_change}")
    if not np.all(np.isfinite(residual)) or not np.all(np.isfinite(j)):
        raise ValueError("Nonfinite residual or Jacobian")
    _, singular, vt = np.linalg.svd(j, full_matrices=False)
    if singular[-1] / singular[0] <= 1e-10:
        raise ValueError(
            "Parameters are not identifiable: scaled Jacobian rank deficient"
        )
    dof = len(residual) - len(values)
    if dof <= 0:
        raise ValueError("No residual degrees of freedom")
    variance = float(residual @ residual / dof)
    score = float(
        np.max(
            np.abs(j.T @ residual)
            / (np.linalg.norm(j, axis=0) * np.linalg.norm(residual))
        )
    )
    if score > 0.01:
        raise ValueError(
            f"Recorded coefficients are not stationary for the replayed objective: {score}"
        )
    matrix = ((vt.T / singular**2) @ vt) * variance * np.outer(scale, scale)
    result = describe_covariance(matrix, values, names, dof)
    result.update(
        residual_variance=variance,
        max_abs_residual_jacobian_cosine=score,
        residual_sum_squares=float(residual @ residual),
        residual_count=len(residual),
        scaled_jacobian_condition_number=float(singular[0] / singular[-1]),
        derivative_step_halving_relative_change=derivative_change,
        method="Residual-scaled local least-squares covariance: SSE/(n-k) times inverse J^T J; five-point scaled central Jacobian, checked at half step.",
    )
    return result


def describe_covariance(matrix, values, names, dof):
    matrix = np.asarray(matrix)
    errors = np.sqrt(np.diag(matrix))
    if not np.all(np.isfinite(errors)) or np.any(errors <= 0):
        raise ValueError("Unusable conditional standard errors")
    multiplier = student_t.ppf(0.975, dof)
    return dict(
        parameters=names,
        matrix=matrix.tolist(),
        standard_errors=dict(zip(names, map(float, errors))),
        correlation_matrix=(matrix / np.outer(errors, errors)).tolist(),
        degrees_of_freedom=float(dof),
        local_student_t_95_intervals={
            n: [float(v - multiplier * s), float(v + multiplier * s)]
            for n, v, s in zip(names, values, errors)
        },
    )


def staged_komabayashi(row):
    """Propagate the cold-stage estimator through thermal fitting on shared rows."""
    from scripts.reproduce_wang_zhou_komabayashi_mgsio3 import (
        observations as source_rows,
    )
    from scripts.reproduce_wang_zhou_komabayashi_mgsio3 import thermal_model

    data = np.asarray(source_rows()[row["record_identifier"]])
    v, p, t = data[:, 1:].T
    values = np.array([r["refit"] for r in row["parameters"]])

    def model(x):
        return thermal_model(*x).pressure(v, t)

    scaled, scale = jacobian(model, values)
    j = scaled
    cold = t == 300
    j[cold, 1:] = 0  # The constrained thermal increment is exactly zero at Tr.
    a, b = j[:, :1], j[:, 1:]
    first = np.zeros((1, len(p)))
    first[:, cold] = np.linalg.solve(a[cold].T @ a[cold], a[cold].T)
    second = np.linalg.solve(b.T @ b, b.T) @ (np.eye(len(p)) - a @ first)
    influence = np.vstack([first, second])
    residual_map = np.eye(len(p)) - j @ influence
    effective_dof = float(np.sum(residual_map**2))
    residual = model(values) - p
    variance = float(residual @ residual / effective_dof)
    covariance = variance * (influence @ influence.T) * np.outer(scale, scale)
    result = describe_covariance(
        covariance, values, [r["parameter"] for r in row["parameters"]], effective_dof
    )
    result.update(
        method="Linear influence propagation through the cold and thermal stages on shared observations; common independent pressure-residual variance estimated with tr((I-JA)(I-JA)^T) effective residual degrees of freedom.",
        residual_count=len(p),
        residual_sum_squares=float(residual @ residual),
        residual_variance=variance,
        estimator_influence=(scale[:, None] * influence).tolist(),
        prediction_jacobian=(j / scale).tolist(),
        stage_sizes=[int(cold.sum()), len(p)],
        stage_uncertainty_propagated=True,
        input_array_sha256=hashlib.sha256(data.tobytes()).hexdigest(),
    )
    return result


def finger_rows(row, record):
    from peritheos.eos.finger1981 import CELL_TO_MOLAR, Finger1981Argon

    rows = table("argon-fcc-finger-1981-table1.csv")
    v = col(rows, "molar_volume_cm3_mol") / CELL_TO_MOLAR
    p = col(rows, "pressure_gpa")
    fit = row["reproduction"]["refits"]["equal_pressure"]
    names = ["K0_prime", "K0_double_prime"]
    values = np.array([fit["K0_prime"], fit["K0_double_prime_gpa_inverse"]])
    comparisons = [
        dict(
            parameter=n,
            published=record["eos"]["parameters"][n],
            published_error=record["parameter_errors"][n],
            refit=float(x),
            difference=float(x - record["eos"]["parameters"][n]),
            relative_difference=float(
                abs(x - record["eos"]["parameters"][n])
                / abs(record["eos"]["parameters"][n])
            ),
            within_combined_2sigma=None,
            similar=True,
        )
        for n, x in zip(names, values)
    ]
    return comparisons, lambda x: Finger1981Argon(*x).pressure(v, 293) - p


def ordinary_residual(document, record, row):
    """Reuse the same scientific row adapters as the audited curve replay."""
    series = observations(document, record, row)
    if len(series.pressure) != row["observations"]:
        raise ValueError("Source selection differs from original reconstruction")
    params = row["parameters"]

    def model(x, volume=None, temperature=None):
        changed = [{**p, "refit": float(value)} for p, value in zip(params, x)]
        eos = executable(document, record, changed, "refit")
        return np.asarray(
            eos.pressure(
                series.volume if volume is None else volume,
                series.temperature if temperature is None else temperature,
            ),
            dtype=float,
        )

    values = np.array([p["refit"] for p in params])
    # Verify we reconstructed the original curve/selection, before estimating errors.
    pressure_residual = model(values) - series.pressure
    expected = row.get("rmse_gpa")
    if expected is not None and not np.isclose(
        np.sqrt(np.mean(pressure_residual**2)), expected, rtol=5e-5, atol=2e-6
    ):
        raise ValueError("Replayed RMS differs from original reconstruction")
    weights = np.ones(len(series.pressure))
    rid = row["record_identifier"]
    weight_note = "Equal weights on the original pressure residuals; observed volumes and fixed coefficients conditioned on."
    if rid == "diamond_datchi_2007_vinet_1":

        def residual(x):
            changed = [{**p, "refit": float(value)} for p, value in zip(params, x)]
            eos = executable(document, record, changed, "refit")
            return (np.asarray(eos.volume(series.pressure)) - series.volume) / 8

        return (
            residual,
            "Original run-2/equal-volume diagnostic in atomic A^3; digitization and calibration systematics are not estimated by residual scatter.",
        )
    if rid == "argon_fcc_ono_2020_vinet":
        rows = table("argon-fcc-ono-2020-table1.csv")
        from scripts.reproduce_argon import vinet

        v = col(rows, "volume_a3")
        derivative = (
            vinet(v + 1e-4, 184.5, 1.07, 8.02) - vinet(v - 1e-4, 184.5, 1.07, 8.02)
        ) / 2e-4
        weights = np.hypot(
            col(rows, "pressure_sigma_gpa"), derivative * col(rows, "volume_sigma_a3")
        )
        weight_note = "Original effective coordinate widths, frozen at published coefficients; residual-scaled covariance conditional on those fixed weights."
    elif "_chidester_2018_" in rid:
        material = "thorianite" if rid.startswith("thorianite") else "tho2_cotunnite"
        rows = [
            r
            for r in table(f"chidester_2018_{material}_pvt.csv")
            if r["included_in_fit"] == "1"
        ]
        v = series.volume
        source = record["eos"]["parameters"]
        x = (source["V0"] / v) ** (1 / 3)
        derivative = -0.5 * source["K0"] * (7 * x**7 - 5 * x**5) / v
        c = record["thermal"]["parameters"]["alpha_KT"]
        factor = document["formula_units_per_cell"] / 0.602214076
        weights = np.sqrt(
            col(rows, "pressure_error_gpa") ** 2
            + (derivative * col(rows, "molar_volume_error_cm3_mol") * factor) ** 2
            + (c * col(rows, "temperature_error_k")) ** 2
            + np.where(series.temperature == 300, (0.03 * series.pressure) ** 2, 0)
        )
        weight_note = "Original effective P/V/T weights frozen at published coefficients, including the audit's Ar 3% cold-pressure sensitivity; diagonal and residual-scaled, with shared scale covariance unresolved."
    elif "_farla_2025_" in rid:
        from scripts.reproduce_rubidium_halides_farla_2025 import (
            SOURCE,
            read_rows,
            reference_pressure,
        )

        compound = rid.split("_", 1)[0]
        rows = [
            r
            for r in read_rows(compound, "observations")
            if r["volume_a3"]
            and (r["source_table"] == "S1" or float(r["temperature_k"]) > 300)
        ]
        v, t = series.volume, series.temperature
        molar, theta, _ = SOURCE[compound]

        def evaluate(v, t):
            return reference_pressure(v, t, molar, theta, values)

        dv = (evaluate(v + 1e-4, t) - evaluate(v - 1e-4, t)) / 2e-4
        dt = (evaluate(v, t + 0.01) - evaluate(v, t - 0.01)) / 0.02
        weights = np.sqrt(
            col(rows, "pressure_uncertainty_gpa") ** 2
            + (dv * col(rows, "volume_uncertainty_a3")) ** 2
            + (dt * col(rows, "temperature_uncertainty_k")) ** 2
        )
        weight_note = "Final frozen effective P/V/T weights at the existing converged fit, without a pressure floor. Conditional on the analyst row mask and source-reported CsCl pressures; calibration replay and benchmark discrepancies remain unresolved."
    if np.any(weights <= 0) or not np.all(np.isfinite(weights)):
        raise ValueError("Missing or zero effective residual width")
    return lambda x: (model(x) - series.pressure) / weights, weight_note


def reconstruct(row):
    rid = row["record_identifier"]
    document = convert_eosmat(
        load_eosmat(
            ROOT / "peritheos/data/materials" / (row["material_identifier"] + ".eosmat")
        ),
        version=3,
    )
    record = next(r for r in document["eos_records"] if r["identifier"] == rid)
    params = copy.deepcopy(row.get("parameters", []))
    note = "Existing reconstructed errors retained; no source confidence convention inferred."
    if rid == "argon_fcc_finger_1981_murnaghan2_debye":
        params, residual = finger_rows(row, record)
        uncertainty = covariance_from_residual(
            residual, [p["refit"] for p in params], [p["parameter"] for p in params]
        )
        note = "All 19 molar-volume rows and the full original Murnaghan/zero-point/Debye model at 293 K; other coefficients fixed."
    elif rid == "mgsio3_post_perovskite_komabayashi_2008_thermal_bm3":
        uncertainty = staged_komabayashi(row)
        note = "Cold-stage K0 uncertainty propagated into the thermal stage, including covariance from reuse of the three cold observations; fixed V0, K0-prime and alpha(300) conditioned on."
    elif rid == "kcl_b2_ma_2024_bm3_mgd":
        fitted = json.loads(
            (ROOT / "docs/data/ma-2024-kcl-reproduction.json").read_text()
        )["fits"]["reported_errors"]
        transform = np.diag([1 / 0.602214076, 1, 1])
        covariance = transform @ np.array(fitted["diagnostic_covariance"]) @ transform
        uncertainty = describe_covariance(
            covariance,
            [p["refit"] for p in params],
            [p["parameter"] for p in params],
            16,
        )
        uncertainty["method"] = (
            "Recovered existing joint acoustic/Walker residual-scaled diagnostic covariance; molar V0 converted to Z=1 conventional-cell units on both covariance axes."
        )
        note = "Joint 11 KT(V) plus eight Walker P(V) observations; same deposited relative weights as the stored diagnostic, not author covariance. Static coefficients only; thermal and calibration qualifications remain."
    elif all(p.get("refit_error") is not None for p in params):
        uncertainty = dict(
            parameters=[p["parameter"] for p in params],
            standard_errors={p["parameter"]: p["refit_error"] for p in params},
            method="Previously reconstructed standard errors retained unchanged; no missing off-diagonal covariance inferred.",
        )
    else:
        residual, note = ordinary_residual(document, record, row)
        uncertainty = covariance_from_residual(
            residual, [p["refit"] for p in params], [p["parameter"] for p in params]
        )
    for p in params:
        p["refit_error"] = uncertainty["standard_errors"][p["parameter"]]
    within = all(
        p.get("published_error") is not None
        and p["published_error"] > 0
        and abs(p["refit"] - p["published"]) <= p["published_error"]
        for p in params
    )
    decision = None
    if not within:
        decision = "At least one central value lies outside its printed source parameter width; retained similarity under this follow-up's conservative printed-width criterion, even where the broader stored error comparison passes."
    combined_parity = rid in COMBINED_ERROR_CASES and all(
        p.get("within_combined_2sigma") is True
        and p.get("similar") is True
        and abs(p["refit"] - p["published"])
        <= 2
        * math.hypot(p["published_error"] / COMBINED_ERROR_CASES[rid], p["refit_error"])
        for p in params
    )
    if combined_parity:
        decision = None
    status = "similar" if decision else "parity"
    decision = (
        decision
        or (
            "Every compared free coefficient passes the combined-error comparison and numerical similarity screen. "
            + (
                "The source explicitly reports two-sigma parameter widths; dividing those widths by two before combining with reconstruction standard errors also passes. "
                if COMBINED_ERROR_CASES[rid] == 2.0
                else "Source confidence is unspecified, so the combined-two-width comparison is diagnostic, not a formal two-sigma confidence statement. "
            )
            + "The printed-source-width-only cutoff is not required for this parameter-agreement scope. Original weights, calibration covariance and author workflow remain unrecovered; shared-data cross-fit covariance is unknown."
            if combined_parity
            else None
        )
        or "Every compared free central value lies inside its printed source error width, with finite reconstructed standard errors under the explicit conditional error model. This is numerical parameter parity, not a formal source-confidence test or recovery of author workflow."
    )
    if rid == "rbcl_b2_farla_2025_bm3_mgd_cscl" and status == "parity":
        decision += (
            " Parameter agreement is conditional on fixed V0/theta0, the documented analyst row mask, source-reported CsCl pressures and frozen diagonal effective P/V/T weights."
            " CsCl calibration replay and source benchmark discrepancies remain unresolved as separate findings; neither calibration parity, publication-digit benchmark parity nor author-workflow reproduction is claimed."
        )
    if rid in PARAMETER_AGREEMENT_QUALIFICATIONS:
        decision += " " + PARAMETER_AGREEMENT_QUALIFICATIONS[rid]
    uncertainty["previous_assessment"] = {
        k: row[k] for k in ("status", "qualification", "reason") if k in row
    }
    uncertainty["assumptions"] = [
        note,
        "Local standard errors condition on the recorded row selection, fixed/adopted parameters and pressure-coordinate reduction; calibration systematics and cross-row correlations are not propagated.",
        "Residual-scaled estimates assume independent residuals with common variance after any frozen weighting. They are not published parameter uncertainties or a recovered author error budget.",
        "Student-t intervals, where supplied, are marginal local approximations, not exact nonlinear or joint confidence regions.",
        "Source error confidence is retained as reported (including Dewaele's 95% intervals); missing confidence remains unknown and no new combined-two-sigma flag is inferred.",
    ]
    return dict(
        record_identifier=rid,
        material_identifier=row["material_identifier"],
        parameters=params,
        covariance=uncertainty,
        weighting=note,
        status=status,
        decision=decision,
        original_source_notes=row.get("source_notes", ""),
    )


def source_bindings(inputs):
    """Bind this follow-up to factual source resources, not generated metadata."""
    resources = {}
    for row in inputs:
        doc = load_eosmat(
            ROOT / "peritheos/data/materials" / (row["material_identifier"] + ".eosmat")
        )
        selected = set(row.get("dataset_identifiers", []))
        for dataset in doc.get("datasets", []):
            if dataset["identifier"] not in selected:
                continue
            resource = dataset.get("resource", {})
            if "path" in resource:
                path = ROOT / "peritheos/data" / resource["path"]
                resources[str(path.relative_to(ROOT))] = hashlib.sha256(
                    path.read_bytes()
                ).hexdigest()
    resources[str(INPUT.relative_to(ROOT))] = hashlib.sha256(
        INPUT.read_bytes()
    ).hexdigest()
    ma = ROOT / "docs/data/ma-2024-kcl-reproduction.json"
    resources[str(ma.relative_to(ROOT))] = hashlib.sha256(ma.read_bytes()).hexdigest()
    resources["peritheos/data/pressure-calibrations.json"] = hashlib.sha256(
        (ROOT / "peritheos/data/pressure-calibrations.json").read_bytes()
    ).hexdigest()
    return dict(sorted(resources.items()))


@lru_cache(maxsize=1)
def saved_followup():
    verify_preserved_sources()
    report = json.loads(OUTPUT.read_text())
    for name, digest in report["source_sha256"].items():
        if hashlib.sha256((ROOT / name).read_bytes()).hexdigest() != digest:
            raise ValueError(
                f"Uncertainty follow-up source changed: {name}; rerun its audit"
            )
    return {r["record_identifier"]: r for r in report["records"]}


def augment_outcome(record, outcome):
    """Apply only the named, independently calculated follow-up assessments."""
    result = saved_followup().get(record["identifier"])
    if result is None:
        return outcome
    if "error" in result:
        raise ValueError(
            "Unresolved uncertainty calculation: " + result["record_identifier"]
        )
    if outcome.get("status") not in {"similar", "parity"}:
        raise ValueError(
            "Original reconstruction failed or changed scope: " + record["identifier"]
        )
    if (
        not outcome.get("parameters")
        and record["identifier"] != "argon_fcc_finger_1981_murnaghan2_debye"
    ):
        raise ValueError(
            "Original coefficient comparison is missing: " + record["identifier"]
        )
    original = {p["parameter"]: p for p in outcome.get("parameters", [])}
    comparisons = []
    for p in result["parameters"]:
        old = original.get(p["parameter"])
        if old and (
            old["published"] != p["published"]
            or old.get("published_error") != p.get("published_error")
            or not math.isclose(old["refit"], p["refit"], rel_tol=2e-5, abs_tol=1e-25)
        ):
            raise ValueError(
                "Changed source or reconstructed coefficients: " + record["identifier"]
            )
        comparisons.append(
            {
                **(old or {}),
                **p,
                "within_reported_error": bool(
                    abs(p["refit"] - p["published"]) <= p["published_error"]
                ),
            }
        )
    updated = copy.deepcopy(outcome)
    # Keep old qualification wording explicitly historical, including all source limits.
    baseline = next(
        r
        for r in json.loads(INPUT.read_text())["records"]
        if r["record_identifier"] == record["identifier"]
    )
    history = " ".join(
        str(baseline[k]) for k in ("qualification", "reason") if baseline.get(k)
    )
    updated.update(
        status=result["status"],
        parameters=comparisons,
        covariance=result["covariance"],
        weighting=result["weighting"],
        qualification=result["decision"]
        + " "
        + result["weighting"]
        + " Earlier audit statement (superseded only as to missing reconstruction uncertainty and the resulting numerical label): "
        + history
        + " See docs/literature-reproductions/uncertainty-reconstruction.md for the complete conditional uncertainty budget.",
    )
    if "reason" in updated:
        updated["reason"] = (
            "Earlier audit statement, before the conditional uncertainty follow-up: "
            + baseline.get("reason", "")
        )
    if result["status"] == "parity":
        updated["parity_basis"] = (
            "conditional_combined_error_parameter_agreement"
            if record["identifier"] in COMBINED_ERROR_CASES
            else "within_reported_parameter_errors"
        )
    if record["identifier"] == "platinum_fei_2004_bm3_300k":
        # These are reconstruction constraints, not a recovered author mask.
        updated.update(
            free_parameters=["K0", "K0_prime"],
            validation_fixed_parameters=["V0"],
            fixed_values={"V0": 60.38},
            objective="unweighted_pressure_residuals",
            selection="All seven Table 2 300 K rows; diagnostic V0 fixed, not a recovered source fixed/free protocol",
        )
    uncertainty = result["covariance"]
    if "matrix" in uncertainty:
        updated["parameter_covariance"] = uncertainty["matrix"]
        updated["covariance_convention"] = uncertainty["method"]
        updated["absolute_sigma"] = False
    if "degrees_of_freedom" in uncertainty:
        updated["degrees_of_freedom"] = uncertainty["degrees_of_freedom"]
    return updated


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--record")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    verify_preserved_sources()
    inputs = json.loads(INPUT.read_text())["records"]
    if args.record and args.record not in {r["record_identifier"] for r in inputs}:
        parser.error("Unknown record in this bounded follow-up")
    results = []
    for row in inputs:
        if args.record and row["record_identifier"] != args.record:
            continue
        try:
            result = reconstruct(row)
            results.append(result)
            print(
                row["record_identifier"],
                "OK",
                result["covariance"]["standard_errors"],
                flush=True,
            )
        except (ValueError, KeyError, RuntimeError) as error:
            print(row["record_identifier"], "FAILED", str(error), flush=True)
            results.append(
                dict(record_identifier=row["record_identifier"], error=str(error))
            )
    if any("error" in r for r in results):
        raise SystemExit(
            "Uncertainty reconstruction failed; source assessments were not updated"
        )
    report = dict(
        format="peritheos.conditional-uncertainty-followup",
        format_version=1,
        generated_with="scripts/reconstruct_similar_eos_uncertainties.py",
        scope="parameter_agreement",
        source_sha256=source_bindings(inputs),
        records=results,
    )
    if args.check:
        saved = json.loads(OUTPUT.read_text())
        if report["source_sha256"] != saved["source_sha256"]:
            raise SystemExit("Source provenance changed")
        stored = {r["record_identifier"]: r for r in saved["records"]}
        for r in results:
            old = stored[r["record_identifier"]]
            for key in (
                "decision",
                "weighting",
                "material_identifier",
                "original_source_notes",
            ):
                if r[key] != old[key]:
                    raise SystemExit("Changed scientific qualification: " + key)
            for key in (
                "method",
                "parameters",
                "assumptions",
                "stage_sizes",
                "stage_uncertainty_propagated",
            ):
                if r["covariance"].get(key) != old["covariance"].get(key):
                    raise SystemExit("Changed uncertainty convention: " + key)
            if r["status"] != old["status"] or r["parameters"] != old["parameters"]:
                # Refit errors are numerical results, while source/central values are literal.
                for p, q in zip(r["parameters"], old["parameters"]):
                    for key in set(p) | set(q):
                        if key == "refit_error":
                            if not math.isclose(
                                p[key], q[key], rel_tol=3e-4, abs_tol=1e-30
                            ):
                                raise SystemExit("Changed reconstructed error")
                        elif p.get(key) != q.get(key):
                            raise SystemExit("Changed source or comparison metadata")
                if r["status"] != old["status"]:
                    raise SystemExit("Changed assessment")
            if "matrix" in r["covariance"] and not np.allclose(
                r["covariance"]["matrix"],
                old["covariance"]["matrix"],
                rtol=3e-4,
                atol=1e-30,
            ):
                raise SystemExit("Changed covariance")
    elif args.record:
        raise SystemExit(
            "Selective calculations are read-only; use --record with --check"
        )
    else:
        OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")


if __name__ == "__main__":
    main()
