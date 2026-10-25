"""Review the remaining 96 similar findings without changing scientific sources.

Repository maintenance only. Reuse existing estimates and fit protocols; recover
identified source errors and calculate the two missing Farla covariances. A
source-width or printed-precision match is not a statistical equivalence test.
Run --integrate after reviewing the generated report; --check is read-only.
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import math
from collections import Counter
from functools import lru_cache

from peritheos import convert_eosmat, load_eosmat, validate_eosmat_document
from peritheos._eosmat_format import format_eosmat
from scripts.reconstruct_similar_eos_uncertainties import (
    ROOT,
    covariance_from_residual,
    ordinary_residual,
    semantic_digest,
)

INPUT = ROOT / "docs/data/similar-parameter-review-inputs.json"
OUTPUT = ROOT / "docs/data/similar-parameter-review.json"
REPORT = ROOT / "docs/literature-reproductions/similar-parameter-review.md"
POLICY = "The remaining-similar parameter-scope review separates numerical parameter agreement, published decimal precision and registered-refit replay from calibration, benchmark, data-agreement, author-workflow and recommendation findings. Published source widths are never inferred from reconstruction errors. A source-width-only central-value comparison may retain unknown reconstruction covariance. Explicit fixed constraints are excluded from free-coefficient recovery. Missing uncertainty is not zero; no favorable covariance rescaling is selected to obtain parity."
TANGE_SOURCE = "https://agupubs.onlinelibrary.wiley.com/doi/full/10.1029/2011JB008988"
TANGE_ERRORS = {
    "bridgmanite_tange_2012_bm3_1": {"theta0": 150.0, "gamma0": 0.09, "q": 0.3},
    "bridgmanite_tange_2012_vinet_2": {"theta0": 140.0, "gamma0": 0.09, "q": 0.3},
}
REFITS = {
    "cu_sun_2010_low_vn_refit",
    "argon_fcc_chen_2010_bm3_digitized_refit",
}


def verify_sources():
    inputs = json.loads(INPUT.read_text())
    for material, expected in inputs["source_document_sha256"].items():
        document = copy.deepcopy(
            load_eosmat(ROOT / "peritheos/data/materials" / (material + ".eosmat"))
        )
        for record in document["eos_records"]:
            record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)
        if semantic_digest(document) != expected:
            raise ValueError("Scientific source document changed: " + material)
    for path, expected in inputs["source_resources_sha256"].items():
        if hashlib.sha256((ROOT / path).read_bytes()).hexdigest() != expected:
            raise ValueError("Source resource changed: " + path)


def numerical_screen(p):
    name = p["parameter"].split(".")[-1]
    difference = abs(p["refit"] - p["published"])
    relative = difference / abs(p["published"]) if p["published"] else math.inf
    if name == "K0_prime":
        return difference <= 1 or relative <= 0.2
    limit = (
        0.05
        if name == "V0"
        else 0.15
        if name == "K0"
        else 0.25
        if name in {"theta0", "gamma0", "q"}
        else 0.30
        if name.startswith("alpha") or name.startswith("dK_dT")
        else 0.20
    )
    return relative <= limit


def normalized_comparisons(row):
    result = copy.deepcopy(row.get("parameters", []))
    for p in result:
        p["published_error"] = p.get("published_error", p.get("published_sigma"))
        p["refit_error"] = p.get("refit_error", p.get("refit_sigma"))
        p["difference"] = p["refit"] - p["published"]
    return result


def review(row):
    rid = row["record_identifier"]
    native = load_eosmat(
        ROOT / "peritheos/data/materials" / (row["material_identifier"] + ".eosmat")
    )
    ps = normalized_comparisons(row)
    notes = []
    covariance = None
    fixed = []
    precision = {}
    scope = "parameter_agreement"
    if rid in TANGE_ERRORS:
        for p in ps:
            if p["parameter"] in TANGE_ERRORS[rid]:
                if p["published_error"] is not None:
                    raise ValueError("Tange source-error recovery is no longer needed")
                p["published_error"] = TANGE_ERRORS[rid][p["parameter"]]
        notes.append(
            "Recovered missing ledger theta0/gamma0/q errors from primary Table 2, whose footnote explicitly identifies standard deviations. Source a=1 is fixed and b maps to q. Existing 42-row errors-in-variables fit and fixed V0 are unchanged. Source: "
            + TANGE_SOURCE
        )
    if rid.startswith(("rbbr_b2_farla_2025_", "rbi_b2_farla_2025_")):
        legacy = next(
            r
            for r in convert_eosmat(native, version=3)["eos_records"]
            if r["identifier"] == rid
        )
        residual, weighting = ordinary_residual(legacy_document(native), legacy, row)
        covariance = covariance_from_residual(
            residual, [p["refit"] for p in ps], [p["parameter"] for p in ps]
        )
        for p in ps:
            p["refit_error"] = covariance["standard_errors"][p["parameter"]]
        notes.append(
            weighting
            + " Conditional local errors use residual scaling and frozen diagonal weights, with fixed V0/theta0, all DAC and heated LVP rows, and source CsCl pressures. Original EoSFit row flags, calibration covariance and exact author workflow remain unrecovered."
        )
    if "_solomatova_2016_" in rid:
        if "_ls_" in rid:
            fixed = ["K0_prime"]
        if rid == "mg0490fe0483ti0027o_solomatova_2016_hs_bm3_reference_1":
            fixed = ["V0"]
        notes.append(
            "Target-branch agreement only, conditional on the existing coupled approximate-volume fit and Table 2 priors. Fixed LS K0-prime=4, and the Fp48 HS V0 anchor where applicable, are constraints, not independently recovered values. Prior regularization supports agreement and does not independently validate the endmembers or full spin-crossover fit. Reconstruction covariance remains unknown; this is a source-width central-value comparison, not combined-error or author-workflow parity."
        )
    if "_marquardt_2009b_ls_bm3_s1_" in rid:
        precision = {"K0": 0.5}
        notes.append(
            "Table S1 prints K0 as integer GPa: half-last-digit tolerance is 0.5 GPa, not an invented standard error. Only K0 is varied in the existing six-row cold LS fit, with each source V0/K0-prime fixed. No fit is preferred; uncertainty, exact source weighting and spin-mixture reproduction remain unestablished. Source tokens: peritheos/data/datasets/mg090fe010o-marquardt-2009-table-s1-ls-models.csv."
        )
    if rid == "rhenium_zha_2004_bm3_log_thermal":
        precision = {"alpha_KT_ref": 0.000005, "dK_dT_V": 0.000005}
        notes.append(
            "Both Table III thermal coefficients agree at five printed decimal places using the existing staged isochore-slope objective and fixed cold curve. Source errors and reconstruction covariance remain unknown. The eight Au pressure reductions and the 0.40849 GPa published-curve RMS remain separate unresolved findings; no measurement-error or author-workflow parity is claimed."
        )
    if rid == "magnesite_vocadlo_1999_bm3_2":
        precision = {"V0": 0.005}
        scope = "source_equation"
        notes.append(
            "Source-equation checkpoint coefficient recovery only: Table 1 pressures are calculated from the BM3 parameterization. Their near-exact refit is circular evidence for author-fit recovery, not an independent fit of measured pressures or the original E-V objective. V0=287.49 A^3 is checked at its 0.01 A^3 printed precision; the K0/K0-prime comparisons retain their source widths. The separately documented E-V fit and Table 1 cell-normalization correction remain distinct."
        )
    if rid in REFITS:
        scope = "independent_refit"
        if not ps:
            legacy = next(
                r
                for r in legacy_document(native)["eos_records"]
                if r["identifier"] == rid
            )
            published = legacy["eos"]["parameters"]
            ps = [
                dict(
                    parameter=n,
                    published=published[n],
                    published_error=None,
                    refit=row["refit"]["parameters"][n],
                    refit_error=None,
                    difference=row["refit"]["parameters"][n] - published[n],
                )
                for n in ("V0", "K0", "K0_prime")
            ]
        notes.append(
            "Replay of the registered independent Peritheos refit's stored point estimates, with numerical relative tolerance 1e-6 and absolute tolerance 1e-8 in the recorded units. This is solver repeatability, not agreement with source-author coefficients or an uncertainty estimate. Source-data limitations, graphical correlations and original-publication findings remain unchanged."
        )
    free = [p for p in ps if p["parameter"] not in fixed]
    source_width = bool(free) and all(
        p["published_error"] is not None
        and p["published_error"] > 0
        and abs(p["difference"]) <= p["published_error"]
        for p in free
    )
    combined = bool(free) and all(
        p["published_error"] is not None
        and p["published_error"] > 0
        and p["refit_error"] is not None
        and p["refit_error"] > 0
        and abs(p["difference"])
        <= 2 * math.hypot(p["published_error"], p["refit_error"])
        and numerical_screen(p)
        for p in free
    )
    printed = (
        bool(free)
        and bool(precision)
        and all(
            abs(p["difference"])
            <= (
                p["published_error"]
                if p["published_error"] is not None and p["published_error"] > 0
                else precision.get(p["parameter"], -1)
            )
            for p in free
        )
    )
    replay = (
        rid in REFITS
        and bool(free)
        and all(
            math.isclose(p["published"], p["refit"], rel_tol=1e-6, abs_tol=1e-8)
            for p in free
        )
    )
    if replay:
        basis = "reproduction_of_stored_peritheos_refit"
        reason = "Registered independent refit point estimates reproduce within the declared numerical solver tolerance."
    elif source_width:
        basis = "within_reported_parameter_errors"
        reason = "Every compared free coefficient lies inside its positive printed source error width."
    elif combined:
        basis = "conditional_combined_error_parameter_agreement"
        reason = "Every compared free coefficient meets both the combined-two-width diagnostic and numerical similarity screen. Source confidence and shared-data covariance are not inferred."
    elif printed:
        basis = "conditional_reported_widths_and_print_precision"
        reason = "Every compared free coefficient meets its printed source error width or the explicitly verified half-last-digit tolerance where source errors are absent. Printing precision is not a confidence interval."
    else:
        basis = None
        if not free:
            reason = "No complete free-coefficient comparison is available."
        elif any(p["published_error"] is None for p in free):
            reason = "At least one source parameter error is unavailable, and the complete comparison does not pass an explicitly verified printed-precision or registered-refit replay criterion."
        elif not all(numerical_screen(p) for p in free):
            reason = "The numerical similarity guard fails for at least one coefficient, and not every coefficient is inside its own source error width. This is a numerical limitation, not an author-workflow cap."
        else:
            reason = "At least one free coefficient fails the available combined-error comparison, or a required reconstruction error remains unavailable. Separate workflow qualifications cannot remove that numerical gap."
    if rid.startswith("palladium_frost_2023_"):
        notes.append(
            "The existing absolute-coordinate-error comparison fails for K0-prime. Residual-scaled covariance passes as a separate sensitivity but would change the uncertainty convention; no favorable rescaling is adopted for this review. Original-fit status and the explicit recommendation/default hold remain."
        )
    if rid == "fes_vi_morard_2026_bm3_mgd":
        notes.append(
            "Conditional staged cold-then-thermal parameter agreement under the stored EosFit diagnostic and reported widths, whose confidence is unspecified. GUI publication digits, cold-fit pressure scales and the full +/-3 GPa residual envelope remain separate qualifications. Use remains limited to the source's 40-150 GPa scope and observed thermal coverage; coefficient agreement does not establish phase stability or independently change the source-use assessment."
        )
    if rid == "naalsi2o6_zhao_1997_bm3_1":
        notes.append(
            "All four varied coefficients are within source widths on all 31 selected hydrostatic rows, with source V0=403 A^3 and K0-prime=5 fixed. The small alpha1 coefficient shifts about 46% but remains inside its printed error; this explicit source-width criterion is narrower in claim than the generic relative-change guard. Original weighting, objective, covariance and confidence remain unspecified."
        )
    limitations = []
    if basis is None:
        for p in free:
            name = p["parameter"]
            if p["published_error"] is None:
                if name in precision:
                    if abs(p["difference"]) > precision[name]:
                        limitations.append(
                            f"{name}: shift {abs(p['difference']):.6g} exceeds printed half-unit {precision[name]:.6g}"
                        )
                else:
                    limitations.append(name + ": source error unavailable")
            elif p["refit_error"] is None:
                if abs(p["difference"]) > p["published_error"]:
                    limitations.append(
                        name
                        + ": outside source width; reconstruction error unavailable"
                    )
            elif abs(p["difference"]) > 2 * math.hypot(
                p["published_error"], p["refit_error"]
            ):
                limitations.append(
                    f"{name}: shift {abs(p['difference']) / math.hypot(p['published_error'], p['refit_error']):.3g} combined widths exceeds 2"
                )
            elif not numerical_screen(p):
                limitations.append(name + ": numerical similarity guard not met")
    return dict(
        record_identifier=rid,
        material_identifier=row["material_identifier"],
        label=str(native.get("name", row["material_identifier"]))
        + ": "
        + str(row.get("label", rid)),
        original_status=row["status"],
        status="parity" if basis else "similar",
        scope=scope,
        parity_basis=basis,
        decision=reason,
        qualifications=notes,
        parameters=ps,
        free_comparison_parameters=[p["parameter"] for p in free],
        fixed_comparison_parameters=fixed,
        printed_half_units=precision,
        covariance=covariance,
        remaining_numerical_limits=limitations,
    )


def legacy_document(native):
    return convert_eosmat(native, version=3)


def calculate():
    verify_sources()
    inputs = json.loads(INPUT.read_text())
    return dict(
        format="peritheos.similar-parameter-scope-review",
        format_version=1,
        input_sha256=hashlib.sha256(INPUT.read_bytes()).hexdigest(),
        policy=POLICY,
        records=[review(row) for row in inputs["records"]],
    )


def check_report(saved, current):
    """Require literal evidence/decisions, allowing only numerical covariance roundoff."""
    left, right = copy.deepcopy(saved), copy.deepcopy(current)
    assert len(left["records"]) == len(right["records"])
    for a, b in zip(left["records"], right["records"]):
        if a["covariance"] is None and b["covariance"] is None:
            continue
        ca, cb = a.pop("covariance"), b.pop("covariance")
        for diagnostic, bound in (
            ("derivative_step_halving_relative_change", 2e-4),
            ("max_abs_residual_jacobian_cosine", 0.01),
        ):
            assert 0 <= ca.pop(diagnostic) < bound
            assert 0 <= cb.pop(diagnostic) < bound

        def close(x, y):
            if isinstance(x, dict):
                assert x.keys() == y.keys()
                for key in x:
                    close(x[key], y[key])
            elif isinstance(x, list):
                assert len(x) == len(y)
                for first, second in zip(x, y):
                    close(first, second)
            elif isinstance(x, float):
                assert math.isclose(x, y, rel_tol=3e-4, abs_tol=1e-25)
            else:
                assert x == y

        close(ca, cb)
        for pa, pb in zip(a["parameters"], b["parameters"]):
            close(pa.pop("refit_error"), pb.pop("refit_error"))
    assert left == right


def render_report(report):
    counts = Counter(r["status"] for r in report["records"])
    scopes = Counter(r["scope"] for r in report["records"] if r["status"] == "parity")
    lines = [
        "# Remaining similar records: parameter-scope review",
        "",
        "Reviewed on 25 October 2026, after the separate 63-record uncertainty follow-up.",
        "",
        f"All **{len(report['records'])} remaining similar records** were checked. **{counts['parity']}** now meet an explicitly scoped parity criterion and **{counts['similar']}** retain similarity.",
        f"The promotions comprise **{scopes['parameter_agreement']} conditional parameter-agreement results**, **{scopes['independent_refit']} registered independent-refit replays**, and **{scopes['source_equation']} source-equation checkpoint check**. The latter three are not source-author parameter-refit parity.",
        "",
        "## Criterion and preserved evidence",
        "",
        POLICY,
        "",
        "The source-width criterion compares central values with each positive printed source width; it is a numerical criterion even where reconstruction covariance remains unavailable. The combined-width criterion requires both finite positive errors and the numerical similarity guard. Verified printed precision is a separate nonstatistical criterion, never a substitute error bar. Fixed values are identified explicitly and are not counted as recovered free parameters.",
        "",
        "Source coefficients, printed errors, observations, calibration reductions, models, fit masks, defaults, source-validation findings and use holds are unchanged. Only reconstruction assessments and the derived ledger/report are updated. Original ledger findings are frozen in [the review inputs](../data/similar-parameter-review-inputs.json); source-document and dataset hashes protect them.",
        "",
        "## Recovered information",
        "",
        f"Tange (2012) Table 2 explicitly reports standard deviations for theta0, gamma0 and b (mapped to q for fixed a=1). These were present in source notes but absent from ledger comparisons. [Primary Table 2]({TANGE_SOURCE}) confirms BM3 errors 150 K, 0.09 and 0.3; Vinet errors 140 K, 0.09 and 0.3. The existing 42-row fits, source coefficients and reconstruction errors are unchanged.",
        "",
        "For Farla RbBr and RbI, missing conditional covariance was calculated at the existing converged coefficients with the original final frozen effective P/V/T weights. Row counts, pressure RMS and local stationarity replay before errors are accepted. Fixed V0/theta0, selected DAC/heated LVP rows, printed CsCl pressures and calibration/workflow qualifications remain unchanged. The full covariance matrices and assumptions are in [the numerical review](../data/similar-parameter-review.json).",
        "",
        "## Promoted scoped comparisons",
        "",
        "| Record | Scope and criterion | Conditions and separate findings |",
        "| --- | --- | --- |",
    ]
    for row in report["records"]:
        if row["status"] != "parity":
            continue
        rid = row["record_identifier"]
        link = f"../reconstruction-reports/{row['material_identifier']}.md#record-{rid}"
        notes = " ".join(row["qualifications"]).replace("|", "\\|")
        label = row["label"].replace("|", "\\|")
        lines.append(
            f"| [{label}]({link}) | {row['scope']}; `{row['parity_basis']}` | {notes} |"
        )
    lines.extend(
        [
            "",
            "## Retained similarity",
            "",
            "These are remaining numerical or uncertainty gaps, rather than a blanket author-workflow cap. Existing scientific limitations remain in each native assessment and source audit. In particular, Frost Pd stays similar under its existing absolute-coordinate-error comparison; residual-scaled covariance is a separately passing sensitivity and is not silently adopted. The explicit recommendation/default hold remains.",
            "",
            "| Record | Remaining comparison limit |",
            "| --- | --- |",
        ]
    )
    for row in report["records"]:
        if row["status"] != "similar":
            continue
        rid = row["record_identifier"]
        link = f"../reconstruction-reports/{row['material_identifier']}.md#record-{rid}"
        limits = "; ".join(row["remaining_numerical_limits"]) or row["decision"]
        label = row["label"].replace("|", "\\|")
        lines.append(f"| [{label}]({link}) | {limits} |")
    lines.extend(
        [
            "",
            "## Reproduction",
            "",
            "Repository maintenance only:",
            "",
            "```bash",
            "uv run python scripts/review_similar_parameter_scope.py --check",
            "uv run pytest -q tests/test_similar_parameter_scope_review.py",
            "```",
            "",
        ]
    )
    return "\n".join(lines)


@lru_cache(maxsize=1)
def saved_review():
    verify_sources()
    report = json.loads(OUTPUT.read_text())
    if report["input_sha256"] != hashlib.sha256(INPUT.read_bytes()).hexdigest():
        raise ValueError("Parameter-scope review input changed")
    return {r["record_identifier"]: r for r in report["records"]}


def augment_outcome(record, outcome):
    case = saved_review().get(record["identifier"])
    if case is None:
        return outcome
    if case["status"] == "similar":
        baseline = next(
            r
            for r in json.loads(INPUT.read_text())["records"]
            if r["record_identifier"] == record["identifier"]
        )
        if outcome["status"] != baseline["status"]:
            raise ValueError(
                "Unpromoted numerical finding changed: " + record["identifier"]
            )
        fresh = {p["parameter"]: p for p in outcome.get("parameters", [])}
        for old in baseline.get("parameters", []):
            current = fresh.get(old["parameter"])
            if (
                current is None
                or any(
                    current.get(key) != old.get(key)
                    for key in ("published", "published_error")
                )
                or not math.isclose(
                    current["refit"], old["refit"], rel_tol=3e-4, abs_tol=1e-25
                )
            ):
                raise ValueError(
                    "Unpromoted comparison changed: " + record["identifier"]
                )
        # Preserve the reviewed snapshot, including its qualifications and diagnostics;
        # insignificant platform roundoff does not create a new scientific finding.
        return {**copy.deepcopy(outcome), **copy.deepcopy(baseline)}
    if outcome["status"] not in {"similar", "parity"}:
        raise ValueError(
            "Original numerical assessment changed scope: " + record["identifier"]
        )
    # Local solves can differ in weak small coefficients across platforms.
    # The named scientific acceptance criterion is checked separately below.
    before = {p["parameter"]: p for p in outcome.get("parameters", [])}
    for p in case["parameters"]:
        original = before.get(p["parameter"])
        if original and (
            original["published"] != p["published"]
            or not math.isclose(
                original["refit"], p["refit"], rel_tol=3e-4, abs_tol=1e-25
            )
        ):
            raise ValueError(
                "Source or fitted point estimate changed: " + record["identifier"]
            )
        if original and p["parameter"] in case["free_comparison_parameters"]:
            difference = abs(original["refit"] - p["published"])
            if case["parity_basis"] == "within_reported_parameter_errors":
                accepted = difference <= p["published_error"]
            elif (
                case["parity_basis"] == "conditional_combined_error_parameter_agreement"
            ):
                accepted = difference <= 2 * math.hypot(
                    p["published_error"], p["refit_error"]
                ) and numerical_screen({**p, "refit": original["refit"]})
            elif (
                case["parity_basis"]
                == "conditional_reported_widths_and_print_precision"
            ):
                width = p["published_error"]
                if width is None:
                    width = case["printed_half_units"][p["parameter"]]
                accepted = difference <= width
            else:
                accepted = math.isclose(
                    original["refit"], p["published"], rel_tol=1e-6, abs_tol=1e-8
                )
            if not accepted:
                raise ValueError(
                    "Fresh reconstruction no longer meets its scoped criterion: "
                    + record["identifier"]
                )
    updated = copy.deepcopy(outcome)
    baseline = next(
        r
        for r in json.loads(INPUT.read_text())["records"]
        if r["record_identifier"] == record["identifier"]
    )
    history = " ".join(
        str(baseline[k]) for k in ("qualification", "reason") if baseline.get(k)
    )
    updated.update(
        status="parity",
        parity_basis=case["parity_basis"],
        parameters=case["parameters"],
        qualification=" ".join([case["decision"], *case["qualifications"]])
        + " Earlier audit wording (superseded only for the named numerical scope): "
        + history,
        comparison_scope="Compared free coefficients: "
        + ", ".join(case["free_comparison_parameters"])
        + "; scope="
        + case["scope"],
    )
    if case["fixed_comparison_parameters"]:
        updated["validation_fixed_parameters"] = case["fixed_comparison_parameters"]
    if case["covariance"] is not None:
        updated["covariance"] = case["covariance"]
        updated["covariance_convention"] = case["covariance"]["method"]
    if case["scope"] == "source_equation":
        updated["fit_kind"] = "source_curve_validation_only"
        updated["reconstruction_kind"] = "published_calculated_pressure_checkpoints"
    return updated


def integrate():
    from scripts.standardize_reconstruction_metadata import assessment_for_record
    from scripts.validate_primary_eos_refits import render_markdown

    ledger_path = ROOT / "docs/data/primary-eos-refits.json"
    ledger = json.loads(ledger_path.read_text())
    documents = {}
    for row in ledger["records"]:
        case = saved_review().get(row["record_identifier"])
        if case is None or case["status"] != "parity":
            continue
        material = row["material_identifier"]
        path = ROOT / "peritheos/data/materials" / (material + ".eosmat")
        if material not in documents:
            documents[material] = load_eosmat(path)
        document = documents[material]
        record = next(
            r
            for r in document["eos_records"]
            if r["identifier"] == row["record_identifier"]
        )
        row.update(augment_outcome(record, row))
        record["evidence"]["reproduction"]["assessment"] = assessment_for_record(
            record, row
        )
    for material, document in documents.items():
        validate_eosmat_document(document)
        path = ROOT / "peritheos/data/materials" / (material + ".eosmat")
        original, formatted = path.read_text(), format_eosmat(document)
        for record in document["eos_records"]:
            case = saved_review().get(record["identifier"])
            if not case or case["status"] != "parity":
                continue
            marker = '      "identifier": "' + record["identifier"] + '",'

            def bounds(text):
                at = text.index(marker)
                return text.rfind("    {", 0, at), text.index("\n    }", at) + len(
                    "\n    }"
                )

            a, b = bounds(original)
            c, d = bounds(formatted)
            original = original[:a] + formatted[c:d] + original[b:]
        path.write_text(original)
    ledger["summary"] = {
        "total": len(ledger["records"]),
        **dict(sorted(Counter(r["status"] for r in ledger["records"]).items())),
    }
    if POLICY not in ledger["policy"]["parity"]:
        ledger["policy"]["parity"] += " " + POLICY
    ledger_path.write_text(
        json.dumps(ledger, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    (ROOT / "docs/primary-eos-refits.md").write_text(render_markdown(ledger))
    verify_sources()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--integrate", action="store_true")
    args = parser.parse_args()
    report = calculate()
    if args.check:
        previous = json.loads(OUTPUT.read_text())
        try:
            check_report(previous, report)
        except AssertionError as error:
            raise SystemExit("Parameter-scope review is stale") from error
        if REPORT.read_text() != render_report(report):
            raise SystemExit("Parameter-scope review documentation is stale")
    else:
        OUTPUT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        REPORT.write_text(render_report(report))
    counts = Counter(r["status"] for r in report["records"])
    print(dict(counts))
    for row in report["records"]:
        if row["status"] == "parity":
            print(row["record_identifier"], row["scope"], row["parity_basis"])
    if args.integrate:
        integrate()


if __name__ == "__main__":
    main()
