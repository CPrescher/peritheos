#!/usr/bin/env python3
"""Generate the paper-level investigation ledger from audited record evidence."""

from __future__ import annotations

import argparse
import json
import re
from collections import Counter, defaultdict
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
AUDIT_PATH = ROOT / "peritheos" / "data" / "primary-source-audit.json"
REFIT_PATH = ROOT / "docs" / "data" / "primary-eos-refits.json"
NONPRODUCTION_PATH = ROOT / "docs" / "data" / "nonproduction-paper-investigations.json"
OUTPUT_PATH = ROOT / "docs" / "paper-investigation-ledger.md"

STATUS_LABELS = {
    "parity": "parity",
    "similar": "similar",
    "parity_not_achieved": "parity not achieved",
    "not_refittable": "direct refit unavailable",
    "source_reconstruction": "source reconstruction",
}

OUTCOME_LABELS = {
    "unreproduced_original_with_independent_refit": "Original unreproduced; independent refit available",
    "reproduced": "Reproduced",
    "partial_reproduction": "Partly reproduced",
    "mixed_reproduction": "Mixed: reproduced and discrepant records",
    "parity_not_achieved": "Coefficient parity not achieved",
    "direct_refit_unavailable": "Direct refit unavailable",
    "source_reconstruction": "Source reconstruction",
    "withheld_unreproduced": "Withheld: could not reproduce",
    "deferred_incomplete_model": "Deferred: incomplete source/model mapping",
    "deferred_batch_boundary": "Deferred: complete family queued for next batch",
}

CITATION_OVERRIDES = {
    "10.1029/jb095ib13p21737": "Mao et al. (1990)",
    "10.2138/am-2002-2-316": "Crichton et al. (2002)",
    "10.2138/am.2006.2347": "Ono et al. (2006)",
    "10.2138/am-2000-11-1229": "Fei et al. (2000)",
    "10.1073/pnas.0609013104": "Fei et al. (2007)",
    "10.1103/physrevlett.74.1371": "Duffy et al. (1995)",
}


PAPER_SCOPE_NOTES = {
    "10.1007/s00269-003-0353-z": (
        "Full published thermal EOS: parity not achieved, not recommended. "
        "Only individually source-checked 300/473/673 K isotherms remain executable. "
        "The original 873 K normalization and continuous source reconstruction are unresolved. "
        "All seven Peritheos refits are historical diagnostics, withheld from executable "
        "selection; numerical replay of those fits is not publication parity or a "
        "validated replacement. Further model development is outside this curation task. "
        "See the [audit](literature-reproductions/brucite.md)."
    ),
    "10.1180/002646100549427": (
        "Fe-bearing phase E volume EOS successfully recovered: parity within "
        "reported uncertainty on all 11 primary rows. Published coefficients "
        "are preserved. Largest weighted-refit shift is 0.123 reported error "
        "widths; maximum curve difference is about 0.0096 GPa. Exact original "
        "workflow availability is separate provenance, not a failed-fit "
        "classification. See the [audit, coefficient differences and residual "
        "plot](literature-reproductions/phase-e-crichton-ross-2000.md)."
    ),
    "10.1029/2021jb023562": (
        "Conditional numerical parity for the three free forsterite coefficients "
        "using all 51 reconstructed observations and printed Table 2 error widths. "
        "The recommended thermal record is "
        "`forsterite_katsura_2022_author_script_bm3_mgd`: full-precision author-script "
        "K0-prime=4.21116896988574, gamma0=0.997258584556838, "
        "q=2.43854091371459, theta0=760 K; V0=289.58 A3 and K0=127.4 GPa are fixed. "
        "The theta0=768 K printed-table alternative is retained. Table errors "
        "are comparison widths only for the script vector; script-specific errors "
        "and exact bootstrap reproduction remain unavailable. The original "
        "2009 combined fit remains a separate unresolved investigation. See the "
        "[source comparison](literature-reproductions/forsterite-katsura-2022-refit-comparison.md)."
    ),
    "10.1103/physrevb.78.104119": (
        "Conditional numerical K0-prime parity on both ruby coordinates using "
        "all 86 observed rows with V0/K0 fixed. Original weights and uncertainty "
        "estimation remain unresolved; no exact author-fit claim. See the "
        "[dedicated audit](literature-reproductions/takemura-2008-gold.md)."
    ),
    "10.1103/physrevb.80.104114": (
        "The published fitted 300 K Au BM3/Vinet and Pt Vinet branches are "
        "bundled with density-derived four-atom volume normalizations. "
        "Their equations are independently reproduced. The Pt Vinet branch "
        "resolves Sakai (2018)'s named reference; "
        "paired calibrant observations remain missing. All 168 Pt Table V "
        "states and all 162 populated Au Table III states are full-model output, "
        "not fitted-isotherm observations. Au liquid markers and blank cells "
        "are retained. The Au thermal pressure reconstruction is numerically "
        "verified against derived Table III output (156 unmarked states, RMS "
        "0.008923 GPa); its coefficients are not the published parameter set. "
        "Both Au and Pt bounded PVT reconstructions are now registered in the "
        "library; Pt additionally retains a derived residual-pressure correction "
        "and agrees with unmarked Table V states within 0.00594 GPa. "
        "These two records are nondefault, explicit-selection diagnostics with "
        "scientific validation status not_reproduced; source reconstruction here "
        "means numerical table-output verification, not validation of the published "
        "analytical PVT EOS. The [pressure-convention audit]"
        "(literature-reproductions/yokoo-2009-pressure-conventions.md) tests the "
        "remaining thermal-increment mismatch without empirical corrections. "
        "The Au/Pt PVT investigation was closed for now on 2026-10-09, pending "
        "author inputs; the published analytical EOS remains unreproduced. "
        "Reopen when the authors' Table III/V pressure code or spreadsheet, or an "
        "equivalent complete specification/correction, resolves the discrepancy. "
        "Required inputs are the adopted phonon parameters and conventions, "
        "electronic-pressure routine and reference conventions, and cold Vc "
        "separately from ambient V0 for absolute pressure. Missing electronic "
        "energy concerns the original shock/caloric calculation, not bounded "
        "pressure evaluation. See the [Au PVT validation]"
        "(literature-reproductions/gold-pvt-validation.md), "
        "[Pt audit](literature-reproductions/yokoo-2009-platinum.md) and "
        "[Au audit](literature-reproductions/yokoo-2009-gold.md)."
    ),
    "10.1073/pnas.0609013104": (
        "Conditional Au numerical parity reproduces K0-prime/q within printed "
        "widths using a staged equal-pressure-weight replay; source weights, "
        "exact MgO reduction and error confidence/covariance remain unresolved. "
        "See the [Au replay](literature-reproductions/fei-2007-gold.md). "
        "The Pt Vinet-MGD equation is independently reproduced, including the "
        "printed variable-exponent Debye law. A 78-row joint diagnostic "
        "reconstructs the upstream Au pressure reduction but does not recover "
        "the published thermal coefficients; source weights and constraints "
        "remain unspecified. The complete published Au PVT equation is independently "
        "verified, with Figure 1 vector-curve agreement at graphical precision "
        "and recovered hot MgO/cold ruby comparisons; see the [Au PVT audit]"
        "(literature-reproductions/gold-pvt-validation.md). See the [Pt audit]"
        "(literature-reproductions/fei-2007-platinum.md)."
    ),
    "10.1063/1.3054331": (
        "The complete Vinet/integrated-q Debye plus electronic pressure scale "
        "is registered using 51 recovered Tsuchiya-Kawamura Pt nodes; the "
        "erratum corrects Au only. Linear interpolation is explicitly a numerical "
        "implementation choice. All 35 source pressure-grid states agree within "
        "0.0071 GPa. Original shock/expansion sampling, weights and joint-fit "
        "objective remain unresolved; equation reconstruction is distinct from "
        "author-fit parity. See the "
        "[Matsui audit](literature-reproductions/matsui-2009-platinum.md)."
    ),
    "10.1029/2024jb028819": (
        "Reproduction covers Ma's final reduced acoustic/Walker inputs and "
        "model grids; source weights/covariance and exact original-lattice "
        "Walker recalibration remain unresolved. High-temperature coverage is "
        "modeled. See the [Ma audit](literature-reproductions/ma-2024-kcl.md)."
    ),
    "10.1103/h4pj-rvxx": (
        "FeS source record accepted unchanged on parameter agreement within "
        "reported uncertainties (maintainer decision, 2026-10-19). Table S1, paired "
        "KCl markers, the selected 167-row author input and the reference cold "
        "datasets are bundled separately. Author-input staged fit reproduction is classified similar; the residual-range discrepancy is retained separately. "
        "Distributed-GUI fits agree with the console; the residual discrepancy "
        "is informational. Exact original workflow remains unverified; see the "
        "[source audit](literature-reproductions/morard-2026-fes.md)."
    ),
    "10.1103/physrevb.81.144110": (
        "DO NOT USE the original-constant reconstruction. Separate Peritheos density refit; "
        "[comparison and explanation](#chen-2010-original-constant-reconstruction-and-independent-refit)."
    ),
    "10.1002/2013jb010898": (
        "Two published B2 FeSi candidates withheld: unresolved thermal normalization; two independent Peritheos refits added; "
        "see [source audit](literature-reproductions/fischer-2014-fesi.md)."
    ),
}

# Dedicated source-protocol audits supersede the older generic regressions.
# Keep the historical record ledger, but do not turn refit reproducibility
# into a claim that every published author fit was reproduced.
PAPER_REPRODUCTION_OVERRIDES = {
    "10.1103/fxgq-96sg": {
        "outcome": "partial_reproduction",
        "summary": "Eight phase fits retain conditional coefficient parity. MgO and NaCl-B2 are not_reproduced and not recommended: the published curves give 1.161 and 0.622 GPa pressure RMS, versus 0.476 and 0.293 GPa for direct fixed-V0 diagnostic refits. Table S1 accounts for all 194 published DAC-2 MgO observations and 131 NaCl-B2 observations (125 above 35 GPa). Exact author fit-input completeness remains unconfirmed: the fitting description identifies this study's data, but data availability mentions other author-held data. See the [source-input audit](literature-reproductions/shen-smith-2026-cu-refits.md).",
    },
    "10.1063/5.0161038": {
        "outcome": "parity_not_achieved",
        "summary": "Final-paper protocol recheck: both Pd refits remain numerically similar; the two-width uncertainty comparison fails with absolute-error covariance but passes with residual-scaled covariance. Author objective, ambient treatment, constraints and covariance convention remain unrecovered. Scientific status not_reproduced denotes an unrecovered original fit, not a demonstrated failed author-protocol replay. Quantitative recommendation and Pd automatic default remain on hold pending clarification. The earlier adjusted-volume RMS comparison is corrected in the dedicated Frost audit.",
    },
    "10.2138/am-2022-8435": {
        "outcome": "partial_reproduction",
        "summary": "Bridgmanite fits retained; published CaSiO3 thermal parameterization not reproduced and not recommended. Its historical _refit identifier denotes Fu's published coefficients, not the separately attributed Peritheos refit.",
    },
    "10.1063/5.0179469": {
        "outcome": "partial_reproduction",
        "summary": "Experimental Pd BM3 fit not reproduced and not recommended; separate 0 K DFT record preserved. Pd automatic default selection is on hold after the separate Frost audit.",
    },
    "10.1016/j.pepi.2013.12.010": {
        "outcome": "partial_reproduction",
        "summary": "Five Fe-Ni thermal records not recommended: about 7.1 GPa RMS on the P7 hot data versus about 1.1 GPa reported. Separate isotherms retain their source-specific qualifications; Table 8 checkpoints do not resolve the thermal-data discrepancy.",
    },
    "10.1029/2009gl038107": {
        "outcome": "withheld_unreproduced",
        "summary": "Wadsleyite default on hold and deferred pending pressure-calibration and MGD mapping review; published coefficients and all 85 source rows preserved, no replacement default selected.",
    },
    "10.1007/s00269-003-0353-z": {
        "outcome": "partial_reproduction",
        "summary": "Three individual isotherms source-checked; full published thermal EOS not reproduced and not recommended. Unresolved source curves and all independent refits withheld as diagnostic evidence",
    },
    "10.1029/95jb03254": {
        "outcome": "partial_reproduction",
        "summary": (
            "1 rounded source fit reproduced; room-temperature parameters within "
            "source errors; pressure-weighted source fits unresolved; independent "
            "thermal refit reproduced. See the "
            "[dedicated audit](literature-reproductions/wang-1996-casio3.md)."
        ),
    },
}


def load_json(path: Path) -> dict[str, Any]:
    """Load one UTF-8 JSON document."""
    return json.loads(path.read_text(encoding="utf-8"))


def paper_key(audit: dict[str, Any]) -> tuple[str, str]:
    """Return a stable paper key, preferring DOI over access URL."""
    doi = audit.get("doi")
    if doi:
        return "doi", str(doi).lower()
    check = audit.get("primary_source_check", {})
    return "url", str(check["access_url"])


def citation_from_labels(labels: list[str]) -> str:
    """Choose a compact author-year citation from record-specific labels."""
    candidates = []
    for label in labels:
        match = re.match(r"(.+?\(\d{4}\))", label)
        candidates.append(match.group(1) if match else label)
    return min(candidates, key=lambda value: (len(value), value.casefold()))


def markdown_cell(value: str) -> str:
    """Escape content that would otherwise split a Markdown table cell."""
    return value.replace("|", "\\|").replace("\n", " ")


def classify(statuses: Counter[str]) -> str:
    """Collapse record-level refit results into one paper disposition."""
    reproduced = statuses["parity"] + statuses["similar"]
    discrepant = statuses["parity_not_achieved"]
    unavailable = statuses["not_refittable"]
    reconstructed = statuses["source_reconstruction"]
    if discrepant and reproduced:
        return "mixed_reproduction"
    if discrepant:
        return "parity_not_achieved"
    if reproduced and (unavailable or reconstructed):
        return "partial_reproduction"
    if reproduced:
        return "reproduced"
    if reconstructed:
        return "source_reconstruction"
    return "direct_refit_unavailable"


def source_link(paper: dict[str, Any]) -> str:
    """Render the primary source link for a paper."""
    if paper.get("doi"):
        target = f"https://doi.org/{paper['doi']}"
    else:
        target = paper["access_url"]
    return f"[{paper['citation']}]({target})"


def compact_counts(counts: Counter[str], labels: dict[str, str]) -> str:
    """Render nonzero status counts in a stable order."""
    return "; ".join(
        f"{count} {labels[key]}" for key in labels if (count := counts.get(key, 0))
    )


def record_list(records: list[dict[str, Any]], status: str | None = None) -> str:
    """Render record identifiers, optionally filtered by refit status."""
    identifiers = [
        row["record_identifier"]
        for row in records
        if status is None or row["status"] == status
    ]
    return ", ".join(f"`{identifier}`" for identifier in sorted(identifiers))


def discrepancy_summary(record: dict[str, Any]) -> str:
    """Summarize published-to-refit changes for a discrepant record."""
    differences = []
    for parameter in record.get("parameters", []):
        if not parameter.get("similar") and not parameter.get("within_combined_2sigma"):
            differences.append(
                f"{parameter['parameter']} {parameter['published']:.6g} -> "
                f"{parameter['refit']:.6g}"
            )
    return (
        "; ".join(differences)
        or record.get("reason")
        or "published coefficients were not recovered"
    )


def build_papers() -> tuple[list[dict[str, Any]], dict[str, int]]:
    """Combine source audit, refit results, and nonproduction investigations."""
    audits = load_json(AUDIT_PATH)["records"]
    refits = {row["record_identifier"]: row for row in load_json(REFIT_PATH)["records"]}
    audit_ids = {row["record"] for row in audits}
    if audit_ids != set(refits):
        missing_audit = sorted(set(refits) - audit_ids)
        missing_refit = sorted(audit_ids - set(refits))
        raise ValueError(
            "audit/refit record mismatch: "
            f"missing audit={missing_audit}, missing refit={missing_refit}"
        )

    grouped: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    audit_by_key: dict[tuple[str, str], list[dict[str, Any]]] = defaultdict(list)
    for audit in audits:
        key = paper_key(audit)
        grouped[key].append(refits[audit["record"]])
        audit_by_key[key].append(audit)

    papers = []
    for key, records in grouped.items():
        source_rows = audit_by_key[key]
        labels = sorted({row["label"] for row in source_rows})
        statuses = Counter(row["status"] for row in records)
        data_statuses = Counter(row["primary_data_status"] for row in records)
        doi = next((row.get("doi") for row in source_rows if row.get("doi")), None)
        access_url = source_rows[0]["primary_source_check"]["access_url"]
        citation = CITATION_OVERRIDES.get(str(doi).lower()) if doi else None
        papers.append(
            {
                "citation": citation or citation_from_labels(labels),
                "doi": doi,
                "access_url": access_url,
                "outcome": (
                    "unreproduced_original_with_independent_refit"
                    if any(
                        row.get("original_publication_reproduction_status")
                        == "not_reproduced"
                        and row.get("reproduction_status")
                        in {
                            "independent_digitized_data_refit_reproduced",
                            "independent_refit_reproduced",
                        }
                        for row in records
                    )
                    else classify(statuses)
                ),
                "statuses": statuses,
                "data_statuses": data_statuses,
                "records": sorted(records, key=lambda row: row["record_identifier"]),
                "withheld_candidates": [
                    {**candidate, "evidence": row["evidence"]}
                    for row in source_rows
                    for candidate in row.get("withheld_candidates", [])
                ],
            }
        )
        override = PAPER_REPRODUCTION_OVERRIDES.get(str(doi).lower())
        if override:
            papers[-1]["outcome"] = override["outcome"]
            papers[-1]["reproduction_summary"] = override["summary"]

    for paper in load_json(NONPRODUCTION_PATH)["papers"]:
        papers.append({**paper, "records": []})

    papers.sort(key=lambda row: row["citation"].casefold())
    outcome_counts = Counter(paper["outcome"] for paper in papers)
    return papers, dict(outcome_counts)


def chen_comparison_section() -> str:
    """Keep both ledgers explicit about the separate meanings of the curves."""
    report = load_json(ROOT / "docs/data/chen-2010-argon-comparison.json")
    rows = [
        "## Chen 2010: original-constant reconstruction and independent refit",
        "",
        "**Original-constant reconstruction: DO NOT USE for scientific predictions or pressure calibration.** "
        "The record `argon_fcc_chen_2010_bm3_reported_constants` is retained for explicit diagnostic inspection. "
        "It is not a verified author equation: standard BM3 is reconstructed from the rounded "
        "P = 2 GPa, density = 2.18 g/cm³, KT = 15.1 GPa and KT′ = 5.4 constraints.",
        "",
        "The separate `argon_fcc_chen_2010_bm3_digitized_refit` fits all 80 distinct digitized "
        "Brillouin-integrated density positions. It does not use the published fitted-line vertices "
        "or printed elastic constants. This numerical refit does **not** reproduce the original EOS.",
        "",
        "![Chen published curve, printed-constant reconstruction and independent refit](data/chen-2010-argon-comparison.png)",
        "",
        "| At the same pressure | Published plotted curve: volume difference from refit | Printed-constant reconstruction: volume difference from refit |",
        "|---|---:|---:|",
    ]
    for row in report["checkpoints"]:
        delta = row["volume_difference_from_refit_percent"]
        rows.append(
            f"| {row['pressure_gpa']} GPa | {delta['published_curve']:+.2f}% | {delta['reported_constants']:+.2f}% |"
        )
    rows.extend(
        [
            "",
            "Differences are 100 × (Vcomparison / Vrefit − 1). The printed-constant reconstruction "
            "has 6.880 GPa pressure RMS error at the original marker coordinates; the refit has "
            "0.111 GPa. The paper's drawn curve closely follows the refit, with a maximum volume "
            "difference of 1.40% over 1.23–26.06 GPa. The alternative pressure-offset convention "
            "also fails the density comparison.",
            "",
            "The 2 GPa constants are not zero-pressure coefficients. The density integration instead "
            "starts at 1.3 GPa. At 2 GPa the drawn curve's slope implies KT ≈ 11 GPa, the refit gives "
            "12.62 GPa, and the paper prints 15.1 GPa. The cause is unresolved: **no error in Chen's "
            "pressure calibration or adiabatic-to-isothermal correction has been established**.",
            "",
            "The refit is provisional and qualified only within the digitized interval at nominal 290 K. "
            "Its extrapolated rho0 = 1.6895 g/cm³ is higher than Chen's 1.52 ± 0.05 and the 1.564 "
            "room-temperature estimate Chen attributes to low-temperature measurements in ref. 22. "
            "We have not independently verified that temperature extrapolation. Parameter covariance "
            "and statistical uncertainties are unavailable; numerical reproducibility is not independent physical validation.",
            "",
            "[Full source audit, equations, fit method and limitations](literature-reproductions/chen-2010-argon.md) · "
            "[Comparison data](data/chen-2010-argon-comparison.json) · "
            "[Refit report](data/chen-2010-argon-refit.json).",
            "",
        ]
    )
    return "\n".join(rows)


def render() -> str:
    """Render the complete Markdown ledger."""
    papers, outcome_counts = build_papers()
    catalog_papers = [paper for paper in papers if paper["records"]]
    nonproduction = [paper for paper in papers if not paper["records"]]
    record_count = sum(len(paper["records"]) for paper in catalog_papers)
    discrepancy_papers = [
        paper
        for paper in catalog_papers
        if any(record["status"] == "parity_not_achieved" for record in paper["records"])
    ]
    discrepancy_record_count = sum(
        record["status"] == "parity_not_achieved"
        for paper in discrepancy_papers
        for record in paper["records"]
    )
    unavailable_papers = [
        paper
        for paper in catalog_papers
        if any(record["status"] == "not_refittable" for record in paper["records"])
    ]
    unavailable_record_count = sum(
        record["status"] == "not_refittable"
        for paper in unavailable_papers
        for record in paper["records"]
    )
    lines = [
        "# Paper investigation ledger",
        "",
        "Generated from the primary-source audit, the record-level refit ledger, and",
        "the explicit nonproduction investigation register. This page answers a",
        "different question from the EOS catalog: it records what happened to every",
        "primary paper that was actually investigated, including papers that did not",
        "produce an executable record.",
        "",
        "## Status definitions",
        "",
        "- **Reproduced:** every executable record from the paper reached `parity` or",
        "  `similar` in the documented independent check.",
        "- **Partly reproduced:** at least one record was reproduced, while another",
        "  has only a source reconstruction or could not be refitted directly.",
        "- **Coefficient parity not achieved:** the refit ran, but at least one",
        "  published coefficient was outside both the uncertainty and numerical",
        "  similarity criteria, or a coupled source-level objective had a demonstrably",
        "  different optimum. These are source-fit discrepancies, not software-run",
        "  failures; the record may remain for faithful published-curve provenance.",
        "- **Direct refit unavailable:** the equation and parameters were audited, but",
        "  independent coefficient recovery was impossible because primary rows, an",
        "  executable calibration, or the original reduction were unavailable or",
        "  circular. Qualified subset or proxy fits may still provide diagnostics.",
        "- **Source reconstruction:** a source calculation or composition of audited",
        "  equations is reproduced without independently fitting the complete EOS.",
        "  Component refits and remaining source-input gaps are documented separately.",
        "- **Original unreproduced; independent refit available:** the original source EOS",
        "  remains unresolved even when a separate Peritheos refit can be reproduced numerically.",
        "- **Withheld/deferred:** investigation did not pass the executable-record",
        "  acceptance gate, so no production EOS was added.",
        "",
        "## Summary",
        "",
        f"The register covers **{len(papers)} primary papers**: "
        f"**{len(catalog_papers)}** support the {record_count} audited catalog records and "
        f"**{len(nonproduction)}** were investigated without adding a production record.",
        "No numerical refit attempt failed before producing a comparison. The adverse",
        "outcomes are instead explicit coefficient discrepancies, unavailable direct",
        "refits, or acceptance-gate holds.",
        "",
        "| Paper-level outcome | Papers |",
        "|---|---:|",
    ]
    for outcome in OUTCOME_LABELS:
        count = outcome_counts.get(outcome, 0)
        if count:
            lines.append(f"| {OUTCOME_LABELS[outcome]} | {count} |")

    lines.extend(["", chen_comparison_section()])

    lines.extend(
        [
            "",
            "## Withheld or deferred papers",
            "",
        ]
    )
    for paper in nonproduction:
        evidence = paper["evidence"]
        lines.extend(
            [
                f"### {source_link(paper)}",
                "",
                f"**Outcome:** {OUTCOME_LABELS[paper['outcome']]} "
                f"({paper['investigation_date']}).",
                "",
                paper["reason"],
                "",
                f"Evidence: [{evidence}]({evidence}).",
                "",
            ]
        )

    lines.extend(
        [
            "",
            "## Papers with coefficient discrepancies",
            "",
            f"These **{len(discrepancy_papers)} papers** account for all "
            f"{discrepancy_record_count} records",
            "classified as `parity_not_achieved`. Papers with other successful records",
            "are marked as mixed in the complete register.",
            "",
            "| Paper | Affected record | Published-to-refit discrepancy |",
            "|---|---|---|",
        ]
    )
    for paper in discrepancy_papers:
        first = True
        for record in paper["records"]:
            if record["status"] != "parity_not_achieved":
                continue
            paper_cell = source_link(paper) if first else ""
            first = False
            lines.append(
                f"| {paper_cell} | `{record['record_identifier']}` | "
                f"{markdown_cell(discrepancy_summary(record))} |"
            )

    lines.extend(
        [
            "",
            "The full data selection, model mapping, residuals, and bounded explanation",
            "for every row above are in the",
            "[primary EOS refit validation](primary-eos-refits.md#detailed-non-parity-investigations).",
            "",
            "## Papers with unavailable direct refits",
            "",
            f"These **{len(unavailable_papers)} papers** contain "
            f"{unavailable_record_count} records for which a",
            "source-faithful coefficient refit could not be performed. A paper can also",
            "have other records that were reproduced.",
            "",
            "See [incomplete datasets and author requests](dataset-requests.md) for",
            "reviewed outreach recommendations, exact requests, and contact status.",
            "Papers not yet listed there remain unassessed for outreach.",
            "",
            "| Paper | Affected records | Why direct refitting was unavailable |",
            "|---|---|---|",
        ]
    )
    for paper in unavailable_papers:
        unavailable = [
            record
            for record in paper["records"]
            if record["status"] == "not_refittable"
        ]
        reasons = []
        for record in unavailable:
            reason = record["reason"]
            if reason not in reasons:
                reasons.append(reason)
        lines.append(
            f"| {source_link(paper)} | {record_list(unavailable)} | "
            f"{markdown_cell(' '.join(reasons))} |"
        )

    lines.extend(
        [
            "",
            "## Complete investigated-paper register",
            "",
            "This is the exhaustive paper-level index. `Bundled` means numerical primary",
            "rows are stored; `bundled indirect` means source coefficient or derived-output",
            "tables are stored but do not form a common direct-fit observation matrix;",
            "`plot only` means observations were digitized;",
            "`parameterization only` means only the published equation/coefficients can",
            "be checked; `partial reconstruction` means a documented subset of a",
            "larger source objective is executable. Record-level links, fit metrics, and evidence locations are in",
            "the primary-source and refit ledgers.",
            "",
            "| Paper | Final disposition | Catalog records | Record-level results | Primary-data form |",
            "|---|---|---:|---|---|",
        ]
    )
    data_labels = {
        "digitized": "digitized",
        "partial_published_table": "partial published table/checkpoints",
        "bundled": "bundled",
        "bundled_indirect": "bundled indirect",
        "combined_primary_data_partially_bundled": "partial primary rows bundled",
        "plot_only": "plot only/digitized",
        "parameterization_only": "parameterization only",
        "theoretical_parameterization_only": "theoretical parameterization only",
        "partial_reconstruction": "partial reconstruction",
        "external_primary_table_refitted": (
            "external primary table (refitted; not redistributed)"
        ),
        "externally_recovered_not_redistributed": (
            "externally recovered rows (refitted; not redistributed)"
        ),
    }
    for paper in papers:
        if paper["records"]:
            result = paper.get("reproduction_summary") or compact_counts(
                paper["statuses"], STATUS_LABELS
            )
            scope_note = PAPER_SCOPE_NOTES.get(str(paper.get("doi", "")).lower())
            if scope_note:
                result += "; " + scope_note
            data_form = compact_counts(paper["data_statuses"], data_labels)
            qualified_final_input = sum(
                "not a reconstruction of every upstream reduction"
                in str(record.get("qualification", ""))
                for record in paper["records"]
            )
            if qualified_final_input:
                data_form += f"; {qualified_final_input} final-input parity, upstream reduction partial"
            conditional_fixed_volume = sum(
                record.get("parity_basis") == "conditional_fixed_V0_diagnostic_only"
                for record in paper["records"]
            )
            if conditional_fixed_volume:
                data_form += (
                    f"; {conditional_fixed_volume} conditional fixed-V0 refits; "
                    "source weights and marker averaging unspecified"
                )
            record_count = len(paper["records"])
        else:
            result = "no production record"
            data_form = "investigation evidence only"
            record_count = 0
        lines.append(
            f"| {source_link(paper)} | {OUTCOME_LABELS[paper['outcome']]} | "
            f"{record_count} | {result} | {data_form} |"
        )

    incomplete_alternatives = [
        paper for paper in catalog_papers if paper.get("withheld_candidates")
    ]
    if incomplete_alternatives:
        lines.extend(["", "## Withheld alternatives in catalog papers", ""])
        for paper in incomplete_alternatives:
            for candidate in paper["withheld_candidates"]:
                closure = ""
                evidence = candidate.get(
                    "pressure_convention_audit", candidate["evidence"]
                )
                if (
                    candidate.get("investigation_status")
                    == "closed_pending_author_inputs"
                ):
                    closure = (
                        f"**Closed for now ({candidate['closure_date']}); "
                        "awaiting author inputs.** "
                    )
                lines.append(
                    f"- {source_link(paper)}, {candidate['candidate']}: "
                    f"{closure}{candidate['reason']} "
                    f"[Audit]({evidence})."
                )

    lines.extend(
        [
            "",
            "## Maintenance",
            "",
            "This page is generated by `scripts/generate_paper_investigation_ledger.py`.",
            "Update record-level evidence first, add nonproduction investigations to",
            "`docs/data/nonproduction-paper-investigations.json`, regenerate this page,",
            "and run the generator with `--check` in validation.",
            "",
        ]
    )
    return "\n".join(lines)


def main() -> None:
    """Write the ledger or check that the committed copy is current."""
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true", help="fail if output is stale")
    args = parser.parse_args()
    rendered = render()
    if args.check:
        current = (
            OUTPUT_PATH.read_text(encoding="utf-8") if OUTPUT_PATH.exists() else ""
        )
        if current != rendered:
            raise SystemExit(f"stale generated file: {OUTPUT_PATH.relative_to(ROOT)}")
        return
    OUTPUT_PATH.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
