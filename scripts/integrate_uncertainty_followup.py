"""Integrate only the 63 reviewed conditional-uncertainty assessments.

Repository maintenance only. Scientific source fields and unrelated assessments
are preserved. Run after reconstruct_similar_eos_uncertainties and its checks.
"""

from __future__ import annotations

import json
from collections import Counter

from peritheos import convert_eosmat, load_eosmat, validate_eosmat_document
from peritheos._eosmat_format import format_eosmat
from scripts.reconstruct_similar_eos_uncertainties import (
    COMBINED_ERROR_POLICY,
    PARITY_POLICY,
    ROOT,
    augment_outcome,
    saved_followup,
)
from scripts.standardize_reconstruction_metadata import assessment_for_record
from scripts.validate_primary_eos_refits import render_markdown


def integrate():
    ledger_path = ROOT / "docs/data/primary-eos-refits.json"
    ledger = json.loads(ledger_path.read_text())
    by_material = {}
    for row in ledger["records"]:
        if row["record_identifier"] not in saved_followup():
            continue
        material = row["material_identifier"]
        if material not in by_material:
            by_material[material] = load_eosmat(
                ROOT / "peritheos/data/materials" / (material + ".eosmat")
            )
        document = by_material[material]
        native = next(
            r
            for r in document["eos_records"]
            if r["identifier"] == row["record_identifier"]
        )
        legacy = next(
            r
            for r in convert_eosmat(document, version=3)["eos_records"]
            if r["identifier"] == row["record_identifier"]
        )
        row.update(augment_outcome(legacy, row))
        native["evidence"]["reproduction"]["assessment"] = assessment_for_record(
            native, row
        )
    counts = Counter(r["status"] for r in ledger["records"])
    ledger["summary"] = {
        "total": len(ledger["records"]),
        **dict(sorted(counts.items())),
    }
    ledger["policy"]["parity"] = ledger["policy"]["parity"].replace(
        "Digitization/run-selection hypotheses, unresolved source-data disagreement and diagnostic source masks can retain similar despite coefficient agreement.",
        "At parameter_agreement scope,"
        + PARITY_POLICY.split("At parameter_agreement scope,", 1)[1],
    )
    if PARITY_POLICY not in ledger["policy"]["parity"]:
        ledger["policy"]["parity"] += " " + PARITY_POLICY
    if COMBINED_ERROR_POLICY not in ledger["policy"]["parity"]:
        ledger["policy"]["parity"] += " " + COMBINED_ERROR_POLICY
    for material, document in by_material.items():
        validate_eosmat_document(document)
        path = ROOT / "peritheos/data/materials" / (material + ".eosmat")
        original, formatted = path.read_text(), format_eosmat(document)
        for r in document["eos_records"]:
            if r["identifier"] not in saved_followup():
                continue
            marker = '      "identifier": "' + r["identifier"] + '",'

            def bounds(text):
                at = text.index(marker)
                return text.rfind("    {", 0, at), text.index("\n    }", at) + len(
                    "\n    }"
                )

            a, b = bounds(original)
            c, d = bounds(formatted)
            original = original[:a] + formatted[c:d] + original[b:]
        path.write_text(original)
    ledger_path.write_text(
        json.dumps(ledger, indent=2, sort_keys=True, allow_nan=False) + "\n"
    )
    (ROOT / "docs/primary-eos-refits.md").write_text(render_markdown(ledger))
    from scripts.report_numbers import value_error

    report = ROOT / "docs/literature-reproductions/uncertainty-reconstruction.md"
    text = report.read_text()
    table = [
        "| Record | Outcome | Reconstructed free coefficients ± conditional standard error |",
        "| --- | --- | --- |",
    ]
    for row in saved_followup().values():
        coefficients = []
        for p in row["parameters"]:
            value, error = value_error(p["refit"], p["refit_error"])
            coefficients.append(f"{p['parameter']}: {value} ± {error}")
        rid = row["record_identifier"]
        link = f"../reconstruction-reports/{row['material_identifier']}.md#record-{rid}"
        table.append(
            f"| [{rid}]({link}) | {row['status']} | {'; '.join(coefficients)} |"
        )
    start, end = (
        "<!-- uncertainty-followup-table -->",
        "<!-- /uncertainty-followup-table -->",
    )
    a, b = text.index(start) + len(start), text.index(end)
    report.write_text(text[:a] + "\n\n" + "\n".join(table) + "\n\n" + text[b:])
    print(json.dumps(ledger["summary"], sort_keys=True))


if __name__ == "__main__":
    integrate()
