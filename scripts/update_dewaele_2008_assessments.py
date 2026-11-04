"""Refresh this paper's 13 catalog and 54 research assessments separately.

Repository maintenance only. Published coefficients and unrelated papers stay
unchanged. PAW records never enter the catalog ledgers through this command.
"""

from __future__ import annotations

import json
from collections import Counter

from peritheos._eosmat_format import format_eosmat
from peritheos.eosmat import validate_eosmat_document
from scripts.move_dewaele_2008_paw_to_research import (
    COLLECTION,
    move,
    relocate_record_links,
    write_json,
)
from scripts.register_dewaele_2008_metals import audit
from scripts.reproduce_dewaele_2008_metals import (
    PARITY_POLICY,
    REPORT,
    ROOT,
    Z,
    reproduce,
)
from scripts.standardize_reconstruction_metadata import assessment_for_record
from scripts.validate_primary_eos_refits import render_markdown


def update():
    report = reproduce()
    outcomes = report["records"]
    updated, documents = {}, {}
    for material in Z:
        for directory in (ROOT / "peritheos/data/materials", COLLECTION):
            path = directory / f"{material}.eosmat"
            if not path.exists():
                continue
            original = path.read_bytes()
            document = json.loads(original)
            for record in document["eos_records"]:
                rid = record["identifier"]
                if rid not in outcomes:
                    continue
                existing = record["evidence"]["scientific_validation"]
                previous_report = record["evidence"]["reproduction"].get(
                    "previous_catalog_report"
                )
                audit(
                    record,
                    outcomes[rid],
                    existing["primary_source_check"]["locations"][-1],
                    existing["status"],
                )
                if previous_report:
                    record["evidence"]["reproduction"]["previous_catalog_report"] = (
                        previous_report
                    )
                if rid in updated:
                    raise ValueError(
                        f"Record appears in catalog and research collection: {rid}"
                    )
                updated[rid] = record
            documents[path] = (original, document)
    assert set(updated) == set(outcomes)
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")

    canonical_rows = {}
    for path in (
        ROOT / "docs/data/primary-eos-refits.json",
        COLLECTION / "reproduction.json",
    ):
        ledger = json.loads(path.read_text())
        for row in ledger["records"]:
            rid = row["record_identifier"]
            if rid not in outcomes:
                continue
            outcome = outcomes[rid]
            validation = updated[rid]["evidence"]["scientific_validation"]
            row.update(outcome)
            row["primary_data_status"] = validation["primary_data_check"]["status"]
            row["primary_data_finding"] = validation["primary_data_check"]["finding"]
            if "qualification" not in outcome:
                row.pop("qualification", None)
            for field in (
                "usage_recommendation",
                "recommendation_reason",
                "recommendation_date",
                "recommendation_audit",
            ):
                row.pop(field, None)
                if field in validation:
                    row[field] = validation[field]
            if validation.get("recommendation_reason"):
                row["qualification"] = (
                    outcome.get("qualification", "")
                    + " "
                    + validation["recommendation_reason"]
                ).strip()
        counts = dict(sorted(Counter(r["status"] for r in ledger["records"]).items()))
        ledger["summary"] = {"total": len(ledger["records"]), **counts}
        if PARITY_POLICY not in ledger["policy"]["parity"]:
            ledger["policy"]["parity"] += " " + PARITY_POLICY
        path.write_text(
            json.dumps(ledger, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
        canonical_rows.update(
            {r["record_identifier"]: r for r in json.loads(path.read_text())["records"]}
        )
        if path.parent != COLLECTION:
            (ROOT / "docs/primary-eos-refits.md").write_text(render_markdown(ledger))
    for path, (original, document) in documents.items():
        for record in document["eos_records"]:
            if record["identifier"] not in updated:
                continue
            record["evidence"]["reproduction"]["assessment"] = assessment_for_record(
                record, canonical_rows[record["identifier"]]
            )
            if path.parent == COLLECTION:
                relocate_record_links(record)
        validate_eosmat_document(document)
        rendered = format_eosmat(document).encode()
        if rendered != original:
            if path.read_bytes() != original:
                raise RuntimeError(f"Material changed during update: {path}")
            path.write_bytes(rendered)

    for path in (
        ROOT / "peritheos/data/primary-source-audit.json",
        COLLECTION / "source-audit.json",
    ):
        source_audit = json.loads(path.read_text())
        for row in source_audit["records"]:
            if row["record"] not in updated:
                continue
            record = updated[row["record"]]
            validation = record["evidence"]["scientific_validation"]
            row.update(
                label=record["label"],
                status=validation["status"],
                note=validation["note"],
                primary_source_check=validation["primary_source_check"],
            )
            for field in (
                "usage_recommendation",
                "reproduction_status",
                "recommendation_date",
                "recommendation_audit",
                "recommendation_reason",
            ):
                if field in validation:
                    row[field] = validation[field]
                else:
                    row.pop(field, None)
        total_key = "records"
        source_audit["summary"] = {
            total_key: len(source_audit["records"]),
            **dict(
                sorted(Counter(r["status"] for r in source_audit["records"]).items())
            ),
        }
        write_json(path, source_audit)
    move()  # Refresh research documentation and checksums, without re-registering it.
    print("Updated 13 catalog and 54 research EOS assessments; no coefficient changes.")


if __name__ == "__main__":
    update()
