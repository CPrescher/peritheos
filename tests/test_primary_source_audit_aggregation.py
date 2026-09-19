"""Regenerating the catalog must preserve ledger-only source evidence."""

import json

from scripts import apply_primary_source_audit as audit


def test_aggregate_preserves_candidate_evidence_and_missing_card_note(
    tmp_path, monkeypatch
):
    materials = tmp_path / "materials"
    materials.mkdir()
    for name in ("qandilite.eosmat", "vanadium_bcc.eosmat", "manifest.json"):
        (materials / name).write_bytes((audit.MATERIALS / name).read_bytes())
    report = tmp_path / "audit.json"
    source_report = json.loads(audit.REPORT.read_text())
    selected = [
        row
        for row in source_report["records"]
        if row["material"] in {"qandilite", "vanadium_bcc"}
    ]
    report.write_text(json.dumps({"records": selected}))
    monkeypatch.setattr(audit, "MATERIALS", materials)
    monkeypatch.setattr(audit, "REPORT", report)
    audit.aggregate_existing_audits()
    regenerated = json.loads(report.read_text())
    by_id = {row["record"]: row for row in regenerated["records"]}
    for previous in selected:
        assert by_id[previous["record"]] == previous
    assert regenerated["summary"]["records"] == 3
    assert (
        len(by_id["vanadium_bcc_crichton_2016_bm3_thermal"]["withheld_candidates"]) == 2
    )
    first = report.read_bytes()
    audit.aggregate_existing_audits()
    assert report.read_bytes() == first
