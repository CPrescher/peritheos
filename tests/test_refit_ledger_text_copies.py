"""The refit ledger's copies of record text come from the current record."""

from __future__ import annotations

import json
from pathlib import Path

from peritheos import get_material_document
from scripts.review_similar_parameter_scope import TEXT_COPIES, saved_review
from scripts.validate_primary_eos_refits import validate_all

ROOT = Path(__file__).resolve().parents[1]
MATERIAL = "coesite"
RECORD = "coesite_levien_1981_bm3_1"


def test_a_reviewed_snapshot_does_not_override_current_record_text():
    # The record is in the reviewed cohort, whose snapshot keeps older text.
    assert saved_review()[RECORD]["status"] == "similar"
    record = next(
        r
        for r in get_material_document(MATERIAL)["eos_records"]
        if r["identifier"] == RECORD
    )
    row = next(
        r
        for r in validate_all([MATERIAL])["records"]
        if r["record_identifier"] == RECORD
    )
    assert row["label"] == record["label"]
    assert row["source_notes"] == record.get("notes")
    saved = next(
        r
        for r in json.loads((ROOT / "docs/data/primary-eos-refits.json").read_bytes())[
            "records"
        ]
        if r["record_identifier"] == RECORD
    )
    for key in TEXT_COPIES:
        assert row.get(key) == saved.get(key)
