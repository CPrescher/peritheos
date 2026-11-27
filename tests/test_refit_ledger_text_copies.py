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


def test_ledger_check_tolerates_roundoff_except_in_source_facts():
    import pytest

    from scripts.validate_primary_eos_refits import check_ledger_values

    saved = {
        "reproduction": {"refits": {"parameters": [211.74492423778685]}},
        "parameters": [{"published": 184.5, "refit": 184.5}],
        "observations": 19,
    }

    def changed(**edits):
        row = json.loads(json.dumps(saved))
        row["reproduction"]["refits"]["parameters"][0] = edits.get(
            "nested", 211.74492423778685
        )
        row["parameters"][0]["published"] = edits.get("published", 184.5)
        row["observations"] = edits.get("observations", 19)
        return row

    # A calculated number under any other key may differ by solver roundoff.
    check_ledger_values(saved, changed(nested=211.74500162186482))
    with pytest.raises(ValueError, match="calculated value differs"):
        check_ledger_values(saved, changed(nested=211.76))
    # Source facts and counts stay exact.
    with pytest.raises(ValueError, match="source/assessment value differs"):
        check_ledger_values(saved, changed(published=184.5000001))
    with pytest.raises(ValueError, match="source/assessment value differs"):
        check_ledger_values(saved, changed(observations=20))
