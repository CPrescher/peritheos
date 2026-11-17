"""Curated ``check`` texts: the per-record source note and the warning reason.

``scientific_validation.primary_source_check.record_note`` replaces the verbatim
``finding`` in ``check.source_note``; ``scientific_validation.warning_reason``
is copied to ``check.warning_reason``. The Rust mirror is tested in
``crates/peritheos/tests/eosmat_audit.rs`` with the same inputs.
"""

import copy
import json

import jsonschema
import pytest

from peritheos import Material
from peritheos._eosmat_audit import (
    RECORD_NOTE_MAX_LENGTH,
    WARNING_REASON_MAX_LENGTH,
    check_block,
    derived_report_links,
    merge_audit,
    split_document,
)
from peritheos._eosmat_format import format_eosmat
from peritheos.errors import EosmatError

NOTE = "Table 2 V0 and K0 match; V0 converted from cm3/mol."
REASON = "The fit omits the two highest-pressure points."
WARNING = "Qualified use only; see the audit."


def curated(note=None, reason=None, *, restricted=True):
    """The Rust test's ``curated`` document, plus its report links."""
    validation = {
        "status": "primary_source_validated",
        "audit_date": "2026-01-02",
        "primary_source_check": {"finding": "Table 2 values match."},
    }
    if restricted:
        validation["usage_recommendation"] = "qualified_use"
    if note is not None:
        validation["primary_source_check"]["record_note"] = note
    if reason is not None:
        validation["warning_reason"] = reason
    return {
        "format_version": 5,
        "identifier": "example",
        "formula": "MgO",
        "eos_records": [
            {
                "identifier": "example_bm3",
                "model": "BM3",
                "temperature_ref": 300,
                "parameters": {"V0": 74.7, "K0": 160, "K0_prime": 4},
                "reference": "Fictional example; no experimental validation.",
                "evidence": {
                    "scientific_validation": validation,
                    "reproduction": {
                        "assessment": {
                            "status": "parity",
                            "scope": "parameter_agreement",
                            "basis": "Independent refit of Table 1.",
                            # Present so the document can be split.
                            **derived_report_links("example", "example_bm3"),
                        }
                    },
                },
            }
        ],
    }


def record_of(document):
    return Material.from_eosmat(document, require_primary_validation=False).eos_records[
        0
    ]


def test_both_fields_are_copied_into_check_and_read_back():
    document = curated(NOTE, REASON)
    derived = check_block(document["eos_records"][0], document)
    # The same stored values the Rust test reads back.
    assert derived == {
        "reproduction": "reproduced",
        "compared": "published_fit",
        "warning": WARNING,
        "warning_reason": REASON,
        "source_check": "verified",
        "source_checked_on": "2026-01-02",
        "source_note": NOTE,
    }
    record = record_of(document)
    assert (record.source_note, record.warning, record.warning_reason) == (
        NOTE,
        WARNING,
        REASON,
    )
    # The finding itself is untouched.
    audit = dict(record.audit)["evidence"]["scientific_validation"]
    assert audit["primary_source_check"]["finding"] == "Table 2 values match."
    # A stored check block carrying both reads back unchanged.
    stored = copy.deepcopy(document)
    stored["eos_records"][0]["check"] = derived
    assert record_of(stored).warning_reason == REASON
    # warning_reason is formatted directly after warning.
    check = json.loads(format_eosmat(stored))["eos_records"][0]["check"]
    assert list(check)[:5] == [
        "reproduction",
        "compared",
        "warning",
        "warning_reason",
        "source_check",
    ]


def test_absent_fields_keep_the_verbatim_finding_and_no_reason():
    document = curated()
    derived = check_block(document["eos_records"][0], document)
    assert derived["source_note"] == "Table 2 values match."
    assert "warning_reason" not in derived
    record = record_of(document)
    assert record.source_note == "Table 2 values match."
    assert record.warning == WARNING
    assert record.warning_reason is None


def test_split_and_merge_keep_both_fields_in_the_sidecar():
    document = curated(NOTE, REASON)
    records, sidecar = split_document(document)
    check = records["eos_records"][0]["check"]
    assert (check["source_note"], check["warning_reason"]) == (NOTE, REASON)
    validation = sidecar["records"]["example_bm3"]["evidence"]["scientific_validation"]
    assert validation["primary_source_check"]["record_note"] == NOTE
    assert validation["warning_reason"] == REASON
    # The record file alone reads the stored check; merged, the evidence agrees.
    alone = record_of(records)
    assert (alone.source_note, alone.warning_reason) == (NOTE, REASON)
    merged = record_of(merge_audit(records, sidecar))
    assert (merged.source_note, merged.warning_reason) == (NOTE, REASON)


@pytest.mark.parametrize(
    ("note", "reason", "restricted", "message"),
    [
        (
            "x" * (RECORD_NOTE_MAX_LENGTH + 1),
            None,
            True,
            "primary_source_check.record_note must be at most 300 characters",
        ),
        (
            None,
            "x" * (WARNING_REASON_MAX_LENGTH + 1),
            True,
            "warning_reason must be at most 140 characters",
        ),
        (" ", None, True, "record_note must be a non-empty string"),
        (None, "", True, "warning_reason must be a non-empty string"),
        (None, 3, True, "warning_reason must be a non-empty string"),
        (
            None,
            "Fictional reason.",
            False,
            "warning_reason is given but the record has no warning",
        ),
    ],
)
def test_invalid_sidecar_texts_are_rejected(note, reason, restricted, message):
    document = curated(note, reason, restricted=restricted)
    with pytest.raises(EosmatError, match=message):
        record_of(document)
    # Split into a record file and a sidecar, the merged document is rejected
    # the same way; the reason is never silently dropped from check.
    if restricted or reason is None:
        return
    records, sidecar = split_document(document)
    assert records["eos_records"][0]["check"]["warning_reason"] == reason
    with pytest.raises(EosmatError, match="only given with a warning"):
        record_of(records)


def test_limits_are_inclusive():
    document = curated("x" * RECORD_NOTE_MAX_LENGTH, "x" * WARNING_REASON_MAX_LENGTH)
    record = record_of(document)
    assert len(record.source_note) == RECORD_NOTE_MAX_LENGTH
    assert len(record.warning_reason) == WARNING_REASON_MAX_LENGTH


@pytest.mark.parametrize(
    ("check", "message"),
    [
        ({"warning": None, "warning_reason": "R."}, "only given with a warning"),
        ({"warning_reason": "R."}, "only given with a warning"),
        (
            {"warning": "Do not use.", "warning_reason": "x" * 141},
            "at most 140 characters",
        ),
        ({"warning": "Do not use.", "warning_reason": ""}, "non-empty string"),
    ],
)
def test_check_warning_reason_needs_a_warning_and_obeys_the_limit(check, message):
    document = curated()
    del document["eos_records"][0]["evidence"]
    document["eos_records"][0]["check"] = check
    with pytest.raises(EosmatError, match=message):
        record_of(document)


def test_schema_describes_both_fields():
    from importlib import resources

    schema = json.loads(
        resources.files("peritheos.data")
        .joinpath("eosmat-v5.schema.json")
        .read_text(encoding="utf-8")
    )
    validator = jsonschema.Draft202012Validator(schema)
    records, sidecar = split_document(curated(NOTE, REASON))
    validator.validate(records)
    validator.validate(merge_audit(records, sidecar))
    too_long = curated(NOTE, REASON)
    too_long["eos_records"][0]["evidence"]["scientific_validation"][
        "warning_reason"
    ] = "x" * 141
    assert not validator.is_valid(too_long)
    records["eos_records"][0]["check"]["warning_reason"] = "x" * 141
    assert not validator.is_valid(records)


def test_overview_and_summary_show_the_reason():
    record = record_of(curated(NOTE, REASON))
    check = record.overview()["check"]
    assert (check["warning"], check["warning_reason"], check["source_note"]) == (
        WARNING,
        REASON,
        NOTE,
    )
    lines = record.summary().splitlines()
    assert f"  Warning: {WARNING} {REASON}" in lines
    assert any(line.endswith(f". {NOTE}") for line in lines)
    plain = record_of(curated())
    assert plain.overview()["check"]["warning_reason"] is None
    assert f"  Warning: {WARNING}" in plain.summary().splitlines()


def test_stored_copies_must_equal_the_inline_evidence():
    """Mirrors ``stored_check_texts_must_equal_the_evidence`` in Rust."""
    for note, reason in ((NOTE, REASON), (None, None)):
        document = curated(note, reason)
        record = document["eos_records"][0]
        record["check"] = check_block(record, document)
        record_of(document)
    cases = [
        # A stored note must be the record_note, else the verbatim finding.
        (curated(NOTE, REASON), "source_note", "Table 2 values match.", "source_note"),
        (curated(None, REASON), "source_note", NOTE, "source_note"),
        # A stored reason must be the evidence's reason, also when it has none.
        (curated(NOTE, REASON), "warning_reason", "Other reason.", "warning_reason"),
        (curated(NOTE, None), "warning_reason", REASON, "warning_reason"),
    ]
    for document, key, value, field in cases:
        record = document["eos_records"][0]
        record["check"] = check_block(record, document)
        record["check"][key] = value
        with pytest.raises(
            EosmatError, match=f"check.{field} differs from the evidence"
        ):
            record_of(document)
    # A standalone file without evidence keeps its own texts.
    document = curated()
    del document["eos_records"][0]["evidence"]
    document["eos_records"][0]["check"] = {
        "warning": "Do not use.",
        "warning_reason": "My own reason.",
        "source_note": "My own note.",
    }
    assert record_of(document).source_note == "My own note."


@pytest.mark.parametrize("blank", [" ", "\t\n", " ", " 　"])
def test_white_space_only_texts_are_rejected_in_check_and_evidence(blank):
    document = curated()
    del document["eos_records"][0]["evidence"]
    for check in (
        {"warning": "Do not use.", "warning_reason": blank},
        {"source_note": blank},
    ):
        document["eos_records"][0]["check"] = check
        with pytest.raises(EosmatError, match="must be a non-empty string"):
            record_of(document)
    for note, reason in ((blank, None), (None, blank)):
        with pytest.raises(EosmatError, match="must be a non-empty string"):
            record_of(curated(note, reason))


def test_only_unicode_white_space_is_blank():
    # U+001F is white space to str.isspace but not Unicode White_Space; like
    # the Rust reader, the library keeps such a text.
    assert "\x1f".isspace()
    record = record_of(curated("\x1f", None))
    assert record.source_note == "\x1f"
