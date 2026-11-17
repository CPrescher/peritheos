"""Repository check: bundled record text within the record text standard.

Reports every bundled record whose ``label``, ``notes``, ``validity.notes``,
``thermal.notes``, ``check.source_note`` or ``check.warning_reason``, and every
material whose ``notes``, exceed the limits of ``docs/record-text-standard.md``,
and every label that is not unique within its material. ``check.source_note``
is the record's curated ``record_note`` when it has one, so writing a note
within the limit (``scripts/apply_check_text.py``) resolves a
``source_note_length`` case. Known cases are listed in
``scripts/record-text-baseline.json``; the test suite fails on a case that is
not listed and on a listed case that no longer occurs, so the baseline only
shrinks. ``--prune`` removes cases that no longer occur; a new case is fixed in
the text (``scripts/apply_record_text.py``), not added to the baseline.

Usage::

    python scripts/check_record_text.py           # report; exit 1 on a change
    python scripts/check_record_text.py --prune   # drop resolved baseline cases

This is not an installed-package tool.
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from peritheos._eosmat_audit import (  # noqa: E402
    RECORD_NOTE_MAX_LENGTH,
    WARNING_REASON_MAX_LENGTH,
)
from scripts.apply_record_text import (  # noqa: E402
    LABEL_MAX_CHARACTERS,
    NOTES_MAX_CHARACTERS,
    NOTES_MAX_SENTENCES,
    THERMAL_NOTES_MAX_CHARACTERS,
    VALIDITY_MAX_NOTES,
    VALIDITY_NOTE_MAX_CHARACTERS,
    sentence_count,
)

BASELINE = ROOT / "scripts/record-text-baseline.json"
MATERIAL_NOTES_MAX_CHARACTERS = 300
SOURCE_NOTE_MAX_CHARACTERS = RECORD_NOTE_MAX_LENGTH
WARNING_REASON_MAX_CHARACTERS = WARNING_REASON_MAX_LENGTH


def record_cases(record: dict) -> list[str]:
    """Names of the limits one record exceeds."""
    cases = []
    label = record.get("label")
    if isinstance(label, str) and len(label) > LABEL_MAX_CHARACTERS:
        cases.append("label_length")
    notes = record.get("notes")
    if isinstance(notes, str) and (
        len(notes) > NOTES_MAX_CHARACTERS or sentence_count(notes) > NOTES_MAX_SENTENCES
    ):
        cases.append("notes_length")
    validity = record.get("validity")
    validity_notes = validity.get("notes") if isinstance(validity, dict) else None
    if isinstance(validity_notes, list) and (
        len(validity_notes) > VALIDITY_MAX_NOTES
        or any(len(note) > VALIDITY_NOTE_MAX_CHARACTERS for note in validity_notes)
    ):
        cases.append("validity_length")
    thermal = record.get("thermal")
    thermal_notes = thermal.get("notes") if isinstance(thermal, dict) else None
    if (
        isinstance(thermal_notes, str)
        and len(thermal_notes) > THERMAL_NOTES_MAX_CHARACTERS
    ):
        cases.append("thermal_length")
    check = record.get("check")
    source_note = check.get("source_note") if isinstance(check, dict) else None
    if isinstance(source_note, str) and len(source_note) > SOURCE_NOTE_MAX_CHARACTERS:
        cases.append("source_note_length")
    reason = check.get("warning_reason") if isinstance(check, dict) else None
    if isinstance(reason, str) and len(reason) > WARNING_REASON_MAX_CHARACTERS:
        cases.append("warning_reason_length")
    return cases


def current_cases(root: Path = ROOT) -> dict[str, list[str]]:
    """Every case in the bundled catalog, keyed ``material`` or ``material/record``."""
    cases: dict[str, list[str]] = {}
    for path in sorted((root / "peritheos/data/materials").glob("*.eosmat")):
        document = json.loads(path.read_bytes())
        material = path.stem
        notes = document.get("notes")
        if isinstance(notes, str) and len(notes) > MATERIAL_NOTES_MAX_CHARACTERS:
            cases[material] = ["material_notes_length"]
        labels = [r.get("label") for r in document.get("eos_records", [])]
        for record in document.get("eos_records", []):
            found = record_cases(record)
            if record.get("label") is not None and labels.count(record["label"]) > 1:
                found.append("label_unique")
            if found:
                cases[f"{material}/{record['identifier']}"] = found
    return cases


def compare(
    current: dict[str, list[str]], baseline: dict[str, list[str]]
) -> tuple[list[str], list[str]]:
    """New cases (not in the baseline) and resolved ones (no longer present)."""
    new = [
        f"{key}: {case}"
        for key, found in current.items()
        for case in found
        if case not in baseline.get(key, [])
    ]
    resolved = [
        f"{key}: {case}"
        for key, listed in baseline.items()
        for case in listed
        if case not in current.get(key, [])
    ]
    return new, resolved


def load_baseline(path: Path = BASELINE) -> dict[str, list[str]]:
    return json.loads(path.read_bytes())["cases"]


def write_baseline(cases: dict[str, list[str]], path: Path = BASELINE) -> None:
    document = {
        "description": (
            "Bundled records and materials whose text exceeded the limits of "
            "docs/record-text-standard.md when the check was introduced. "
            "scripts/check_record_text.py fails on any case not listed here "
            "and on a listed case that no longer occurs (remove it with "
            "--prune). Never add a case; fix the text instead."
        ),
        "cases": dict(sorted(cases.items())),
    }
    path.write_text(json.dumps(document, indent=2, ensure_ascii=False) + "\n")


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument(
        "--prune", action="store_true", help="drop resolved cases from the baseline"
    )
    arguments = parser.parse_args(argv)
    current = current_cases()
    baseline = load_baseline()
    new, resolved = compare(current, baseline)
    if arguments.prune:
        kept = {
            key: [case for case in listed if case in current.get(key, [])]
            for key, listed in baseline.items()
        }
        write_baseline({key: listed for key, listed in kept.items() if listed})
        print(f"Removed {len(resolved)} resolved cases from the baseline.")
        resolved = []
    for line in new:
        print(f"NEW {line}")
    for line in resolved:
        print(f"RESOLVED {line} (run --prune)")
    total = sum(map(len, current.values()))
    print(f"{total} cases in the catalog; {len(new)} new; {len(resolved)} resolved.")
    return 1 if new or resolved else 0


if __name__ == "__main__":
    sys.exit(main())
