"""Repository maintenance: apply a reviewed rewrite of bundled record text.

Free text in bundled records (``label``, ``notes``, ``validity.notes`` and
``thermal.notes``) follows ``docs/record-text-standard.md``. This tool writes a
reviewed rewrite into ``peritheos/data/materials/<id>.eosmat`` and keeps the
complete original of every replaced field verbatim in the record's
``superseded_record_text`` block of ``peritheos/data/audits/<id>.audit.json``,
so no wording is lost. A later rewrite of the same field keeps the first
original. Nothing else in either file changes; both are written through the
repository formatter, and the audit-split receipt binds the revised files
(``scripts/split_eosmat_audits.py --verify`` still recovers every pre-split
file byte for byte). The registered refit ledger
(``docs/data/primary-eos-refits.json``) copies each record's label and notes;
those two copies are updated and ``docs/primary-eos-refits.md`` is re-rendered
from the saved ledger, without recomputing anything.

Input files (JSON), any number::

    {"<material>": {"<record>": {
        "before": {"notes": str|null, "validity_notes": [str]|null,
                   "thermal_notes": str|null},
        "after": {same keys},
        ...other keys are ignored}}}

``null`` (or an empty list) in ``after`` removes the field; a ``validity``
object left empty is removed. Labels come only from ``--labels`` (CSV columns
``material, record, current_label, proposed_label, flag``); only rows flagged
``mechanical`` are applied.

A record is refused, and reported, when a ``before`` text or the CSV's current
label differs from the file, or when one of these guards fails:

``numbers``
    an after text contains a number (decimals, or integers of two or more
    digits) found neither in the record's current text fields nor in its
    structured fields;
``notes_length``
    notes over 400 characters or 3 sentences;
``validity_length``
    a validity note over 160 characters, or more than 2 notes;
``thermal_length``
    thermal notes over 250 characters;
``label_length``
    a label over 80 characters;
``label_unique``
    two records of one material would share a label;
``forbidden_text``
    a repository path (``docs/``, ``.json``, ``scripts/``) or the words
    ``LitCurate`` or ``Audit 20`` in an after text.

``--allow RECORD:GUARD`` overrides one guard for one named record; every
override is listed in the output. A refused record is left unchanged; the
others are applied and the exit status is 1.

Usage::

    python scripts/apply_record_text.py BATCH.json ... [--labels CSV] [--check]
    python scripts/apply_record_text.py --verify [--revision c26e8d81]

``--check`` reports what would change and every refusal without writing.
``--verify`` confirms, for every record holding a ``superseded_record_text``
block, that restoring the stored originals gives the record's text at the
revision, that no other record's text differs from that revision, and that
the ledger's copies equal the catalog. This is not an installed-package tool.
"""

from __future__ import annotations

import argparse
import copy
import csv
import datetime
import json
import re
import subprocess
import sys
from collections.abc import Callable, Iterable, Mapping
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from peritheos._eosmat_audit import (  # noqa: E402
    SUPERSEDED_TEXT_FIELDS,
    merge_audit,
    validate_superseded_record_text,
)
from peritheos._eosmat_format import format_audit, format_eosmat  # noqa: E402
from peritheos.eosmat import validate_eosmat_document  # noqa: E402
from scripts.split_eosmat_audits import (  # noqa: E402
    RECEIPT,
    bind_record_text_revision,
    record_text,
    restore_record_text,
    set_record_text,
)

DEFAULT_REVISION = "c26e8d81"
#: The registered refit ledger copies each record's label and notes.
LEDGER = "docs/data/primary-eos-refits.json"
LEDGER_MARKDOWN = "docs/primary-eos-refits.md"
#: Fields an input file may rewrite; labels come only from the label CSV.
JSON_FIELDS = ("notes", "validity_notes", "thermal_notes")
NOTES_MAX_CHARACTERS = 400
NOTES_MAX_SENTENCES = 3
VALIDITY_NOTE_MAX_CHARACTERS = 160
VALIDITY_MAX_NOTES = 2
THERMAL_NOTES_MAX_CHARACTERS = 250
LABEL_MAX_CHARACTERS = 80
FORBIDDEN_TEXT = ("docs/", ".json", "scripts/", "LitCurate", "Audit 20")
GUARDS = (
    "numbers",
    "notes_length",
    "validity_length",
    "thermal_length",
    "label_length",
    "label_unique",
    "forbidden_text",
)
#: Abbreviations whose full stop does not end a sentence.
ABBREVIATIONS = (
    "et al.",
    "e.g.",
    "i.e.",
    "cf.",
    "vs.",
    "approx.",
    "ca.",
    "Fig.",
    "Figs.",
    "Eq.",
    "Eqs.",
    "Ref.",
    "Refs.",
    "No.",
    "Tab.",
    "Suppl.",
)
_NUMBER = re.compile(r"(?<![0-9.])[0-9]+(?:\.[0-9]+)?")


def _absent(value: Any) -> bool:
    return value is None or value == "" or value == []


def _normal(value: Any) -> Any:
    """Absent text (``None``, ``""`` or ``[]``) compares as ``None``."""
    return None if _absent(value) else value


def sentence_count(text: str) -> int:
    """Sentences in ``text``; abbreviations and decimals do not end one."""
    masked = text
    for abbreviation in ABBREVIATIONS:
        masked = masked.replace(abbreviation, abbreviation.replace(".", "\0"))
    # A single capital initial ("A. B. Smith") does not end a sentence.
    masked = re.sub(r"(?<![A-Za-z])([A-Z])\.(?=\s+[A-Z][a-z]*\.?)", "\\1\0", masked)
    ends = re.findall(r"[.!?]+[\"')\]]*(?=\s|$)", masked)
    tail = re.split(r"[.!?]+[\"')\]]*(?=\s|$)", masked)[-1].strip()
    return len(ends) + (1 if tail else 0)


def text_numbers(text: str) -> list[str]:
    """Number tokens the guard checks: decimals and integers of 2+ digits."""
    return [t for t in _NUMBER.findall(text) if "." in t or len(t) >= 2]


def _values(value: Any) -> Iterable[float]:
    if isinstance(value, bool):
        return
    if isinstance(value, (int, float)):
        yield abs(float(value))
    elif isinstance(value, str):
        for token in _NUMBER.findall(value):
            yield float(token)
    elif isinstance(value, Mapping):
        for key, item in value.items():
            yield from _values(key)
            yield from _values(item)
    elif isinstance(value, list):
        for item in value:
            yield from _values(item)


def structured_fields(record: Mapping[str, Any]) -> dict[str, Any]:
    """The record without its four free-text fields and its ``check`` block."""
    result = copy.deepcopy(dict(record))
    for key in ("label", "notes", "check", "superseded_record_text"):
        result.pop(key, None)
    for key in ("validity", "thermal"):
        if isinstance(result.get(key), dict):
            result[key].pop("notes", None)
    return result


def allowed_numbers(
    record: Mapping[str, Any], stored: Mapping[str, Any] | None = None
) -> set[float]:
    """Numbers in the record's current text, structured fields and stored originals."""
    texts = [record_text(dict(record), name) for name in SUPERSEDED_TEXT_FIELDS]
    originals = [v for k, v in (stored or {}).items() if k in SUPERSEDED_TEXT_FIELDS]
    return (
        set(_values(texts))
        | set(_values(originals))
        | set(_values(structured_fields(record)))
    )


@dataclass
class Refusal:
    material: str
    record: str
    reason: str
    detail: str
    guard: bool = False

    def line(self) -> str:
        kind = "guard " if self.guard else ""
        return (
            f"REFUSED {self.material}/{self.record}: {kind}{self.reason}: {self.detail}"
        )


@dataclass
class Plan:
    """What a run would change, refuse and override."""

    changes: dict[str, dict[str, dict[str, Any]]] = field(default_factory=dict)
    refusals: list[Refusal] = field(default_factory=list)
    allowed: list[str] = field(default_factory=list)
    notices: list[str] = field(default_factory=list)


def read_batches(paths: Iterable[Path], plan: Plan) -> dict[str, dict[str, dict]]:
    """Merge input files; a record named twice with different content is refused."""
    merged: dict[str, dict[str, dict]] = {}
    conflicts: set[tuple[str, str]] = set()
    for path in paths:
        data = json.loads(Path(path).read_text(encoding="utf-8"))
        if not isinstance(data, dict):
            raise ValueError(f"{path}: top level must be an object of materials")
        for material, records in data.items():
            if not isinstance(records, dict):
                raise ValueError(f"{path}: {material} must map record identifiers")
            for record, entry in records.items():
                if not isinstance(entry, dict) or not isinstance(
                    entry.get("before"), dict
                ):
                    raise ValueError(f"{path}: {material}/{record} needs before/after")
                if not isinstance(entry.get("after"), dict):
                    raise ValueError(f"{path}: {material}/{record} needs before/after")
                if "label" in entry["after"] or "label" in entry["before"]:
                    plan.notices.append(
                        f"{material}/{record}: label in {path.name} ignored; labels "
                        "come only from --labels"
                    )
                item = {
                    "before": {
                        k: entry["before"][k]
                        for k in JSON_FIELDS
                        if k in entry["before"]
                    },
                    "after": {
                        k: entry["after"][k] for k in JSON_FIELDS if k in entry["after"]
                    },
                }
                previous = merged.setdefault(material, {}).get(record)
                if previous is not None and previous != item:
                    conflicts.add((material, record))
                merged[material][record] = item
    for material, record in sorted(conflicts):
        plan.refusals.append(
            Refusal(material, record, "conflict", "named twice with different text")
        )
        del merged[material][record]
    return merged


def read_labels(
    path: Path | None,
) -> tuple[dict[str, dict[str, tuple[str, str]]], dict]:
    """Mechanical label rows by material and record, and a flag count."""
    labels: dict[str, dict[str, tuple[str, str]]] = {}
    counts: dict[str, int] = {}
    if path is None:
        return labels, counts
    with open(path, newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            flag = row["flag"].strip()
            counts[flag] = counts.get(flag, 0) + 1
            if flag != "mechanical":
                continue
            if row["proposed_label"] == row["current_label"]:
                continue
            labels.setdefault(row["material"], {})[row["record"]] = (
                row["current_label"],
                row["proposed_label"],
            )
    return labels, counts


def _receipt(root: Path) -> dict:
    path = root / RECEIPT.relative_to(ROOT)
    return json.loads(path.read_bytes()) if path.exists() else {"materials": {}}


def _guards(
    record: Mapping[str, Any],
    new: Mapping[str, Any],
    stored: Mapping[str, Any] | None = None,
) -> list[tuple[str, str]]:
    """Failed guards (name, detail) for the changed fields ``new``."""
    failed: list[tuple[str, str]] = []
    allowed = allowed_numbers(record, stored)
    texts: list[tuple[str, str]] = []
    for name, value in new.items():
        if _absent(value):
            continue
        for text in value if isinstance(value, list) else [value]:
            texts.append((name, text))
    unknown = sorted(
        {
            token
            for _, text in texts
            for token in text_numbers(text)
            if float(token) not in allowed
        }
    )
    if unknown:
        failed.append(("numbers", f"not in the record: {', '.join(unknown)}"))
    notes = new.get("notes")
    if not _absent(notes):
        if len(notes) > NOTES_MAX_CHARACTERS:
            failed.append(("notes_length", f"{len(notes)} characters"))
        elif sentence_count(notes) > NOTES_MAX_SENTENCES:
            failed.append(("notes_length", f"{sentence_count(notes)} sentences"))
    validity = new.get("validity_notes")
    if not _absent(validity):
        if len(validity) > VALIDITY_MAX_NOTES:
            failed.append(("validity_length", f"{len(validity)} notes"))
        longest = max(len(note) for note in validity)
        if longest > VALIDITY_NOTE_MAX_CHARACTERS:
            failed.append(("validity_length", f"a note of {longest} characters"))
    thermal = new.get("thermal_notes")
    if not _absent(thermal) and len(thermal) > THERMAL_NOTES_MAX_CHARACTERS:
        failed.append(("thermal_length", f"{len(thermal)} characters"))
    label = new.get("label")
    if not _absent(label) and len(label) > LABEL_MAX_CHARACTERS:
        failed.append(("label_length", f"{len(label)} characters"))
    found = sorted(
        {word for _, text in texts for word in FORBIDDEN_TEXT if word in text}
    )
    if found:
        failed.append(("forbidden_text", ", ".join(repr(word) for word in found)))
    return failed


def plan_changes(
    batches: Mapping[str, Mapping[str, dict]],
    labels: Mapping[str, Mapping[str, tuple[str, str]]],
    *,
    allow: Iterable[str] = (),
    root: Path = ROOT,
    plan: Plan | None = None,
) -> Plan:
    """Check every requested change against the files and the guards."""
    plan = Plan() if plan is None else plan
    overrides = set(allow)
    bindings = _receipt(root).get("materials", {})
    for material in sorted(set(batches) | set(labels)):
        requested = {
            **dict.fromkeys(labels.get(material, {})),
            **batches.get(material, {}),
        }
        record_path = root / f"peritheos/data/materials/{material}.eosmat"
        if not record_path.exists():
            for record in sorted(requested):
                plan.refusals.append(
                    Refusal(material, record, "unknown_material", "no file")
                )
            continue
        document = json.loads(record_path.read_bytes())
        sidecar = json.loads(
            (root / f"peritheos/data/audits/{material}.audit.json").read_bytes()
        )
        records = {r["identifier"]: r for r in document["eos_records"]}
        order = [r["identifier"] for r in document["eos_records"]]
        bound = set(bindings.get(material, {}).get("superseded_record_text", []))
        proposals: dict[str, dict[str, Any]] = {}
        for identifier in sorted(requested):
            record = records.get(identifier)
            if record is None:
                plan.refusals.append(
                    Refusal(material, identifier, "unknown_record", "not in the file")
                )
                continue
            reasons = []
            new: dict[str, Any] = {}
            entry = batches.get(material, {}).get(identifier)
            if entry is not None:
                for name in JSON_FIELDS:
                    if name not in entry["after"]:
                        continue
                    if name not in entry["before"]:
                        reasons.append(
                            ("before_missing", f"no before value for {name}")
                        )
                        continue
                    current = record_text(record, name)
                    if _normal(entry["before"][name]) != _normal(current):
                        reasons.append(
                            ("before_mismatch", f"{name} differs from the file")
                        )
                        continue
                    if _normal(entry["after"][name]) != _normal(current):
                        new[name] = _normal(entry["after"][name])
            label = labels.get(material, {}).get(identifier)
            if label is not None:
                if label[0] != record.get("label"):
                    reasons.append(
                        ("label_mismatch", "CSV current_label differs from the file")
                    )
                elif label[1] != record.get("label"):
                    new["label"] = label[1]
            stored = (
                sidecar["records"].get(identifier, {}).get("superseded_record_text")
            )
            if (
                new
                and stored is not None
                and material in bindings
                and order.index(identifier) not in bound
            ):
                reasons.append(
                    (
                        "predates_split",
                        "the stored originals predate the audit split; re-split first",
                    )
                )
            if (
                new
                and "thermal_notes" in new
                and not isinstance(record.get("thermal"), dict)
            ):
                reasons.append(("no_thermal", "the record has no thermal object"))
            for reason, detail in reasons:
                plan.refusals.append(Refusal(material, identifier, reason, detail))
            if reasons or not new:
                continue
            refused = False
            for guard, detail in _guards(record, new, stored):
                if f"{identifier}:{guard}" in overrides:
                    plan.allowed.append(f"{material}/{identifier}:{guard} ({detail})")
                    continue
                plan.refusals.append(Refusal(material, identifier, guard, detail, True))
                refused = True
            if not refused:
                proposals[identifier] = new
        # Labels must stay unique within the material once the changes apply.
        while True:
            final = {
                identifier: proposals.get(identifier, {}).get(
                    "label", record.get("label")
                )
                for identifier, record in records.items()
            }
            clashes = [
                identifier
                for identifier, change in proposals.items()
                if "label" in change
                and list(final.values()).count(change["label"]) > 1
                and f"{identifier}:label_unique" not in overrides
            ]
            if not clashes:
                break
            for identifier in clashes:
                plan.refusals.append(
                    Refusal(
                        material,
                        identifier,
                        "label_unique",
                        f"label {proposals[identifier]['label']!r} is not unique",
                        True,
                    )
                )
                del proposals[identifier]
        for identifier, change in proposals.items():
            if "label" in change and list(final.values()).count(change["label"]) > 1:
                plan.allowed.append(f"{material}/{identifier}:label_unique")
        if proposals:
            plan.changes[material] = {i: proposals[i] for i in order if i in proposals}
    return plan


def revise_record(record: dict, entry: dict, new: Mapping[str, Any], date: str) -> None:
    """Write ``new`` text into ``record``; keep each replaced original in ``entry``.

    An original already stored is never overwritten. A field rewritten back
    to its stored original drops that original; an emptied block is removed.
    ``superseded_on`` is the date of the latest rewrite that stored an
    original.
    """
    block = dict(entry.get("superseded_record_text", {}))
    added = False
    for name, value in new.items():
        if name not in block:
            block[name] = record_text(record, name)
            added = True
        elif _normal(block[name]) == _normal(value):
            del block[name]
        set_record_text(record, name, value)
    if added or "superseded_on" not in block:
        block["superseded_on"] = date
    ordered = {
        k: block[k] for k in (*SUPERSEDED_TEXT_FIELDS, "superseded_on") if k in block
    }
    entry.pop("superseded_record_text", None)
    if any(k in ordered for k in SUPERSEDED_TEXT_FIELDS):
        validate_superseded_record_text(ordered, str(record.get("identifier")))
        entry["superseded_record_text"] = ordered


def apply_plan(plan: Plan, date: str, root: Path = ROOT) -> None:
    """Write the planned changes, bind them in the receipt and verify."""
    revised: dict[str, list[int]] = {}
    copies: dict[str, dict[str, Any]] = {}
    bindings = _receipt(root).get("materials", {})
    for material, changes in plan.changes.items():
        record_path = root / f"peritheos/data/materials/{material}.eosmat"
        audit_path = root / f"peritheos/data/audits/{material}.audit.json"
        record_text_before = record_path.read_text(encoding="utf-8")
        audit_text_before = audit_path.read_text(encoding="utf-8")
        document = json.loads(record_text_before)
        sidecar = json.loads(audit_text_before)
        if (
            format_eosmat(document) != record_text_before
            or format_audit(sidecar) != audit_text_before
        ):
            raise ValueError(f"{material}: files are not in canonical layout")
        order = [r["identifier"] for r in document["eos_records"]]
        entries = sidecar["records"]
        for record in document["eos_records"]:
            identifier = record["identifier"]
            if identifier not in changes:
                continue
            entry = entries.setdefault(identifier, {})
            revise_record(record, entry, changes[identifier], date)
            if not entry:
                del entries[identifier]
        sidecar["records"] = {i: entries[i] for i in order if i in entries}
        record_bytes = format_eosmat(document)
        audit_bytes = format_audit(sidecar)
        merged = merge_audit(json.loads(record_bytes), json.loads(audit_bytes))
        validate_eosmat_document(json.loads(record_bytes))
        validate_eosmat_document(merged)
        record_path.write_text(record_bytes, encoding="utf-8")
        audit_path.write_text(audit_bytes, encoding="utf-8")
        bound = set(bindings.get(material, {}).get("superseded_record_text", []))
        for identifier in changes:
            index = order.index(identifier)
            if "superseded_record_text" in sidecar["records"].get(identifier, {}):
                bound.add(index)
            else:
                bound.discard(index)
        revised[material] = sorted(bound)
        for record in document["eos_records"]:
            if record["identifier"] in changes:
                copies[record["identifier"]] = {
                    "label": record.get("label"),
                    "source_notes": record.get("notes"),
                }
    bind_record_text_revision(revised, root)
    sync_ledger(copies, root)


def sync_ledger(copies: Mapping[str, Mapping[str, Any]], root: Path = ROOT) -> None:
    """Copy revised labels and notes into the registered refit ledger.

    ``docs/data/primary-eos-refits.json`` keeps each record's ``label`` and
    ``notes`` (as ``source_notes``), and ``validate_primary_eos_refits.py``
    requires them to equal the catalog. Only these two fields of the revised
    records change; ``docs/primary-eos-refits.md`` is re-rendered from the
    saved ledger by that script's own renderer. Nothing is recomputed.
    """
    path = root / LEDGER
    if not copies or not path.exists():
        return
    ledger = json.loads(path.read_text(encoding="utf-8"))
    for entry in ledger["records"]:
        if entry["record_identifier"] in copies:
            entry.update(copies[entry["record_identifier"]])
    path.write_text(
        json.dumps(ledger, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    from scripts.validate_primary_eos_refits import render_markdown

    (root / LEDGER_MARKDOWN).write_text(render_markdown(ledger), encoding="utf-8")


def _git_file(revision: str, relative: str, root: Path) -> bytes | None:
    result = subprocess.run(
        ["git", "show", f"{revision}:{relative}"],
        cwd=root,
        capture_output=True,
        check=False,
    )
    return result.stdout if result.returncode == 0 else None


def verify(
    revision: str = DEFAULT_REVISION,
    root: Path = ROOT,
    read: Callable[[str], bytes | None] | None = None,
) -> tuple[int, list[str]]:
    """Restore stored originals and compare the text with ``revision``.

    Returns the number of records with a ``superseded_record_text`` block and
    the failures. A record present at the revision without a block must have
    unchanged text. ``read`` returns a repository file at the revision
    (default: ``git show``).
    """
    if read is None:

        def read(relative: str) -> bytes | None:
            return _git_file(revision, relative, root)

    failures: list[str] = []
    checked = 0
    ledger_path = root / LEDGER
    ledger = (
        {
            entry["record_identifier"]: entry
            for entry in json.loads(ledger_path.read_bytes())["records"]
        }
        if ledger_path.exists()
        else {}
    )
    for path in sorted((root / "peritheos/data/materials").glob("*.eosmat")):
        material = path.stem
        current = json.loads(path.read_bytes())
        sidecar = json.loads(
            (root / f"peritheos/data/audits/{material}.audit.json").read_bytes()
        )
        old_bytes = read(f"peritheos/data/materials/{material}.eosmat")
        old = (
            {}
            if old_bytes is None
            else {r["identifier"]: r for r in json.loads(old_bytes)["eos_records"]}
        )
        for record in current["eos_records"]:
            identifier = record["identifier"]
            block = sidecar["records"].get(identifier, {}).get("superseded_record_text")
            restored = copy.deepcopy(record)
            if block is not None:
                checked += 1
                validate_superseded_record_text(block, identifier)
                copy_in_ledger = ledger.get(identifier)
                if copy_in_ledger is not None and (
                    copy_in_ledger.get("label") != record.get("label")
                    or copy_in_ledger.get("source_notes") != record.get("notes")
                ):
                    failures.append(
                        f"{material}/{identifier}: label or notes differ in {LEDGER}"
                    )
                restored["superseded_record_text"] = block
                restore_record_text(restored)
            if identifier not in old:
                if block is not None:
                    failures.append(f"{material}/{identifier}: not in {revision}")
                continue
            for name in SUPERSEDED_TEXT_FIELDS:
                if record_text(restored, name) != record_text(old[identifier], name):
                    what = "restored" if block is not None else "unrecorded change of"
                    failures.append(
                        f"{material}/{identifier}: {what} {name} differs from {revision}"
                    )
    return checked, failures


def _size(value: Any) -> int:
    if _absent(value):
        return 0
    return len(value) if isinstance(value, str) else sum(map(len, value))


def _describe(
    material: str, identifier: str, new: Mapping[str, Any], record: dict
) -> str:
    parts = [
        f"{name} {_size(record_text(record, name))}->{_size(value)} characters"
        for name, value in new.items()
    ]
    return f"{material}/{identifier}: " + ", ".join(parts)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batches", nargs="*", type=Path, help="input JSON files")
    parser.add_argument("--labels", type=Path, help="label CSV (mechanical rows apply)")
    parser.add_argument("--check", action="store_true", help="dry run")
    parser.add_argument(
        "--allow",
        action="append",
        default=[],
        metavar="RECORD:GUARD",
        help=f"override one guard for one record ({', '.join(GUARDS)})",
    )
    parser.add_argument("--date", default=datetime.date.today().isoformat())
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--revision", default=DEFAULT_REVISION)
    parser.add_argument("--report", type=Path, help="write a JSON report")
    arguments = parser.parse_args(argv)
    if arguments.verify:
        checked, failures = verify(arguments.revision)
        for failure in failures:
            print(f"FAILED {failure}")
        print(
            f"Verified {checked} records with superseded_record_text against "
            f"{arguments.revision}: {'ok' if not failures else f'{len(failures)} failures'}."
        )
        return 1 if failures else 0
    for item in arguments.allow:
        record, _, guard = item.partition(":")
        if not record or guard not in GUARDS:
            parser.error(f"--allow {item!r}: expected RECORD:GUARD, GUARD in {GUARDS}")
    plan = Plan()
    batches = read_batches(arguments.batches, plan)
    labels, counts = read_labels(arguments.labels)
    plan_changes(batches, labels, allow=arguments.allow, plan=plan)
    for notice in plan.notices:
        print(f"NOTE {notice}")
    verb, tag = (
        ("Would change", "WOULD-CHANGE") if arguments.check else ("Changed", "CHANGED")
    )
    lines = []
    for material, changes in plan.changes.items():
        document = json.loads(
            (ROOT / f"peritheos/data/materials/{material}.eosmat").read_bytes()
        )
        records = {r["identifier"]: r for r in document["eos_records"]}
        for identifier, new in changes.items():
            lines.append(_describe(material, identifier, new, records[identifier]))
    if not arguments.check:
        apply_plan(plan, arguments.date)
    for line in lines:
        print(f"{tag} {line}")
    for allowed in plan.allowed:
        print(f"ALLOWED {allowed}")
    for refusal in plan.refusals:
        print(refusal.line())
    fields = sum(len(c) for m in plan.changes.values() for c in m.values())
    print(
        f"{verb} {sum(len(m) for m in plan.changes.values())} records ({fields} fields) "
        f"in {len(plan.changes)} materials; {len(plan.refusals)} refusals; "
        f"{len(plan.allowed)} overrides."
        + (f" Label CSV flags: {counts}." if counts else "")
    )
    if arguments.report:
        arguments.report.write_text(
            json.dumps(
                {
                    "changes": plan.changes,
                    "refusals": [r.__dict__ for r in plan.refusals],
                    "allowed": plan.allowed,
                    "notices": plan.notices,
                    "label_flags": counts,
                },
                indent=2,
                ensure_ascii=False,
            )
            + "\n",
            encoding="utf-8",
        )
    return 1 if plan.refusals else 0


if __name__ == "__main__":
    sys.exit(main())
