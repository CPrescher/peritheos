"""Repository maintenance: write the curated check texts of bundled records.

Two short texts in a record's audit sidecar entry
(``peritheos/data/audits/<id>.audit.json``) are copied into the record's
``check`` block in ``peritheos/data/materials/<id>.eosmat``
(``docs/record-text-standard.md``):

``source_note``
    written to ``evidence.scientific_validation.primary_source_check.record_note``
    and shown as ``check.source_note`` in place of the verbatim ``finding``.
    ``finding`` itself is never changed; ``null`` removes the note, so
    ``check.source_note`` is the verbatim ``finding`` again.
``warning_reason``
    written to ``evidence.scientific_validation.warning_reason`` and copied to
    ``check.warning_reason``; ``null`` removes it. ``check.warning`` wording
    is unchanged.

Input files (JSON), any number::

    {"<material>": {"<record>": {"source_note": str|null,
                                 "warning_reason": str|null,
                                 ...other keys are ignored}}}

Either key may be absent. The record's ``check`` block is regenerated from the
merged evidence (only ``source_note`` and ``warning_reason`` may change), both
files are written through the repository formatter and the audit-split
receipt binds the revision: the value each field had before its first
post-split change is kept as a ``check_text`` entry, so
``scripts/split_eosmat_audits.py --verify`` still recovers every pre-split
file byte for byte.

A record is refused, and reported, when its material or record is unknown,
its material's files are not in canonical layout (``not_canonical``; the
dry run reports this too), a value is neither text nor ``null``, the record
has no ``primary_source_check`` (for ``source_note``), or a
``warning_reason`` is given for a record whose ``check.warning`` is ``null``
(the library rejects such a record, so this cannot be overridden). These
guards refuse it too:

``source_note_length``
    a source note over 300 characters;
``warning_reason_length``
    a warning reason over 140 characters;
``numbers``
    a number (decimals, numbers with an exponent, or integers of two or more
    digits; sign, mantissa and exponent read as one number, ``−`` as minus;
    the value of a printed form such as ``58.868(2)``; also directly after
    ``p.``, ``pp.``, ``Eq.``, ``Eqs.``, ``Fig.``, ``Table.`` or ``ca.``) that
    neither equals, sign included, a number in the record's ``finding`` and
    other sidecar evidence, its record fields or its stored
    ``superseded_record_text``, nor is a correct rounding of one at the
    precision written with at least two significant digits (``6.9`` for
    ``6.880``, ``23.89`` for ``23.89043``);
``forbidden_text``
    a repository path (``docs/``, ``.json``, ``scripts/``) or the words
    ``LitCurate`` or ``Audit 20``;
``shared_note``
    a source note identical to the note another record of the same material
    shows after the run (its new or existing ``record_note``, or its
    verbatim ``finding``), so a shared note is a conscious exception.

``--allow RECORD:GUARD`` overrides one guard for one named record; every
override is listed in the output. A refused record is left unchanged; the
others are applied and the exit status is 1. Every file and the receipt
binding are computed and verified before anything is written, so an error
leaves the tree untouched.

Usage::

    python scripts/apply_check_text.py BATCH.json ... --check   # dry run
    python scripts/apply_check_text.py BATCH.json ...           # apply
    python scripts/apply_check_text.py --verify

``--verify`` confirms that every bundled ``check`` block equals its
derivation from the merged sidecar (so ``source_note`` and ``warning_reason``
are current), that every curated text is within its limit and that every
material holding one still recovers its pre-split file exactly. This is not
an installed-package tool.
"""

from __future__ import annotations

import argparse
import copy
import json
import sys
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from peritheos._eosmat_audit import (  # noqa: E402
    RECORD_NOTE_MAX_LENGTH,
    WARNING_REASON_MAX_LENGTH,
    check_block,
    merge_audit,
    source_note,
    usage_warning,
)
from peritheos._eosmat_format import format_audit, format_eosmat  # noqa: E402
from peritheos.eosmat import validate_eosmat_document  # noqa: E402
from scripts.apply_record_text import (  # noqa: E402
    FORBIDDEN_TEXT,
    Plan,
    Refusal,
    _receipt,
    _values,
    known_number,
    text_numbers,
)
from scripts.split_eosmat_audits import (  # noqa: E402
    CHECK_TEXT_FIELDS,
    RECEIPT,
    check_text,
    check_text_receipt,
    pre_split_bytes,
    receipt_text,
    set_check_text,
)

#: Input keys and the sidecar field each one writes.
INPUT_FIELDS = {"source_note": "record_note", "warning_reason": "warning_reason"}
LIMITS = {
    "source_note": RECORD_NOTE_MAX_LENGTH,
    "warning_reason": WARNING_REASON_MAX_LENGTH,
}
GUARDS = (
    "source_note_length",
    "warning_reason_length",
    "numbers",
    "forbidden_text",
    "shared_note",
)
_ABSENT = object()


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
                if not isinstance(entry, dict):
                    raise ValueError(f"{path}: {material}/{record} must be an object")
                item = {k: entry[k] for k in INPUT_FIELDS if k in entry}
                if not item:
                    plan.notices.append(
                        f"{material}/{record}: no source_note or warning_reason in "
                        f"{Path(path).name}"
                    )
                    continue
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


def _paths(root: Path, material: str) -> tuple[Path, Path]:
    return (
        root / f"peritheos/data/materials/{material}.eosmat",
        root / f"peritheos/data/audits/{material}.audit.json",
    )


def _read(root: Path, material: str) -> tuple[dict, dict]:
    """A material's record file and sidecar; both must be in canonical layout."""
    record_path, audit_path = _paths(root, material)
    record_text = record_path.read_text(encoding="utf-8")
    audit_text = audit_path.read_text(encoding="utf-8")
    document = json.loads(record_text)
    sidecar = json.loads(audit_text)
    if format_eosmat(document) != record_text or format_audit(sidecar) != audit_text:
        raise ValueError(f"{material}: files are not in canonical layout")
    return document, sidecar


def allowed_numbers(record: Mapping[str, Any]) -> set[float]:
    """Numbers in a merged record: evidence, record fields and stored originals.

    The ``check`` block and the curated texts themselves do not count.
    """
    source = copy.deepcopy(dict(record))
    source.pop("check", None)
    for field in CHECK_TEXT_FIELDS:
        set_check_text(source, field, None)
    return set(_values(source))


def _guards(
    record: Mapping[str, Any], new: Mapping[str, str | None]
) -> list[tuple[str, str]]:
    """Failed guards (name, detail) for the new texts, without ``shared_note``."""
    failed: list[tuple[str, str]] = []
    texts = [(key, value) for key, value in new.items() if value is not None]
    for key, value in texts:
        if len(value) > LIMITS[key]:
            failed.append(
                (f"{key}_length", f"{len(value)} characters (limit {LIMITS[key]})")
            )
    allowed = allowed_numbers(record)
    unknown = sorted(
        {
            token
            for _, value in texts
            for token in text_numbers(value)
            if not known_number(token, allowed, rounding=True)
        }
    )
    if unknown:
        failed.append(("numbers", f"not in the record: {', '.join(unknown)}"))
    found = sorted(
        {word for _, value in texts for word in FORBIDDEN_TEXT if word in value}
    )
    if found:
        failed.append(("forbidden_text", ", ".join(repr(word) for word in found)))
    return failed


def plan_changes(
    batches: Mapping[str, Mapping[str, dict]],
    *,
    allow: Iterable[str] = (),
    root: Path = ROOT,
    plan: Plan | None = None,
) -> Plan:
    """Check every requested change against the files and the guards."""
    plan = Plan() if plan is None else plan
    overrides = set(allow)

    def guard(material: str, identifier: str, name: str, detail: str) -> bool:
        """Record a failed guard; ``True`` when it refuses the record."""
        if f"{identifier}:{name}" in overrides:
            plan.allowed.append(f"{material}/{identifier}:{name} ({detail})")
            return False
        plan.refusals.append(Refusal(material, identifier, name, detail, True))
        return True

    for material in sorted(batches):
        requested = batches[material]
        record_path, audit_path = _paths(root, material)
        if not record_path.exists() or not audit_path.exists():
            for identifier in sorted(requested):
                plan.refusals.append(
                    Refusal(material, identifier, "unknown_material", "no file")
                )
            continue
        try:
            document, sidecar = _read(root, material)
        except ValueError as error:
            # The real run would rewrite the whole layout; refuse instead.
            for identifier in sorted(requested):
                plan.refusals.append(
                    Refusal(material, identifier, "not_canonical", str(error))
                )
            continue
        merged = merge_audit(document, sidecar)
        records = {r["identifier"]: r for r in merged["eos_records"]}
        order = [r["identifier"] for r in merged["eos_records"]]
        proposals: dict[str, dict[str, str | None]] = {}
        for identifier in sorted(requested):
            record = records.get(identifier)
            if record is None:
                plan.refusals.append(
                    Refusal(material, identifier, "unknown_record", "not in the file")
                )
                continue
            reasons: list[tuple[str, str]] = []
            new: dict[str, str | None] = {}
            for key, field in INPUT_FIELDS.items():
                value = requested[identifier].get(key, _ABSENT)
                if value is _ABSENT:
                    continue
                if value is not None and (not isinstance(value, str) or not value):
                    reasons.append(("not_text", f"{key} must be text or null"))
                    continue
                if value is not None and value != value.strip():
                    reasons.append(
                        ("not_text", f"{key} has leading or trailing white space")
                    )
                    continue
                if value != check_text(record, field):
                    new[key] = value
            if new.get("source_note") is not None and not isinstance(
                record.get("evidence", {})
                .get("scientific_validation", {})
                .get("primary_source_check"),
                dict,
            ):
                reasons.append(("no_source_check", "the record has no source check"))
            if (
                new.get("warning_reason") is not None
                and usage_warning(record, merged) is None
            ):
                reasons.append(
                    (
                        "no_warning",
                        "the record has no warning (check.warning is null); a "
                        "warning_reason is only given with a warning",
                    )
                )
            for reason, detail in reasons:
                plan.refusals.append(Refusal(material, identifier, reason, detail))
            if reasons or not new:
                continue
            refused = False
            for name, detail in _guards(record, new):
                refused |= guard(material, identifier, name, detail)
            if not refused:
                proposals[identifier] = new
        # A source note must differ from every other record's shown note.
        while True:
            shown: dict[str, str | None] = {}
            for identifier, record in records.items():
                trial = copy.deepcopy(record)
                if "source_note" in proposals.get(identifier, {}):
                    set_check_text(
                        trial, "record_note", proposals[identifier]["source_note"]
                    )
                shown[identifier] = source_note(trial)
            clashes = [
                identifier
                for identifier, change in proposals.items()
                if change.get("source_note") is not None
                and f"{identifier}:shared_note" not in overrides
                and any(
                    other != identifier and note == change["source_note"]
                    for other, note in shown.items()
                )
            ]
            if not clashes:
                break
            for identifier in clashes:
                others = sorted(
                    other
                    for other, note in shown.items()
                    if other != identifier
                    and note == proposals[identifier]["source_note"]
                )
                guard(material, identifier, "shared_note", f"same as {others}")
                del proposals[identifier]
        for identifier, change in proposals.items():
            if change.get("source_note") is not None and any(
                other != identifier and note == change["source_note"]
                for other, note in shown.items()
            ):
                plan.allowed.append(f"{material}/{identifier}:shared_note")
        if proposals:
            plan.changes[material] = {i: proposals[i] for i in order if i in proposals}
    return plan


def apply_plan(plan: Plan, root: Path = ROOT) -> None:
    """Write the planned changes, regenerate ``check``, bind and verify.

    Every material's new files and the receipt binding are computed and
    verified in memory first; files are written only when the whole plan is
    valid, so a failure leaves the tree untouched.
    """
    receipt = _receipt(root)
    bindings = receipt.get("materials", {})
    revised: dict[str, list] = {}
    files: dict[str, tuple[bytes, bytes]] = {}
    for material, changes in plan.changes.items():
        document, sidecar = _read(root, material)
        order = [r["identifier"] for r in document["eos_records"]]
        originals = {
            index: dict(values)
            for index, values in bindings.get(material, {}).get("check_text", [])
        }
        for identifier, new in changes.items():
            entry = sidecar["records"].get(identifier)
            if entry is None:
                raise ValueError(f"{material}/{identifier}: no sidecar entry")
            index = order.index(identifier)
            stored = originals.setdefault(index, {})
            for key, value in new.items():
                field = INPUT_FIELDS[key]
                stored.setdefault(field, check_text(entry, field))
                set_check_text(entry, field, value)
                if stored[field] == value:
                    del stored[field]
            if not stored:
                del originals[index]
        merged = merge_audit(document, sidecar)
        complete = {r["identifier"]: r for r in merged["eos_records"]}
        for record in document["eos_records"]:
            identifier = record["identifier"]
            if identifier not in changes:
                continue
            old = record["check"]
            derived = check_block(
                {k: v for k, v in complete[identifier].items() if k != "check"},
                merged,
                audit=old.get("audit"),
            )
            kept = ("source_note", "warning_reason")
            if {k: v for k, v in old.items() if k not in kept} != {
                k: v for k, v in derived.items() if k not in kept
            }:
                raise ValueError(
                    f"{material}/{identifier}: the check block differs from its "
                    "evidence beyond source_note and warning_reason"
                )
            record["check"] = derived
        record_text = format_eosmat(document)
        audit_text = format_audit(sidecar)
        validate_eosmat_document(json.loads(record_text))
        validate_eosmat_document(
            merge_audit(json.loads(record_text), json.loads(audit_text))
        )
        files[material] = (record_text.encode("utf-8"), audit_text.encode("utf-8"))
        revised[material] = [[index, originals[index]] for index in sorted(originals)]
    # Raises before anything is written if a pre-split file would be lost.
    new_receipt = check_text_receipt(receipt, revised, files) if bindings else None
    for material, (record_bytes, audit_bytes) in files.items():
        record_path, audit_path = _paths(root, material)
        record_path.write_bytes(record_bytes)
        audit_path.write_bytes(audit_bytes)
    if new_receipt is not None:
        (root / RECEIPT.relative_to(ROOT)).write_text(
            receipt_text(new_receipt), encoding="utf-8"
        )


def verify(root: Path = ROOT) -> tuple[int, list[str]]:
    """Every ``check`` equals its derivation; curated texts recover losslessly.

    Returns the number of records holding a curated text and the failures.
    """
    failures: list[str] = []
    curated = 0
    bound = _receipt(root).get("materials", {})
    for path in sorted((root / "peritheos/data/materials").glob("*.eosmat")):
        material = path.stem
        document = json.loads(path.read_bytes())
        audit_path = root / f"peritheos/data/audits/{material}.audit.json"
        merged = merge_audit(document, json.loads(audit_path.read_bytes()))
        try:
            validate_eosmat_document(merged)
        except Exception as error:  # noqa: BLE001 - reported, not raised
            failures.append(f"{material}: {error}")
            continue
        holds = False
        for record, complete in zip(document["eos_records"], merged["eos_records"]):
            identifier = record["identifier"]
            if any(check_text(complete, f) is not None for f in CHECK_TEXT_FIELDS):
                curated += 1
                holds = True
            check = record.get("check")
            if not isinstance(check, dict):
                continue
            derived = check_block(
                {k: v for k, v in complete.items() if k != "check"},
                merged,
                audit=check.get("audit"),
            )
            if check != derived:
                keys = sorted(
                    k
                    for k in set(check) | set(derived)
                    if check.get(k) != derived.get(k)
                )
                failures.append(
                    f"{material}/{identifier}: check differs from its evidence "
                    f"in {keys}"
                )
        if (holds or "check_text" in bound.get(material, {})) and material in bound:
            try:
                pre_split_bytes(material, root=root)
            except ValueError as error:
                failures.append(f"{material}: {error}")
    return curated, failures


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("batches", nargs="*", type=Path, help="input JSON files")
    parser.add_argument("--check", action="store_true", help="dry run")
    parser.add_argument(
        "--allow",
        action="append",
        default=[],
        metavar="RECORD:GUARD",
        help=f"override one guard for one record ({', '.join(GUARDS)})",
    )
    parser.add_argument("--verify", action="store_true")
    parser.add_argument("--report", type=Path, help="write a JSON report")
    arguments = parser.parse_args(argv)
    if arguments.verify:
        curated, failures = verify()
        for failure in failures:
            print(f"FAILED {failure}")
        print(
            f"Verified every check block; {curated} records hold a curated check "
            f"text: {'ok' if not failures else f'{len(failures)} failures'}."
        )
        return 1 if failures else 0
    for item in arguments.allow:
        record, _, name = item.partition(":")
        if not record or name not in GUARDS:
            parser.error(f"--allow {item!r}: expected RECORD:GUARD, GUARD in {GUARDS}")
    plan = Plan()
    batches = read_batches(arguments.batches, plan)
    plan_changes(batches, allow=arguments.allow, plan=plan)
    for notice in plan.notices:
        print(f"NOTE {notice}")
    verb, tag = (
        ("Would change", "WOULD-CHANGE") if arguments.check else ("Changed", "CHANGED")
    )
    if not arguments.check:
        apply_plan(plan)
    for material, changes in plan.changes.items():
        for identifier, new in changes.items():
            parts = [
                f"{key} {'removed' if value is None else f'{len(value)} characters'}"
                for key, value in new.items()
            ]
            print(f"{tag} {material}/{identifier}: {', '.join(parts)}")
    for allowed in plan.allowed:
        print(f"ALLOWED {allowed}")
    for refusal in plan.refusals:
        print(refusal.line())
    fields = sum(len(c) for m in plan.changes.values() for c in m.values())
    print(
        f"{verb} {sum(len(m) for m in plan.changes.values())} records ({fields} fields) "
        f"in {len(plan.changes)} materials; {len(plan.refusals)} refusals; "
        f"{len(plan.allowed)} overrides."
    )
    if arguments.report:
        arguments.report.write_text(
            json.dumps(
                {
                    "changes": plan.changes,
                    "refusals": [r.__dict__ for r in plan.refusals],
                    "allowed": plan.allowed,
                    "notices": plan.notices,
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
