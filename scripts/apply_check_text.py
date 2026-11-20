"""Repository maintenance: write the curated check texts of bundled records.

Two short texts and a coded usage recommendation in a record's audit sidecar
entry (``peritheos/data/audits/<id>.audit.json``) are copied into, or
generate, the record's ``check`` block in
``peritheos/data/materials/<id>.eosmat`` (``docs/record-text-standard.md``):

``source_note``
    written to ``evidence.scientific_validation.primary_source_check.record_note``
    and shown as ``check.source_note`` in place of the verbatim ``finding``.
    ``finding`` itself is never changed; ``null`` removes the note, so
    ``check.source_note`` is the verbatim ``finding`` again.
``warning_reason``
    written to ``evidence.scientific_validation.warning_reason`` and copied to
    ``check.warning_reason``; ``null`` removes it. ``check.warning`` wording
    is unchanged.
``usage_recommendation``
    one of the coded values of
    :data:`peritheos._eosmat_audit.USAGE_RECOMMENDATION_WARNINGS`, written to
    ``evidence.scientific_validation.usage_recommendation``; ``check.warning``
    is regenerated from it (fixed wording). ``null`` returns the field to its
    stored original. Free text replaced by a code is kept verbatim as
    ``scientific_validation.recommendation_reason`` (see below).

Input files (JSON), any number::

    {"<material>": {"<record>": {"source_note": str|null,
                                 "warning_reason": str|null,
                                 "usage_recommendation": str|null,
                                 ...other keys are ignored}}}

Any key may be absent. The record's ``check`` block is regenerated from the
merged evidence (only ``source_note``, ``warning_reason`` and, with
``usage_recommendation``, ``warning`` may change), both files are written
through the repository formatter and the audit-split receipt binds the
revision: the value each field had before its first post-split change is kept
as a ``check_text`` entry, so ``scripts/split_eosmat_audits.py --verify``
still recovers every pre-split file byte for byte.

What a usage recommendation carries with it, in both directions:

- Replacing free text by a code moves the text, verbatim, to
  ``recommendation_reason`` when that field is absent (an equal
  ``recommendation_reason`` is left as is). A different, existing
  ``recommendation_reason`` is refused (``recommendation_reason_present``): the
  curator merges the two texts by hand first. Nothing is ever shortened.
- ``evidence.reproduction.assessment.qualifications`` is rebuilt with
  ``scripts/standardize_reconstruction_metadata.py``'s
  ``assessment_for_record`` from the record and its row of the registered
  refit ledger (``docs/data/primary-eos-refits.json``), so the
  ``Usage recommendation: <text>`` entry becomes the moved text plus
  ``Usage recommendation: <code>``, exactly as in records that already carry
  a code. A change that would alter any other part of the assessment is
  refused (``beyond_fields``).
- The ledger row's copies (``usage_recommendation``, ``recommendation_reason``
  and the ``recommendation_reason`` appended to ``qualification``, as
  ``validate_primary_eos_refits.py`` writes them) are updated and
  ``docs/primary-eos-refits.md`` is re-rendered from the saved ledger;
  nothing is recomputed.
- The receipt keeps the originals of ``usage_recommendation``,
  ``recommendation_reason`` and ``assessment_qualifications``; ``null``, or
  setting the original code again, puts all three back.

A record is refused, and reported, when its material or record is unknown,
its material's files are not in canonical layout (``not_canonical``; the
dry run reports this too), a value is neither text nor ``null``, the record
has no ``primary_source_check`` (for ``source_note``), a
``usage_recommendation`` is not a coded value (``not_coded``) or the record
has no ledger row (``no_ledger_row``), the resulting record would keep a
``warning_reason`` without a warning (``no_warning``), or a change gives the
record its first warning without a ``warning_reason`` (``reason_required``).
These cannot be overridden (the library rejects a reason without a warning).
These guards refuse it too:

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

Every changed record's resulting ``check.warning`` is printed, in the dry run
too. ``--verify`` confirms that every bundled ``check`` block equals its
derivation from the merged sidecar (so ``source_note``, ``warning`` and
``warning_reason`` are current), that every material holding a bound field
still recovers its pre-split file exactly and that the ledger's
``usage_recommendation`` and ``recommendation_reason`` copies equal the
sidecars. This is not an installed-package tool.
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
    USAGE_RECOMMENDATION_WARNINGS,
    WARNING_REASON_MAX_LENGTH,
    check_block,
    merge_audit,
    source_note,
    usage_warning,
    warning_reason,
)
from peritheos._eosmat_format import format_audit, format_eosmat  # noqa: E402
from peritheos.eosmat import validate_eosmat_document  # noqa: E402
from scripts.apply_record_text import (  # noqa: E402
    FORBIDDEN_TEXT,
    LEDGER,
    LEDGER_MARKDOWN,
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
    RECOMMENDATION_FIELDS,
    check_text,
    check_text_receipt,
    pre_split_bytes,
    receipt_text,
    set_check_text,
)
from scripts.standardize_reconstruction_metadata import (  # noqa: E402
    assessment_for_record,
)

#: Text input keys and the sidecar field each one writes.
INPUT_FIELDS = {"source_note": "record_note", "warning_reason": "warning_reason"}
#: The input key (and sidecar field) of the coded usage recommendation.
USAGE_KEY = "usage_recommendation"
INPUT_KEYS = (*INPUT_FIELDS, USAGE_KEY)
#: Accepted usage recommendations: the codes with a fixed ``check.warning``.
CODED_RECOMMENDATIONS = tuple(USAGE_RECOMMENDATION_WARNINGS)
#: Ledger row fields copied from ``scientific_validation``.
LEDGER_FIELDS = ("usage_recommendation", "recommendation_reason")
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


class _Refused(Exception):
    """A change the files cannot take; reported as a refusal."""

    def __init__(self, reason: str, detail: str) -> None:
        super().__init__(f"{reason}: {detail}")
        self.reason = reason
        self.detail = detail


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
                item = {k: entry[k] for k in INPUT_KEYS if k in entry}
                if not item:
                    plan.notices.append(
                        f"{material}/{record}: no source_note, warning_reason or "
                        f"usage_recommendation in {Path(path).name}"
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


def _ledger(root: Path) -> dict | None:
    """The registered refit ledger, or ``None`` when the root has none."""
    path = root / LEDGER
    return json.loads(path.read_text(encoding="utf-8")) if path.exists() else None


def _ledger_rows(ledger: Mapping[str, Any] | None) -> dict[str, dict]:
    if ledger is None:
        return {}
    return {row["record_identifier"]: row for row in ledger["records"]}


def _originals(binding: Mapping[str, Any]) -> dict[int, dict]:
    """A material's ``check_text`` originals by record index."""
    return {index: dict(values) for index, values in binding.get("check_text", [])}


def _is_code(value: Any) -> bool:
    """A coded value as the library reads it: lower_snake_case words."""
    return (
        isinstance(value, str) and value.replace("_", "").isalpha() and value.islower()
    )


def _validation(entry: Mapping[str, Any]) -> dict:
    validation = entry.get("evidence", {}).get("scientific_validation")
    if not isinstance(validation, dict):
        raise _Refused("no_source_check", "the record has no scientific_validation")
    return validation


def ledger_row(
    row: Mapping[str, Any], before: Mapping[str, Any], after: Mapping[str, Any]
) -> dict:
    """``row`` with its copies of ``before`` (a ``scientific_validation``) replaced.

    ``usage_recommendation`` and ``recommendation_reason`` are copied from
    ``after`` (absent: removed). A changed ``recommendation_reason`` is taken off
    the end of ``qualification`` and the new one appended, as
    ``scripts/validate_primary_eos_refits.py`` composes the row. A row with
    sorted keys keeps them sorted.
    """
    result = copy.deepcopy(dict(row))
    for key in LEDGER_FIELDS:
        if row.get(key) != before.get(key):
            raise _Refused("ledger_differs", f"{LEDGER} {key} differs from the sidecar")
        if key in after:
            result[key] = copy.deepcopy(after[key])
        else:
            result.pop(key, None)
    old = before.get("recommendation_reason")
    new = after.get("recommendation_reason")
    if old != new:
        qualification = result.get("qualification", "")
        if old:
            if qualification == old:
                qualification = ""
            elif qualification.endswith(" " + old):
                qualification = qualification[: -len(old) - 1]
            else:
                raise _Refused(
                    "ledger_differs",
                    f"{LEDGER} qualification does not end with recommendation_reason",
                )
        if new:
            qualification = f"{qualification} {new}".strip()
        if qualification:
            result["qualification"] = qualification
        else:
            result.pop("qualification", None)
    if list(row) == sorted(row):
        result = dict(sorted(result.items()))
    return result


def _revise_recommendation(
    entry: dict,
    record: Mapping[str, Any],
    row: Mapping[str, Any] | None,
    value: str | None,
    stored: dict,
) -> dict:
    """Set the coded usage recommendation of a sidecar ``entry``; see the module doc.

    ``value`` ``None`` (or the stored original code) puts the stored originals
    back. Returns the record's new ledger row.
    """
    if row is None:
        raise _Refused("no_ledger_row", f"no row in {LEDGER}")
    assessment = record.get("evidence", {}).get("reproduction", {}).get("assessment")
    if not isinstance(assessment, dict) or assessment != assessment_for_record(
        record, row
    ):
        raise _Refused(
            "assessment_not_derived",
            "the stored assessment is not derived from the ledger row",
        )
    before = copy.deepcopy(_validation(entry))
    current = before.get(USAGE_KEY)
    original = stored.get(USAGE_KEY, current)
    target = original if value is None else value
    if target == original and USAGE_KEY in stored:
        # Back to the original: every bound recommendation field returns.
        for field in RECOMMENDATION_FIELDS:
            if field in stored:
                set_check_text(entry, field, stored.pop(field))
    else:
        stored.setdefault(USAGE_KEY, current)
        if isinstance(current, str) and not _is_code(current):
            # The replaced free text is kept verbatim.
            reason = before.get("recommendation_reason")
            if reason is None:
                stored.setdefault("recommendation_reason", None)
                set_check_text(entry, "recommendation_reason", current)
            elif reason != current:
                raise _Refused(
                    "recommendation_reason_present",
                    "the record has a different recommendation_reason; merge the "
                    "free-text usage_recommendation into it first",
                )
        set_check_text(entry, USAGE_KEY, target)
    after = _validation(entry)
    new_row = ledger_row(row, before, after)
    back = ledger_row(new_row, after, before)
    if list(back.items()) != list(row.items()):
        raise _Refused("not_reversible", f"the {LEDGER} row would not restore exactly")
    trial = {**record, "evidence": entry["evidence"]}
    derived = assessment_for_record(trial, new_row)
    changed = sorted(
        key
        for key in set(derived) | set(assessment)
        if key != "qualifications" and derived.get(key) != assessment.get(key)
    )
    if changed:
        raise _Refused(
            "beyond_fields", f"the assessment's {', '.join(changed)} would change"
        )
    qualifications = derived.get("qualifications")
    if USAGE_KEY in stored:
        stored.setdefault(
            "assessment_qualifications", check_text(entry, "assessment_qualifications")
        )
        set_check_text(entry, "assessment_qualifications", qualifications)
    elif qualifications != check_text(entry, "assessment_qualifications"):
        raise _Refused(
            "not_reversible", "the restored assessment qualifications differ"
        )
    for field in RECOMMENDATION_FIELDS:
        if field in stored and stored[field] == check_text(entry, field):
            del stored[field]
    return new_row


def revise_entry(
    entry: dict,
    record: Mapping[str, Any],
    row: Mapping[str, Any] | None,
    new: Mapping[str, Any],
    stored: dict,
) -> dict | None:
    """Apply one record's change to its sidecar ``entry`` in place.

    ``record`` is the merged record before the change, ``row`` its ledger row
    and ``stored`` its ``check_text`` originals (updated in place). Returns the
    new ledger row when the usage recommendation changes, else ``None``.
    Raises :class:`_Refused`.
    """
    new_row = None
    if USAGE_KEY in new:
        new_row = _revise_recommendation(entry, record, row, new[USAGE_KEY], stored)
    for key, field in INPUT_FIELDS.items():
        if key not in new:
            continue
        value = new[key]
        stored.setdefault(field, check_text(entry, field))
        set_check_text(entry, field, value)
        if stored[field] == value:
            del stored[field]
    return new_row


def _revise_material(
    material: str,
    document: dict,
    sidecar: dict,
    changes: Mapping[str, Mapping[str, Any]],
    originals: dict[int, dict],
    rows: Mapping[str, Mapping[str, Any]],
) -> tuple[dict[str, dict], dict[str, Any]]:
    """Apply ``changes`` to a material's files and originals, in place.

    Regenerates each changed record's ``check`` block. Returns the new ledger
    rows and each changed record's resulting ``check.warning``. Raises
    :class:`_Refused` for a change the files cannot take.
    """
    before = {r["identifier"]: r for r in merge_audit(document, sidecar)["eos_records"]}
    order = [r["identifier"] for r in document["eos_records"]]
    new_rows: dict[str, dict] = {}
    for identifier, new in changes.items():
        entry = sidecar["records"].get(identifier)
        if entry is None:
            raise _Refused("no_sidecar_entry", f"{material}/{identifier}")
        index = order.index(identifier)
        stored = originals.setdefault(index, {})
        row = revise_entry(entry, before[identifier], rows.get(identifier), new, stored)
        if row is not None:
            new_rows[identifier] = row
        if not stored:
            del originals[index]
    merged = merge_audit(document, sidecar)
    complete = {r["identifier"]: r for r in merged["eos_records"]}
    warnings: dict[str, Any] = {}
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
        kept = {"source_note", "warning_reason"}
        if USAGE_KEY in changes[identifier]:
            kept.add("warning")
        if {k: v for k, v in old.items() if k not in kept} != {
            k: v for k, v in derived.items() if k not in kept
        }:
            raise _Refused(
                "beyond_fields",
                f"the check block differs from its evidence beyond {sorted(kept)}",
            )
        record["check"] = derived
        warnings[identifier] = derived["warning"]
    return new_rows, warnings


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
    texts = [
        (key, value)
        for key, value in new.items()
        if key in INPUT_FIELDS and value is not None
    ]
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
    bindings = _receipt(root).get("materials", {})
    rows: dict[str, dict] | None = None

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
        originals = _originals(bindings.get(material, {}))
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
            stored = originals.get(order.index(identifier), {})
            usage = requested[identifier].get(USAGE_KEY, _ABSENT)
            if usage is not _ABSENT:
                if usage is not None and usage not in CODED_RECOMMENDATIONS:
                    reasons.append(
                        (
                            "not_coded",
                            "usage_recommendation must be null or one of "
                            + ", ".join(CODED_RECOMMENDATIONS),
                        )
                    )
                else:
                    current = check_text(record, USAGE_KEY)
                    target = usage if usage is not None else stored.get(USAGE_KEY)
                    if usage is None and USAGE_KEY not in stored:
                        target = current
                    if target != current:
                        new[USAGE_KEY] = usage
            if new.get("source_note") is not None and not isinstance(
                record.get("evidence", {})
                .get("scientific_validation", {})
                .get("primary_source_check"),
                dict,
            ):
                reasons.append(("no_source_check", "the record has no source check"))
            warning = usage_warning(record, merged)
            warning_after = warning
            if not reasons and USAGE_KEY in new:
                if rows is None:
                    rows = _ledger_rows(_ledger(root))
                try:
                    _, warnings = _revise_material(
                        material,
                        copy.deepcopy(document),
                        copy.deepcopy(sidecar),
                        {identifier: new},
                        copy.deepcopy(originals),
                        rows,
                    )
                    warning_after = warnings[identifier]
                except _Refused as error:
                    reasons.append((error.reason, error.detail))
            reason_after = new.get("warning_reason", warning_reason(record))
            if not reasons and reason_after is not None and warning_after is None:
                reasons.append(
                    (
                        "no_warning",
                        "the record has no warning (check.warning is null); a "
                        "warning_reason is only given with a warning",
                    )
                )
            if (
                not reasons
                and warning is None
                and warning_after is not None
                and reason_after is None
            ):
                reasons.append(
                    (
                        "reason_required",
                        "this gives the record its first warning; give a "
                        "warning_reason with it",
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


def build_changes(
    plan: Plan, root: Path = ROOT
) -> tuple[dict[str, tuple[bytes, bytes]], dict | None, dict[Path, str], dict]:
    """Compute and verify every file a plan writes, without writing.

    Returns the new record and sidecar bytes per material, the new receipt
    (``None``: no split material changes), the ledger files by path and each
    changed record's resulting ``check.warning`` by material and record.
    Raises ``ValueError`` if anything would not verify, including a pre-split
    file that would no longer be recovered exactly.
    """
    receipt = _receipt(root)
    bindings = receipt.get("materials", {})
    ledger = _ledger(root)
    rows = _ledger_rows(ledger)
    revised: dict[str, list] = {}
    files: dict[str, tuple[bytes, bytes]] = {}
    new_rows: dict[str, dict] = {}
    warnings: dict[str, dict[str, Any]] = {}
    for material, changes in plan.changes.items():
        document, sidecar = _read(root, material)
        originals = _originals(bindings.get(material, {}))
        try:
            material_rows, warnings[material] = _revise_material(
                material, document, sidecar, changes, originals, rows
            )
        except _Refused as error:
            raise ValueError(f"{material}: {error}") from None
        new_rows.update(material_rows)
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
    ledger_texts: dict[Path, str] = {}
    if new_rows and ledger is not None:
        ledger["records"] = [
            new_rows.get(row["record_identifier"], row) for row in ledger["records"]
        ]
        from scripts.validate_primary_eos_refits import render_markdown

        ledger_texts = {
            root / LEDGER: json.dumps(ledger, indent=2, ensure_ascii=False) + "\n",
            root / LEDGER_MARKDOWN: render_markdown(ledger),
        }
    return files, new_receipt, ledger_texts, warnings


def apply_plan(plan: Plan, root: Path = ROOT) -> dict[str, dict[str, Any]]:
    """Write the planned changes, regenerate ``check``, bind and verify.

    Every material's new files, the receipt binding and the ledger copies are
    computed and verified in memory first (:func:`build_changes`); files are
    written only when the whole plan is valid, so a failure leaves the tree
    untouched. Returns each changed record's resulting ``check.warning``.
    """
    files, new_receipt, ledger_texts, warnings = build_changes(plan, root)
    for material, (record_bytes, audit_bytes) in files.items():
        record_path, audit_path = _paths(root, material)
        record_path.write_bytes(record_bytes)
        audit_path.write_bytes(audit_bytes)
    if new_receipt is not None:
        (root / RECEIPT.relative_to(ROOT)).write_text(
            receipt_text(new_receipt), encoding="utf-8"
        )
    for path, text in ledger_texts.items():
        path.write_text(text, encoding="utf-8")
    return warnings


def verify(root: Path = ROOT) -> tuple[int, list[str]]:
    """Every ``check`` equals its derivation; curated texts recover losslessly.

    Also compares the ledger's ``usage_recommendation`` and
    ``recommendation_reason`` copies with the sidecars. Returns the number of
    records holding a curated text and the failures.
    """
    from scripts.catalog_revisions import check_bindings

    failures: list[str] = check_bindings(root)
    curated = 0
    bound = _receipt(root).get("materials", {})
    rows = _ledger_rows(_ledger(root))
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
            row = rows.get(identifier)
            if row is not None:
                validation = complete.get("evidence", {}).get(
                    "scientific_validation", {}
                )
                for key in LEDGER_FIELDS:
                    if row.get(key) != validation.get(key):
                        failures.append(
                            f"{material}/{identifier}: {key} differs in {LEDGER}"
                        )
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


def _describe(key: str, value: Any) -> str:
    if key == USAGE_KEY:
        return f"{key} {'original' if value is None else value}"
    return f"{key} {'removed' if value is None else f'{len(value)} characters'}"


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
    try:
        warnings = build_changes(plan)[3] if arguments.check else apply_plan(plan)
    except ValueError as error:
        print(f"FAILED {error}; nothing was written.")
        return 1
    for material, changes in plan.changes.items():
        for identifier, new in changes.items():
            parts = ", ".join(_describe(key, value) for key, value in new.items())
            warning = json.dumps(warnings[material][identifier], ensure_ascii=False)
            print(f"{tag} {material}/{identifier}: {parts}; warning: {warning}")
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
                    "warnings": warnings,
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
