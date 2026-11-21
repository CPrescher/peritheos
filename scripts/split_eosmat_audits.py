"""Repository maintenance: split the bundled v5 catalog into records and audits.

Each ``peritheos/data/materials/<id>.eosmat`` keeps the material header, the
``default_records`` map, the EOS records with a short ``check`` block, and
datasets. The audit trail moves to ``peritheos/data/audits/<id>.audit.json``.
The receipt ``docs/data/eosmat-v5-audit-split.json`` records, per material,
every change between the pre-split file and ``record file + sidecar`` merged:
the added ``check`` blocks, the per-record default flags replaced by
``default_records``, the removed generic statements and the key-order layout
of the pre-split file. ``pre_split_bytes`` uses it to recover each pre-split
file byte for byte.

The whole-catalog split ran once (``python scripts/split_eosmat_audits.py``)
and refuses to run again over an existing receipt. Afterwards, curators who
add a new bundled material as a complete inline-evidence v5 file run
``python scripts/split_eosmat_audits.py --add [MATERIAL ...]``. It splits
every named material, or every bundled file that still holds inline audit
blocks or has no sidecar, with the same mechanical procedure and adds its
binding to the receipt. It refuses a material that already has a binding
(the one earlier replacement, of pyrope, is kept under ``superseded``).
``--verify`` re-checks the current catalog. This is not an installed-package
tool; standalone EOSMAT files keep inline evidence and need no sidecar.

A later text revision (``scripts/apply_record_text.py``) replaces free text
(``label``, ``notes``, ``validity.notes``, ``thermal.notes``) and keeps every
replaced original in the record's ``superseded_record_text`` sidecar block. Its
binding lists those records under ``superseded_record_text``; recovery first
puts the stored originals back and removes the block (:func:`restore`), then
undoes the split as before.

The two curated ``check`` texts (``scripts/apply_check_text.py``) live in the
sidecar evidence: ``scientific_validation.primary_source_check.record_note``
and ``scientific_validation.warning_reason``. The same tool sets a coded
``scientific_validation.usage_recommendation`` (the source of
``check.warning``); with it may change ``scientific_validation.
recommendation_reason`` (the replaced prose, kept verbatim) and the
assessment's ``qualifications`` (``assessment_qualifications``). Written after
the split, their binding lists ``check_text`` entries ``[index, {field:
original}]`` holding each field's value before the first post-split change
(``null``: absent). Recovery puts those originals back before anything else;
``finding`` itself is never changed.

A change of scientific content after the split (an added EOS record, set or
corrected structured fields) goes through ``scripts/apply_catalog_revision.py``
and the ledger ``docs/data/catalog-revisions.json``. Re-splitting with
``--add`` would move the binding to ``superseded``, where nothing reads it, so
the original pre-split file and the bound text originals would no longer be
recovered; ``--add`` therefore refuses every material that has a binding. A
catalog revision keeps the binding; an added record's index joins
``added_check``, ``catalog_revision_sha256`` names the revised pre-split
reference that restoring the current files gives and
``catalog_revision_content_sha256`` pins the latest revision's content.
``--verify`` walks that reference back to ``before_sha256``
(``scripts/catalog_revisions.py``).
"""

from __future__ import annotations

import argparse
import copy
import hashlib
import json
import sys
from collections.abc import Callable, Iterable, Mapping
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from peritheos._eosmat_audit import (  # noqa: E402
    AUDIT_RECORD_KEYS,
    GENERIC_QUALIFICATIONS,
    GENERIC_VALIDATION_NOTES,
    SUPERSEDED_TEXT_FIELDS,
    audit_is_external,
    merge_audit,
    source_note,
    split_document,
    usage_warning,
    warning_reason,
)
from peritheos._eosmat_format import (  # noqa: E402
    V4_LAYOUT,
    V5_LAYOUT,
    Layout,
    format_audit,
    format_eosmat,
)
from peritheos.eosmat import validate_eosmat_document  # noqa: E402

MATERIALS = ROOT / "peritheos/data/materials"
AUDITS = ROOT / "peritheos/data/audits"
RECEIPT = ROOT / "docs/data/eosmat-v5-audit-split.json"
SOURCE_COMMIT = "4d89f0fa"


def _sha(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


GENERIC = {
    "qualifications": list(GENERIC_QUALIFICATIONS),
    "note": list(GENERIC_VALIDATION_NOTES),
}
_NOTE = ("evidence", "scientific_validation", "note")
_QUALIFICATIONS = ("evidence", "reproduction", "assessment", "qualifications")


#: Key-order layouts a pre-split file may use. Files split before the v5
#: reading-order revision used the v4 layout; ``layout`` names it in a binding.
LAYOUTS: dict[str | None, Layout] = {None: V5_LAYOUT, "v4": V4_LAYOUT}
_DEFAULT_FLAGS = ("default", "default_for")


def _render(document: dict, binding: dict) -> bytes:
    """The pre-split file bytes, in the layout the binding records."""
    return format_eosmat(document, layout=LAYOUTS[binding.get("layout")]).encode(
        "utf-8"
    )


def _node(document: Any, path: list) -> Any:
    for key in path:
        document = document[key]
    return document


def record_text(record: dict, field: str) -> Any:
    """Current value of one free-text field (``None`` when absent).

    ``field`` is one of :data:`SUPERSEDED_TEXT_FIELDS`; ``validity_notes`` and
    ``thermal_notes`` name ``validity.notes`` and ``thermal.notes``.
    """
    if field in ("label", "notes"):
        return copy.deepcopy(record.get(field))
    container = record.get("validity" if field == "validity_notes" else "thermal")
    if not isinstance(container, dict):
        return None
    return copy.deepcopy(container.get("notes"))


def set_record_text(record: dict, field: str, value: Any) -> None:
    """Write one free-text field; ``None``, ``""`` or ``[]`` removes it.

    A note keeps its position; a new note is appended (``notes`` is the last
    key of every bundled ``validity`` and ``thermal`` object). A ``validity``
    object left empty is removed; a missing one is created.
    """
    if field not in SUPERSEDED_TEXT_FIELDS:
        raise ValueError(f"Unknown record text field {field!r}")
    absent = value is None or value == "" or value == []
    if field in ("label", "notes"):
        if absent:
            record.pop(field, None)
        else:
            record[field] = copy.deepcopy(value)
        return
    key = "validity" if field == "validity_notes" else "thermal"
    container = record.get(key)
    if absent:
        if isinstance(container, dict):
            container.pop("notes", None)
            if key == "validity" and not container:
                del record[key]
        return
    if container is None and key == "validity":
        container = record["validity"] = {}
    if not isinstance(container, dict):
        raise ValueError(f"Record {record.get('identifier')!r} has no {key} object")
    container["notes"] = copy.deepcopy(value)


def restore_record_text(record: dict) -> None:
    """Put a record's ``superseded_record_text`` originals back; drop the block."""
    block = record.pop("superseded_record_text")
    for field in SUPERSEDED_TEXT_FIELDS:
        if field in block:
            set_record_text(record, field, block[field])


#: Curated sidecar texts copied into ``check`` (``source_note`` and
#: ``warning_reason``): ``record_note`` sits in
#: ``scientific_validation.primary_source_check`` directly after ``finding``,
#: ``warning_reason`` in ``scientific_validation`` directly after
#: ``usage_recommendation`` (or last).
CHECK_TEXT_FIELDS = ("record_note", "warning_reason")
#: A coded usage recommendation (the source of ``check.warning``) and what is
#: derived from it, bound the same way: ``usage_recommendation`` and
#: ``recommendation_reason`` in ``scientific_validation`` (a new
#: ``usage_recommendation`` last, a new ``recommendation_reason`` directly after
#: ``warning_reason`` or ``usage_recommendation``), and
#: ``assessment_qualifications``, the list
#: ``reproduction.assessment.qualifications`` (new: directly after ``basis``).
RECOMMENDATION_FIELDS = (
    "usage_recommendation",
    "recommendation_reason",
    "assessment_qualifications",
)
#: Every field a ``check_text`` binding entry may hold, in restore order.
CHECK_TEXT_BINDING_FIELDS = CHECK_TEXT_FIELDS + RECOMMENDATION_FIELDS
_CHECK_TEXT_AFTER: dict[str, tuple[str, ...]] = {
    "record_note": ("finding",),
    "warning_reason": ("usage_recommendation",),
    "usage_recommendation": (),
    "recommendation_reason": ("warning_reason", "usage_recommendation"),
    "assessment_qualifications": ("basis",),
}
_CHECK_TEXT_KEY = {"assessment_qualifications": "qualifications"}


def _check_text_container(record: dict, field: str, create: bool) -> dict | None:
    """The object holding ``field`` in a native or flattened v3 record."""
    evidence = record.get("evidence", record)
    if not isinstance(evidence, dict):
        return None
    if field == "assessment_qualifications":
        reproduction = evidence.get("reproduction")
        assessment = (
            reproduction.get("assessment") if isinstance(reproduction, dict) else None
        )
        if not isinstance(assessment, dict):
            if not create:
                return None
            raise ValueError(f"Record {record.get('identifier')!r} has no assessment")
        return assessment
    validation = evidence.get("scientific_validation")
    if not isinstance(validation, dict):
        if not create:
            return None
        raise ValueError(f"Record {record.get('identifier')!r} has no source check")
    if field != "record_note":
        return validation
    source = validation.get("primary_source_check")
    if not isinstance(source, dict):
        if not create:
            return None
        raise ValueError(
            f"Record {record.get('identifier')!r} has no primary_source_check"
        )
    return source


def check_text(record: dict, field: str) -> Any:
    """Current value of one bound check field (``None`` when absent).

    ``field`` is one of :data:`CHECK_TEXT_BINDING_FIELDS`.
    """
    if field not in CHECK_TEXT_BINDING_FIELDS:
        raise ValueError(f"Unknown check text field {field!r}")
    container = _check_text_container(record, field, create=False)
    if container is None:
        return None
    return copy.deepcopy(container.get(_CHECK_TEXT_KEY.get(field, field)))


def set_check_text(record: dict, field: str, value: Any) -> None:
    """Write one bound check field into a record's evidence; ``None`` removes it.

    A new value is placed directly after ``finding`` (``record_note``),
    ``usage_recommendation`` (``warning_reason``), ``warning_reason`` or else
    ``usage_recommendation`` (``recommendation_reason``) or ``basis``
    (``assessment_qualifications``), or last; an existing one keeps its
    position. ``finding`` is never touched. A stored ``check`` block is not
    changed here (see :func:`refresh_check_text`).
    """
    if field not in CHECK_TEXT_BINDING_FIELDS:
        raise ValueError(f"Unknown check text field {field!r}")
    container = _check_text_container(record, field, create=value is not None)
    if container is None:
        return
    key = _CHECK_TEXT_KEY.get(field, field)
    if value is None:
        container.pop(key, None)
        return
    value = copy.deepcopy(value)
    if key in container:
        container[key] = value
        return
    items = list(container.items())
    keys = [name for name, _ in items]
    after = next((a for a in _CHECK_TEXT_AFTER[field] if a in keys), None)
    position = len(items) if after is None else keys.index(after) + 1
    items.insert(position, (key, value))
    container.clear()
    container.update(items)


def refresh_check_text(record: dict, material: dict | None = None) -> None:
    """Set a stored ``check`` block's ``source_note`` and ``warning_reason``.

    Both are copied from the record's evidence exactly as
    :func:`peritheos._eosmat_audit.check_block` derives them; given the
    ``material`` (for a ``recommendation_hold``), ``warning`` is derived too.
    A record without a ``check`` block is left alone.
    """
    check = record.get("check")
    if not isinstance(check, dict):
        return
    if material is not None:
        check["warning"] = usage_warning(record, material)
    for key, value in (
        ("source_note", source_note(record)),
        ("warning_reason", warning_reason(record)),
    ):
        if value is None:
            check.pop(key, None)
        else:
            check[key] = value


def restore_check_text(
    record: dict, originals: dict, material: dict | None = None
) -> None:
    """Put the pre-revision bound check fields back and refresh ``check``."""
    for field in CHECK_TEXT_BINDING_FIELDS:
        if field in originals:
            set_check_text(record, field, originals[field])
    refresh_check_text(record, material)


def restore(document: dict, binding: dict) -> dict:
    """Undo the split changes listed in one material's receipt binding.

    ``superseded_record_text`` names records whose free text was revised after
    the split; their stored originals are restored and the block removed
    before anything else. ``check_text`` entries ``[index, {field:
    original}]`` then restore the bound check fields written after the split
    (curated check texts and a coded usage recommendation with its derived
    ``recommendation_reason`` and assessment qualifications).
    ``added_check`` names records whose ``check`` block was added.
    ``added_default_records`` marks a ``default_records`` map made from the
    per-record flags; ``default_flags`` restores each record's original
    ``default`` / ``default_for`` keys as ``[index, {key: value}]``.
    ``generic_notes`` / ``generic_qualifications`` restore removed generic
    statements by index into the receipt's ``removed_generic_statements``
    (qualifications as ``[position, statement]`` pairs in original order).
    ``object_orders`` restore key order where it affects the rendered file.
    """
    result = copy.deepcopy(document)
    records = result["eos_records"]
    for index in binding.get("superseded_record_text", []):
        restore_record_text(records[index])
    for index, originals in binding.get("check_text", []):
        restore_check_text(records[index], originals, result)
    for index in binding.get("added_check", []):
        del records[index]["check"]
    if binding.get("added_default_records"):
        del result["default_records"]
    for index, flags in binding.get("default_flags", []):
        records[index].update(flags)
    for index, statement in binding.get("generic_notes", []):
        validation = records[index]["evidence"]["scientific_validation"]
        validation["note"] = GENERIC["note"][statement]
    for index, inserts in binding.get("generic_qualifications", []):
        assessment = records[index]["evidence"]["reproduction"]["assessment"]
        values = assessment.setdefault("qualifications", [])
        for position, statement in inserts:
            values.insert(position, GENERIC["qualifications"][statement])
    for ordering in binding.get("object_orders", []):
        node = _node(result, ordering["path"])
        values = {key: node[key] for key in ordering["keys"]}
        node.clear()
        node.update(values)
    return result


def _differences(before: dict, merged: dict) -> dict:
    """Describe the split as explicit, allowed restore operations."""
    binding: dict[str, Any] = {
        "added_check": [],
        "generic_notes": [],
        "generic_qualifications": [],
    }
    if "default_records" in merged and "default_records" not in before:
        binding["added_default_records"] = True
    flags: list = []
    if len(before["eos_records"]) != len(merged["eos_records"]):
        raise ValueError("Record count changed")
    for index, (old, new) in enumerate(
        zip(before["eos_records"], merged["eos_records"])
    ):
        if "check" in new and "check" not in old:
            binding["added_check"].append(index)
        removed = {k: old[k] for k in _DEFAULT_FLAGS if k in old and k not in new}
        if removed:
            flags.append([index, removed])
        try:
            note = _node(old, list(_NOTE))
        except (KeyError, TypeError):
            note = None
        if note in GENERIC["note"]:
            binding["generic_notes"].append([index, GENERIC["note"].index(note)])
        try:
            qualifications = _node(old, list(_QUALIFICATIONS))
        except (KeyError, TypeError):
            qualifications = []
        inserts = [
            [position, GENERIC["qualifications"].index(text)]
            for position, text in enumerate(qualifications)
            if text in GENERIC["qualifications"]
        ]
        if inserts:
            binding["generic_qualifications"].append([index, inserts])
    if flags:
        binding["default_flags"] = flags
    return binding


def _order_candidates(before: Any, after: Any, path: tuple = ()) -> list[dict]:
    result: list[dict] = []
    if isinstance(before, dict) and isinstance(after, dict):
        if list(before) != list(after) and set(before) == set(after):
            result.append({"path": list(path), "keys": list(before)})
        for key in before:
            if key in after:
                result.extend(_order_candidates(before[key], after[key], (*path, key)))
    elif isinstance(before, list) and isinstance(after, list):
        for index, (a, b) in enumerate(zip(before, after)):
            result.extend(_order_candidates(a, b, (*path, index)))
    return result


def merged_document(material: str, root: Path = ROOT) -> dict:
    """Current record file plus its sidecar, as one inline v5 document."""
    record_file = json.loads(
        (root / f"peritheos/data/materials/{material}.eosmat").read_bytes()
    )
    sidecar = json.loads(
        (root / f"peritheos/data/audits/{material}.audit.json").read_bytes()
    )
    return merge_audit(record_file, sidecar)


def pre_split_bytes(material: str, *, root: Path | None = None) -> bytes:
    """Recover the exact pre-split v5 file, verifying every bound hash.

    For a material with catalog revisions this is the pre-split reference
    after its latest revision (``catalog_revision_sha256``);
    :func:`scripts.catalog_revisions.prior_revision_bytes` recovers the
    original file from it.
    """
    root = ROOT if root is None else root
    relative = f"peritheos/data/materials/{material}.eosmat"
    data = (root / relative).read_bytes()
    receipt_path = root / RECEIPT.relative_to(ROOT)
    if not receipt_path.exists():
        return data
    binding = json.loads(receipt_path.read_bytes())["materials"].get(material)
    if binding is None:
        return data
    from scripts.catalog_revisions import check_head

    # A binding naming a catalog revision needs exactly that ledger head.
    check_head(material, binding, root)
    audit = (root / f"peritheos/data/audits/{material}.audit.json").read_bytes()
    return recover_pre_split(material, data, audit, binding)


def recover_pre_split(
    material: str, record_bytes: bytes, audit_bytes: bytes, binding: dict
) -> bytes:
    """Recover a pre-split file from given file bytes and their binding."""
    relative = f"peritheos/data/materials/{material}.eosmat"
    if _sha(record_bytes) != binding["after_record_sha256"]:
        raise ValueError(f"EOSMAT audit split record file changed: {relative}")
    if _sha(audit_bytes) != binding["after_audit_sha256"]:
        raise ValueError(f"EOSMAT audit split sidecar changed: {material}")
    merged = merge_audit(json.loads(record_bytes), json.loads(audit_bytes))
    original = _render(restore(merged, binding), binding)
    # After a catalog revision the files restore to the revised reference;
    # scripts/catalog_revisions.py walks it back to before_sha256.
    expected = binding.get("catalog_revision_sha256", binding["before_sha256"])
    if _sha(original) != expected:
        raise ValueError(f"EOSMAT audit split predecessor changed: {relative}")
    return original


def pre_split_view(document: dict, root: Path = ROOT) -> dict:
    """Project a complete bundled document onto its pre-split content.

    For frozen historical fingerprints only. Accepts the native or the
    v3-shaped layout, matches records by identifier, removes ``check``,
    restores the original per-record default flags in place of
    ``default_records`` and re-inserts the generic statements recorded in the
    receipt. Text revised after the split is restored from the
    ``superseded_record_text`` blocks and ``check_text`` originals the
    receipt binds. Documents without a
    receipt binding are returned unchanged.
    """
    receipt_path = root / RECEIPT.relative_to(ROOT)
    material = document.get("identifier")
    if not receipt_path.exists() or not isinstance(material, str):
        return copy.deepcopy(document)
    binding = json.loads(receipt_path.read_bytes())["materials"].get(material)
    record_path = root / f"peritheos/data/materials/{material}.eosmat"
    if binding is None or not record_path.exists():
        return copy.deepcopy(document)
    order = [
        r["identifier"] for r in json.loads(record_path.read_bytes())["eos_records"]
    ]
    result = copy.deepcopy(document)
    records = {r.get("identifier"): r for r in result.get("eos_records", [])}

    def evidence(index: int) -> dict | None:
        record = records.get(order[index])
        if record is None:
            return None
        return record.get("evidence", record)

    for index in binding.get("superseded_record_text", []):
        record = records.get(order[index])
        if record is not None and "superseded_record_text" in record:
            restore_record_text(record)
    for index, originals in binding.get("check_text", []):
        record = records.get(order[index])
        if record is not None:
            restore_check_text(record, originals)
    for record in records.values():
        record.pop("check", None)
    # Bindings record every pre-split default flag; exports and the
    # default_records map spell the same selection differently.
    if "added_check" in binding:
        result.pop("default_records", None)
        for record in records.values():
            for key in _DEFAULT_FLAGS:
                record.pop(key, None)
        for index, flags in binding.get("default_flags", []):
            record = records.get(order[index])
            if record is not None:
                record.update(flags)
    for index, statement in binding.get("generic_notes", []):
        node = evidence(index)
        if node is not None and "note" not in node.get("scientific_validation", {}):
            node["scientific_validation"]["note"] = GENERIC["note"][statement]
    for index, inserts in binding.get("generic_qualifications", []):
        node = evidence(index)
        if node is None:
            continue
        assessment = node["reproduction"]["assessment"]
        values = assessment.setdefault("qualifications", [])
        for position, statement in inserts:
            text = GENERIC["qualifications"][statement]
            if text not in values:
                values.insert(position, text)
    return result


def split_material(original: bytes) -> tuple[bytes, bytes, dict]:
    """Split one complete inline v5 file; return record, sidecar and binding.

    The input must be canonically formatted in the current v5 layout or, for
    files written before the v5 reading-order revision, the v4 layout; the
    binding names the layout so the file is recovered byte for byte.
    """
    before = json.loads(original)
    validate_eosmat_document(before)
    layout = next(
        (
            name
            for name, candidate in LAYOUTS.items()
            if format_eosmat(before, layout=candidate).encode("utf-8") == original
        ),
        False,
    )
    if layout is False:
        raise ValueError("Pre-split file is not in canonical EOSMAT layout")
    record_file, sidecar = split_document(before)
    record_bytes = format_eosmat(record_file).encode("utf-8")
    audit_bytes = format_audit(sidecar).encode("utf-8")
    # Reload the rendered files so the proof covers the bytes actually written.
    merged = merge_audit(json.loads(record_bytes), json.loads(audit_bytes))
    validate_eosmat_document(json.loads(record_bytes))
    validate_eosmat_document(merged)
    binding: dict[str, Any] = {
        "before_sha256": _sha(original),
        "after_record_sha256": _sha(record_bytes),
        "after_audit_sha256": _sha(audit_bytes),
        **({"layout": layout} if layout is not None else {}),
        **_differences(before, merged),
    }
    restored = restore(merged, binding)
    if restored != before:
        raise ValueError("Unexpected audit split changes beyond the receipt")
    # Keep only the key orders that the rendered file actually depends on.
    binding["object_orders"] = _order_candidates(before, restored)
    for ordering in list(binding["object_orders"]):
        trial = {
            **binding,
            "object_orders": [o for o in binding["object_orders"] if o is not ordering],
        }
        if _render(restore(merged, trial), binding) == original:
            binding = trial
    restored = restore(merged, binding)
    if restored != before or _render(restored, binding) != original:
        raise ValueError("Record file and sidecar do not reconstruct the source")
    return record_bytes, audit_bytes, binding


def split_catalog() -> None:
    from scripts.scientific_source_snapshot import scientific_sha256

    if RECEIPT.exists():
        raise ValueError("The audit split already has a receipt; use --verify")
    before_scientific = scientific_sha256()
    receipt: dict[str, Any] = {
        "schema": "peritheos.eosmat-v5-audit-split.v1",
        "source_commit": SOURCE_COMMIT,
        "scope": (
            "Layout only: per-record audit blocks moved to peritheos/data/audits; "
            "a check headline added; identifier-derived report links omitted from "
            "sidecars and regenerated on load; the generic statements listed here "
            "removed. Coefficients, errors, covariance, datasets, sources and "
            "scientific classifications are unchanged."
        ),
        "removed_generic_statements": {
            "evidence.reproduction.assessment.qualifications": list(
                GENERIC_QUALIFICATIONS
            ),
            "evidence.scientific_validation.note": list(GENERIC_VALIDATION_NOTES),
        },
        "before_scientific_sha256": before_scientific,
        "materials": {},
    }
    prepared = []
    for path in sorted(MATERIALS.glob("*.eosmat")):
        original = path.read_bytes()
        record_bytes, audit_bytes, binding = split_material(original)
        material = path.stem
        receipt["materials"][material] = binding
        prepared.append((path, original, record_bytes, audit_bytes))
    if any(path.read_bytes() != original for path, original, _, _ in prepared):
        raise ValueError("A source file changed while preparing the audit split")
    AUDITS.mkdir(exist_ok=True)
    init = AUDITS / "__init__.py"
    if not init.exists():
        init.write_text('"""Bundled EOSMAT v5 audit sidecars, one per material."""\n')
    for path, _, record_bytes, audit_bytes in prepared:
        path.write_bytes(record_bytes)
        (AUDITS / f"{path.stem}.audit.json").write_bytes(audit_bytes)
    after_scientific = scientific_sha256()
    if after_scientific != before_scientific:
        raise ValueError("The audit split changed the scientific snapshot")
    receipt["after_scientific_sha256"] = after_scientific
    RECEIPT.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
    verify_catalog()
    print(f"Split {len(prepared)} materials; every pre-split file recovers exactly.")


def verify_catalog(root: Path = ROOT) -> int:
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    records = sorted(
        p.stem for p in (root / "peritheos/data/materials").glob("*.eosmat")
    )
    audits = sorted(
        p.name.removesuffix(".audit.json")
        for p in (root / "peritheos/data/audits").glob("*.audit.json")
    )
    if records != audits:
        raise ValueError("Every material needs exactly one audit sidecar")
    from scripts.catalog_revisions import prior_revision_bytes

    for material, binding in receipt["materials"].items():
        # Revised materials are walked back to the original pre-split file.
        original = prior_revision_bytes(
            pre_split_bytes(material, root=root), material, root
        )
        if _sha(original) != binding["before_sha256"]:
            raise ValueError(f"EOSMAT audit split predecessor not reached: {material}")
    from scripts.catalog_revisions import check_bindings

    failures = check_bindings(root)
    if failures:
        raise ValueError("; ".join(failures))
    return len(receipt["materials"])


RECORD_TEXT_REVISION = (
    "Free-text revisions after the split (label, notes, validity.notes, "
    "thermal.notes) keep every replaced original verbatim in the record's "
    "superseded_record_text sidecar block. A binding's superseded_record_text "
    "lists those record indices; recovery restores the originals and removes "
    "the block first. after_record_sha256 and after_audit_sha256 bind the "
    "revised files; before_sha256 values are unchanged."
)


def receipt_text(receipt: dict) -> str:
    """The receipt file's text."""
    return json.dumps(receipt, indent=2, ensure_ascii=False) + "\n"


def _revised_receipt(
    receipt: dict,
    revised: Mapping[str, Any],
    files: Mapping[str, tuple[bytes, bytes]],
    update: Callable[[dict, Any], None],
    statement: tuple[str, str],
) -> dict | None:
    """A copy of ``receipt`` binding revised files, verified in memory.

    ``files`` maps each revised material to its new record and sidecar bytes.
    ``update`` writes one material's revision into its binding. Every changed
    binding must still recover its pre-split file exactly from those bytes;
    ``None`` means no split material was revised. Nothing is written.
    """
    result = copy.deepcopy(receipt)
    changed = []
    for material, revision in sorted(revised.items()):
        binding = result["materials"].get(material)
        if binding is None:
            continue
        update(binding, revision)
        record_bytes, audit_bytes = files[material]
        binding["after_record_sha256"] = _sha(record_bytes)
        binding["after_audit_sha256"] = _sha(audit_bytes)
        changed.append(material)
    if not changed:
        return None
    result.setdefault(*statement)
    for material in changed:
        recover_pre_split(material, *files[material], result["materials"][material])
    return result


def _disk_files(materials: Iterable[str], root: Path) -> dict[str, tuple[bytes, bytes]]:
    return {
        material: (
            (root / f"peritheos/data/materials/{material}.eosmat").read_bytes(),
            (root / f"peritheos/data/audits/{material}.audit.json").read_bytes(),
        )
        for material in materials
    }


def _set_record_text_indices(binding: dict, indices: list[int]) -> None:
    if indices:
        binding["superseded_record_text"] = sorted(set(indices))
    else:
        binding.pop("superseded_record_text", None)


def record_text_receipt(
    receipt: dict,
    revised: Mapping[str, list[int]],
    files: Mapping[str, tuple[bytes, bytes]],
) -> dict | None:
    """``receipt`` binding free-text revisions of ``files`` (see below)."""
    return _revised_receipt(
        receipt,
        revised,
        files,
        _set_record_text_indices,
        ("record_text_revision", RECORD_TEXT_REVISION),
    )


def bind_record_text_revision(revised: dict[str, list[int]], root: Path = ROOT) -> None:
    """Bind text revisions of split materials in the receipt and verify them.

    ``revised`` maps a material to the complete list of indices of records
    whose ``superseded_record_text`` block was added after the split (an empty
    list removes the entry). The after-hashes are refreshed from the written
    files, and the unchanged pre-split file must still be recovered exactly
    before the receipt is written.
    """
    receipt_path = root / RECEIPT.relative_to(ROOT)
    receipt = json.loads(receipt_path.read_bytes())
    files = _disk_files([m for m in revised if m in receipt["materials"]], root)
    result = record_text_receipt(receipt, revised, files)
    if result is not None:
        receipt_path.write_text(receipt_text(result), encoding="utf-8")


CHECK_TEXT_REVISION = (
    "Curated check texts written after the split "
    "(scientific_validation.primary_source_check.record_note and "
    "scientific_validation.warning_reason, copied into check.source_note and "
    "check.warning_reason) and coded usage recommendations "
    "(scientific_validation.usage_recommendation, the source of check.warning, "
    "with the recommendation_reason and assessment qualifications derived from "
    "it, field assessment_qualifications) are bound as check_text entries "
    "[index, {field: original}], original null meaning absent; recovery "
    "restores them before removing the check blocks. finding is never changed. "
    "after_record_sha256 and after_audit_sha256 bind the revised files; "
    "before_sha256 values are unchanged."
)


def _set_check_text_entries(binding: dict, entries: list) -> None:
    if entries:
        binding["check_text"] = sorted(
            [[index, dict(originals)] for index, originals in entries],
            key=lambda entry: entry[0],
        )
    else:
        binding.pop("check_text", None)


def check_text_receipt(
    receipt: dict,
    revised: Mapping[str, list],
    files: Mapping[str, tuple[bytes, bytes]],
) -> dict | None:
    """``receipt`` binding curated check-text revisions of ``files``.

    Computed and verified in memory (see :func:`bind_check_text_revision`);
    ``None`` when no split material is revised.
    """
    return _revised_receipt(
        receipt,
        revised,
        files,
        _set_check_text_entries,
        ("check_text_revision", CHECK_TEXT_REVISION),
    )


def bind_check_text_revision(revised: dict[str, list], root: Path = ROOT) -> None:
    """Bind curated check-text revisions of split materials in the receipt.

    ``revised`` maps a material to its complete ``check_text`` list of
    ``[index, {field: original}]`` (an empty list removes the entry). The
    after-hashes are refreshed from the written files, and the unchanged
    pre-split file must still be recovered exactly before the receipt is
    written.
    """
    receipt_path = root / RECEIPT.relative_to(ROOT)
    receipt = json.loads(receipt_path.read_bytes())
    files = _disk_files([m for m in revised if m in receipt["materials"]], root)
    result = check_text_receipt(receipt, revised, files)
    if result is not None:
        receipt_path.write_text(receipt_text(result), encoding="utf-8")


def holds_inline_audit(material: str, root: Path = ROOT) -> bool:
    """Whether a bundled file still needs the split: inline audit or no sidecar."""
    document = json.loads(
        (root / f"peritheos/data/materials/{material}.eosmat").read_bytes()
    )
    records = [r for r in document.get("eos_records", []) if isinstance(r, dict)]
    inline = [r["identifier"] for r in records if set(r) & set(AUDIT_RECORD_KEYS)]
    external = [r["identifier"] for r in records if audit_is_external(r)]
    if inline and external:
        raise ValueError(
            f"{material}: records {inline} hold inline audit blocks while "
            f"{external} point to a sidecar; write one complete inline file"
        )
    sidecar = root / f"peritheos/data/audits/{material}.audit.json"
    return bool(inline) or not sidecar.exists()


def split_added_materials(
    materials: list[str] | None = None,
    *,
    root: Path = ROOT,
    check_snapshot: bool = True,
) -> list[str]:
    """Split complete inline v5 files added after the catalog-wide split.

    Each file is treated exactly like the original catalog: it must be a
    valid, canonically formatted v5 document, its split must reconstruct it
    byte for byte, and the scientific snapshot must not change. Its receipt
    binding is added. A material that already has a binding is refused: change
    it with the text tools or ``scripts/apply_catalog_revision.py``. Without
    names, every bundled file holding inline audit blocks or lacking a sidecar
    is split; a named material must be such a file.
    """
    from scripts.scientific_source_snapshot import scientific_sha256

    receipt_path = root / RECEIPT.relative_to(ROOT)
    if not receipt_path.exists():
        raise ValueError("No audit-split receipt; run the catalog-wide split first")
    receipt = json.loads(receipt_path.read_bytes())
    material_dir = root / "peritheos/data/materials"
    audit_dir = root / "peritheos/data/audits"
    available = sorted(p.stem for p in material_dir.glob("*.eosmat"))
    if materials is None:
        selected = [m for m in available if holds_inline_audit(m, root)]
    else:
        unknown = sorted(set(materials) - set(available))
        if unknown:
            raise ValueError(f"No bundled material file for {unknown}")
        already = [m for m in materials if not holds_inline_audit(m, root)]
        if already:
            raise ValueError(
                f"{already} are already split; --add splits a new material "
                "written as one complete inline-evidence file"
            )
        selected = sorted(set(materials))
    if not selected:
        return []
    before_scientific = scientific_sha256(root) if check_snapshot else None
    prepared = []
    for material in selected:
        path = material_dir / f"{material}.eosmat"
        original = path.read_bytes()
        record_bytes, audit_bytes, binding = split_material(original)
        prepared.append((material, path, original, record_bytes, audit_bytes, binding))
    if any(path.read_bytes() != original for _, path, original, *_ in prepared):
        raise ValueError("A source file changed while preparing the audit split")
    existing = [m for m in selected if m in receipt["materials"]]
    if existing:
        # A new binding would supersede the existing one, whose original
        # pre-split file, bound text originals and catalog revisions would no
        # longer be recovered.
        raise ValueError(
            f"{existing} are bundled materials with an audit-split binding; --add "
            "is only for new materials. Change text with "
            "scripts/apply_record_text.py, check texts and warnings with "
            "scripts/apply_check_text.py, and records or fields with "
            "scripts/apply_catalog_revision.py"
        )
    for material, path, _, record_bytes, audit_bytes, binding in prepared:
        receipt["materials"][material] = binding
        path.write_bytes(record_bytes)
        (audit_dir / f"{material}.audit.json").write_bytes(audit_bytes)
    if check_snapshot and scientific_sha256(root) != before_scientific:
        raise ValueError("The audit split changed the scientific snapshot")
    receipt["materials"] = dict(sorted(receipt["materials"].items()))
    if "superseded" in receipt:
        receipt["superseded"] = dict(sorted(receipt["superseded"].items()))
    receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
    verify_catalog(root)
    return selected


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    action = parser.add_mutually_exclusive_group()
    action.add_argument("--verify", action="store_true")
    action.add_argument(
        "--add",
        nargs="*",
        metavar="MATERIAL",
        help="split added inline-evidence files (default: every such file)",
    )
    arguments = parser.parse_args()
    if arguments.verify:
        print(f"Verified {verify_catalog()} materials against the split receipt.")
    elif arguments.add is not None:
        try:
            added = split_added_materials(arguments.add or None)
        except ValueError as error:
            print(f"REFUSED {error}")
            sys.exit(1)
        print(f"Split {len(added)} added materials: {', '.join(added) or 'none'}.")
    else:
        split_catalog()


if __name__ == "__main__":
    main()


__all__ = [
    "CHECK_TEXT_BINDING_FIELDS",
    "CHECK_TEXT_FIELDS",
    "RECOMMENDATION_FIELDS",
    "bind_check_text_revision",
    "bind_record_text_revision",
    "check_text",
    "merged_document",
    "pre_split_bytes",
    "record_text",
    "refresh_check_text",
    "restore_check_text",
    "restore_record_text",
    "set_check_text",
    "set_record_text",
    "split_added_materials",
    "split_material",
]
