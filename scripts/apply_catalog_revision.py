"""Repository maintenance: change the scientific content of a bundled material.

One request file (JSON) describes one revision of one material, of one kind:

``add_record``
    append a complete inline EOS record (record fields and its audit
    ``evidence``, without ``check``). The record file gets the record and a
    ``check`` block derived from the evidence; the sidecar gets its audit
    entry. The material's ``default_records`` is unchanged.
``set_fields``
    set, correct or remove structured fields of existing records, each named
    by a path with its expected current value; the record's ``check`` block is
    regenerated. Identifier, label, notes, ``validity.notes``,
    ``thermal.notes`` (use ``scripts/apply_record_text.py``), ``check``,
    default flags and audit blocks are refused.

``remove_record`` and material-level fields (``default_records``) are not
supported. Request format::

    {"id": "2026-11-20-qandilite-lv-2016-confidence",
     "material": "qandilite", "kind": "set_fields", "date": "2026-11-20",
     "decided_by": "<name>", "reason": "<why>",
     "source": {"citation": "<reference>", "location": "<table, page>"},
     "rights": {"reviewed_by": "<name>", "date": "2026-11-20",
                "statement": "<what was reviewed>",
                "same_basis": true},             # add_record only
     "changes": [{"record": "<record id>", "path": ["parameter_error_confidence"],
                  "expected_absent": true,       # or "expected": <value>
                  "value": 0.6826894921370859}]} # or "remove": true
    # add_record: "record": {...complete inline record...} instead of "changes"

``rights`` is required for every revision because the rights review
(``docs/data/eos-structure-rights-policy.json``) pins the whole reviewed
material: an ``add_record`` statement says that the added record is covered by
the same source and redistribution basis as the material's existing records
(``same_basis: true``) or, with ``same_basis: false``, what differs
(``differences``); a ``set_fields`` statement says that the change adds no
third-party content. The tool never writes that statement itself. It stores
it in the revision, moves the material's pin to the revised pre-split
reference, lists the revision under the policy's ``catalog_revisions`` and
regenerates the review outputs with ``scripts/review_eos_structure_rights.py``.

Everything is computed and checked in a staging copy first: the files restore
to the revised pre-split reference through the audit-split binding, the
revision reverses to the previous reference and, through every earlier
revision, to the original pre-split file and the unchanged v5-migration
predecessor; the scientific snapshot changes by exactly the declared
operations; the reverse operations restore the complete document. Then the
files are written; if a check on the written tree fails, every file is put
back. ``--undo ID`` reverses the most recent revision of a material exactly
(record file, sidecar, receipt, ledger, rights policy and review outputs byte
for byte) and refuses an older one or one whose files changed since.

Usage::

    python scripts/apply_catalog_revision.py REQUEST.json --check  # dry run
    python scripts/apply_catalog_revision.py REQUEST.json          # apply
    python scripts/apply_catalog_revision.py --undo ID [--check]
    python scripts/apply_catalog_revision.py --verify

Afterwards regenerate the derived files the tool names (see
``REGENERATE``). This is not an installed-package tool.
"""

from __future__ import annotations

import argparse
import copy
import datetime
import hashlib
import json
import re
import shutil
import sys
import tempfile
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import dataclass, field
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
    bundled_audit_reference,
    check_block,
    compact_entry,
    expand_entry,
    merge_audit,
)
from peritheos._eosmat_format import format_audit, format_eosmat  # noqa: E402
from peritheos.eosmat import validate_eosmat_document  # noqa: E402
from scripts.catalog_revisions import (  # noqa: E402
    ABSENT,
    LEDGER,
    SPLIT_RECEIPT_STATEMENT,
    apply_operations,
    fingerprint,
    ledger_text,
    prior_revision_snapshot,
    read_ledger,
    read_path,
    reverse_operations,
    verify_ledger,
)
from scripts.check_record_text import record_cases  # noqa: E402
from scripts.review_eos_structure_rights import (  # noqa: E402
    GROUPS,
    INDEX,
    POLICY,
    REVIEW,
    build_review,
)
from scripts.scientific_source_snapshot import scientific_snapshot  # noqa: E402
from scripts.split_eosmat_audits import (  # noqa: E402
    LAYOUTS,
    RECEIPT,
    pre_split_bytes,
    receipt_text,
    record_text,
    recover_pre_split,
)

RECEIPT_PATH = RECEIPT.relative_to(ROOT).as_posix()
COMMON_KEYS = (
    "id",
    "material",
    "kind",
    "date",
    "decided_by",
    "reason",
    "source",
    "rights",
)
KIND_KEYS = {"add_record": "record", "set_fields": "changes"}
CHANGE_KEYS = ("record", "path", "expected", "expected_absent", "value", "remove")
#: Record fields a ``set_fields`` revision may not change: the identity, the
#: free text bound by ``apply_record_text.py``, the generated ``check``
#: block, default flags (the material's ``default_records`` decides) and the
#: audit blocks of the sidecar.
PROTECTED_FIELDS = (
    "identifier",
    "label",
    "notes",
    "check",
    "default",
    "default_for",
    *AUDIT_RECORD_KEYS,
)
PROTECTED_PATHS = (("validity", "notes"), ("thermal", "notes"))
_REVISION_ID = re.compile(r"[a-z0-9][a-z0-9._-]*\Z")
_RECORD_ID = re.compile(r"[a-z0-9][a-z0-9_]*\Z")
#: Generated files a revision can make stale, in the order to regenerate them.
#: The tool writes the record file, sidecar, receipts and rights review itself.
REGENERATE = (
    (
        "python scripts/validate_primary_eos_refits.py",
        "an added record's row in the registered refit ledger "
        "(docs/data/primary-eos-refits.json, docs/primary-eos-refits.md); it "
        "refits every record and rewrites the whole ledger, so review its diff",
    ),
    (
        "python -m scripts.generate_reconstruction_reports",
        "reconstruction report pages, plots and index totals",
    ),
    ("python -m scripts.generate_material_summaries", "material summaries"),
    ("python -m scripts.generate_catalog_coverage", "catalog coverage"),
)
#: Files no generator writes; an added record must be reflected by hand.
UPDATE_BY_HAND = (
    "peritheos/data/materials/manifest.json: eos_records and the "
    "scientific_validation and pressure_calibration counts",
    "peritheos/data/primary-source-audit.json: the record's entry and summary",
    "README.md: the catalog count badge",
    "tests that assert catalog counts or a material's record count",
)


class Refused(ValueError):
    """A request the catalog cannot take; nothing is written."""

    def __init__(self, reason: str, detail: str) -> None:
        super().__init__(f"{reason}: {detail}")
        self.reason = reason
        self.detail = detail


@dataclass
class Plan:
    """Files a revision (or its undo) writes, and what to show the operator."""

    revision: dict
    undo: bool
    files: dict[str, bytes] = field(default_factory=dict)
    shown: list[str] = field(default_factory=list)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _paths(material: str) -> tuple[str, str]:
    return (
        f"peritheos/data/materials/{material}.eosmat",
        f"peritheos/data/audits/{material}.audit.json",
    )


def _same(left: Any, right: Any) -> bool:
    if left is ABSENT or right is ABSENT:
        return left is right
    return json.dumps(left, sort_keys=True) == json.dumps(right, sort_keys=True)


def _text(value: Any) -> bool:
    return isinstance(value, str) and bool(value.strip())


def _date(value: Any, name: str) -> None:
    try:
        datetime.date.fromisoformat(value)
    except (TypeError, ValueError):
        raise Refused("request", f"{name} must be a YYYY-MM-DD date") from None


# Request --------------------------------------------------------------------


def check_request(request: Any) -> None:
    """Refuse a request whose shape is wrong, incomplete or carries extra keys."""
    if not isinstance(request, dict):
        raise Refused("request", "a request is one JSON object")
    kind = request.get("kind")
    if kind not in KIND_KEYS:
        raise Refused("request", f"kind must be one of {sorted(KIND_KEYS)}")
    allowed = {*COMMON_KEYS, KIND_KEYS[kind]}
    extra = sorted(set(request) - allowed)
    if extra:
        raise Refused("undeclared_change", f"{kind} does not take {extra}")
    missing = sorted(allowed - set(request))
    if missing:
        raise Refused("request", f"missing {missing}")
    if not isinstance(request["id"], str) or not _REVISION_ID.match(request["id"]):
        raise Refused("request", "id must be lower-case letters, digits, '.', '-', '_'")
    for key in ("material", "decided_by", "reason"):
        if not _text(request[key]):
            raise Refused("request", f"{key} must be non-empty text")
    _date(request["date"], "date")
    source = request["source"]
    if (
        not isinstance(source, dict)
        or set(source) != {"citation", "location"}
        or not all(_text(v) for v in source.values())
    ):
        raise Refused("request", "source needs exactly citation and location (text)")
    rights = request["rights"]
    if not isinstance(rights, dict):
        raise Refused("rights", "rights must be an object")
    required = {"reviewed_by", "date", "statement"}
    optional = {"same_basis", "differences"}
    if kind == "add_record":
        required.add("same_basis")
        if rights.get("same_basis") is False:
            required.add("differences")
    if not required <= set(rights) or set(rights) - required - optional:
        raise Refused(
            "rights",
            f"{kind} needs a rights statement with {sorted(required)}",
        )
    for key in ("reviewed_by", "statement", "differences"):
        if key in rights and not _text(rights[key]):
            raise Refused("rights", f"rights.{key} must be non-empty text")
    _date(rights["date"], "rights.date")
    if "same_basis" in rights and not isinstance(rights["same_basis"], bool):
        raise Refused("rights", "rights.same_basis must be true or false")
    if kind == "set_fields":
        changes = request["changes"]
        if not isinstance(changes, list) or not changes:
            raise Refused("request", "changes must be a non-empty list")
        for change in changes:
            if not isinstance(change, dict) or set(change) - set(CHANGE_KEYS):
                raise Refused(
                    "undeclared_change", f"a change takes only {list(CHANGE_KEYS)}"
                )
    elif not isinstance(request["record"], dict):
        raise Refused("request", "record must be an object")


def read_request(path: Path) -> dict:
    request = json.loads(Path(path).read_text(encoding="utf-8"))
    check_request(request)
    return request


# Files ------------------------------------------------------------------------


def _read(root: Path, material: str) -> tuple[dict, dict]:
    record_path, audit_path = (root / p for p in _paths(material))
    if not record_path.exists() or not audit_path.exists():
        raise Refused("unknown_material", f"{material} has no record file and sidecar")
    record_text_ = record_path.read_text(encoding="utf-8")
    audit_text = audit_path.read_text(encoding="utf-8")
    document, sidecar = json.loads(record_text_), json.loads(audit_text)
    if format_eosmat(document) != record_text_ or format_audit(sidecar) != audit_text:
        raise Refused("not_canonical", f"{material}: files are not in canonical layout")
    return document, sidecar


def _json_text(value: Any) -> bytes:
    return (json.dumps(value, indent=2, ensure_ascii=False) + "\n").encode("utf-8")


@contextmanager
def _stage(root: Path) -> Iterator[Path]:
    """A scratch copy of everything a revision reads or writes."""
    with tempfile.TemporaryDirectory(prefix="catalog-revision-") as name:
        stage = Path(name)
        for relative in ("peritheos/data/materials", "peritheos/data/audits"):
            shutil.copytree(root / relative, stage / relative)
        datasets = root / "peritheos/data/datasets"
        if datasets.exists():
            (stage / "peritheos/data/datasets").symlink_to(
                datasets.resolve(), target_is_directory=True
            )
        for relative in (RECEIPT_PATH, LEDGER, POLICY, REVIEW, INDEX, GROUPS):
            if (root / relative).exists():
                (stage / relative).parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(root / relative, stage / relative)
        yield stage


_WRITTEN = (LEDGER, RECEIPT_PATH, POLICY, REVIEW, INDEX, GROUPS)


def _changed_files(stage: Path, root: Path, material: str) -> dict[str, bytes]:
    files = {}
    for relative in (*_paths(material), *_WRITTEN):
        new = (stage / relative).read_bytes() if (stage / relative).exists() else None
        old = (root / relative).read_bytes() if (root / relative).exists() else None
        if new is not None and new != old:
            files[relative] = new
    return files


def _review(stage: Path) -> dict[str, dict[str, str]]:
    """Regenerate the rights review in the stage; return the output hashes."""
    hashes = {}
    try:
        outputs = build_review(stage)
    except (ValueError, KeyError, TypeError) as error:
        raise Refused("rights_review", f"the rights review cannot be built: {error!r}")
    for name, content in outputs.items():
        before = (stage / name).read_bytes() if (stage / name).exists() else b""
        (stage / name).parent.mkdir(parents=True, exist_ok=True)
        (stage / name).write_bytes(content)
        hashes[name] = {"before": _sha(before), "after": _sha(content)}
    return hashes


def _snapshot_parts(snapshot: dict, material: str) -> tuple[list, dict]:
    """The material's record projections and the snapshot without them."""
    rest = copy.deepcopy(snapshot)
    records: list = []
    for document in rest["documents"]:
        if document["identifier"] == material:
            records = document.pop("eos_records")
    return records, rest


# Building a revision -----------------------------------------------------


def _identifiers(root: Path) -> set[str]:
    names: set[str] = set()
    for path in (root / "peritheos/data/materials").glob("*.eosmat"):
        for record in json.loads(path.read_bytes()).get("eos_records", []):
            names.add(record.get("identifier"))
            names.update(record.get("aliases") or [])
    return names


def _add_record(
    material: str, document: dict, sidecar: dict, merged: dict, record: Any, root: Path
) -> list[dict]:
    """Append a complete inline record to the record file and the sidecar."""
    identifier = record.get("identifier")
    if not isinstance(identifier, str) or not _RECORD_ID.match(identifier):
        raise Refused("record", "identifier must be lower_snake_case")
    for key, reason in (
        ("check", "check is generated from the evidence"),
        ("superseded_record_text", "a new record has no superseded text"),
        ("default", "default flags would change the default selection"),
        ("default_for", "default flags would change the default selection"),
    ):
        if key in record:
            raise Refused(
                "undeclared_change" if key.startswith("default") else "record",
                reason,
            )
    taken = _identifiers(root)
    for name in [identifier, *(record.get("aliases") or [])]:
        if name in taken:
            raise Refused("identifier_taken", f"{name} is already used in the catalog")
    evidence = record.get("evidence")
    validation = (
        evidence.get("scientific_validation") if isinstance(evidence, dict) else None
    )
    source = (
        validation.get("primary_source_check") if isinstance(validation, dict) else None
    )
    if not isinstance(source, dict) or not _text(source.get("finding")):
        raise Refused(
            "audit",
            "a bundled record needs evidence.scientific_validation with a status "
            "and primary_source_check.finding",
        )
    if not _text(validation.get("status")):
        raise Refused("audit", "evidence.scientific_validation.status is required")
    reproduction = evidence.get("reproduction")
    assessment = (
        reproduction.get("assessment") if isinstance(reproduction, dict) else None
    )
    if not isinstance(assessment, dict) or not all(
        _text(assessment.get(key)) for key in ("status", "scope", "basis")
    ):
        raise Refused(
            "audit",
            "a bundled record needs evidence.reproduction.assessment with status, "
            "scope and basis",
        )
    if not isinstance(evidence.get("parameter_provenance"), dict):
        raise Refused(
            "audit",
            "a bundled record needs evidence.parameter_provenance (where each value "
            "comes from); the rights review binds it",
        )
    if set(assessment.get("qualifications") or []) & set(GENERIC_QUALIFICATIONS) or (
        validation.get("note") in GENERIC_VALIDATION_NOTES
    ):
        raise Refused(
            "generic_statement",
            "the generic statements removed by the audit split are stated once in "
            "docs/eos-labels.md, not per record",
        )
    inline = copy.deepcopy(record)
    audit = {key: inline.pop(key) for key in AUDIT_RECORD_KEYS if key in inline}
    audit = expand_entry(audit, material, identifier)
    complete = {**inline, **audit}
    check = check_block(
        complete, merged, audit=bundled_audit_reference(material, identifier)
    )
    cases = record_cases({**complete, "check": check})
    labels = [r.get("label") for r in document["eos_records"]]
    if record.get("label") is not None and record["label"] in labels:
        cases.append("label_unique")
    if cases:
        raise Refused(
            "record_text",
            f"outside the record text standard: {cases} (docs/record-text-standard.md)",
        )
    document["eos_records"].append({**inline, "check": check})
    sidecar["records"][identifier] = compact_entry(audit, material, identifier)
    # The complete record is filled in from the written files.
    return [{"op": "add_record", "record": {"identifier": identifier}}]


def _set_fields(document: dict, changes: list[dict]) -> list[dict]:
    """Apply declared field changes to the record file; return the operations."""
    records = {r["identifier"]: r for r in document["eos_records"]}
    forward: list[dict] = []
    for change in changes:
        identifier = change.get("record")
        if identifier not in records:
            raise Refused("unknown_record", f"{identifier!r} is not in the material")
        path = change.get("path")
        if (
            not isinstance(path, list)
            or not path
            or not all(isinstance(k, str) and k for k in path)
        ):
            raise Refused("path", "path is a non-empty list of field names")
        if path[0] in PROTECTED_FIELDS or tuple(path[:2]) in PROTECTED_PATHS:
            raise Refused(
                "protected_field",
                f"{'.'.join(path)}: identifier, text, check, default flags and audit "
                "blocks are not changed by a catalog revision",
            )
        for other in forward:
            shorter = min(len(other["path"]), len(path))
            if (
                other["record"] == identifier
                and other["path"][:shorter] == path[:shorter]
            ):
                raise Refused("overlapping_paths", f"{identifier}/{'.'.join(path)}")
        if ("expected" in change) == bool(change.get("expected_absent")):
            raise Refused("expected", "give exactly one of expected or expected_absent")
        if ("value" in change) == bool(change.get("remove")):
            raise Refused("value", "give exactly one of value or remove")
        try:
            current = read_path(records[identifier], path)
        except ValueError as error:
            raise Refused("path", str(error)) from None
        expected = ABSENT if change.get("expected_absent") else change["expected"]
        if not _same(current, expected):
            raise Refused(
                "expected_mismatch",
                f"{identifier}/{'.'.join(path)} is "
                f"{'absent' if current is ABSENT else json.dumps(current)}",
            )
        after_present = "value" in change
        if _same(current, change["value"] if after_present else ABSENT):
            raise Refused("no_change", f"{identifier}/{'.'.join(path)} already has it")
        operation = {
            "op": "set",
            "record": identifier,
            "path": list(path),
            "before_present": current is not ABSENT,
            "before": None if current is ABSENT else current,
            "after_present": after_present,
            "after": copy.deepcopy(change["value"]) if after_present else None,
        }
        apply_operations(document["eos_records"], [operation], "document")
        forward.append(operation)
    return forward


def _refresh_checks(document: dict, sidecar: dict, identifiers: set[str]) -> None:
    merged = merge_audit(document, sidecar)
    complete = {r["identifier"]: r for r in merged["eos_records"]}
    for record in document["eos_records"]:
        if record["identifier"] in identifiers:
            record["check"] = check_block(
                {
                    k: v
                    for k, v in complete[record["identifier"]].items()
                    if k != "check"
                },
                merged,
                audit=record["check"].get("audit"),
            )


def _set_records(operations: list[dict]) -> set[str]:
    return {o["record"] for o in operations if o["op"] == "set"}


def _with_checks(document: dict, source: dict, identifiers: set[str]) -> dict:
    """``document`` with the named records' generated ``check`` from ``source``."""
    result = copy.deepcopy(document)
    checks = {
        r["identifier"]: r.get("check")
        for r in source["eos_records"]
        if r["identifier"] in identifiers
    }
    for record in result["eos_records"]:
        if record["identifier"] in checks:
            record["check"] = copy.deepcopy(checks[record["identifier"]])
    return result


def _revised_binding(
    binding: dict,
    document: dict,
    forward: list[dict],
    files: tuple[bytes, bytes],
    reference: bytes,
) -> dict:
    """The audit-split binding of the revised files (see catalog_revisions)."""
    result = copy.deepcopy(binding)
    result["after_record_sha256"] = _sha(files[0])
    result["after_audit_sha256"] = _sha(files[1])
    order = [r["identifier"] for r in document["eos_records"]]
    for operation in forward:
        if operation["op"] == "add_record":
            # Its check block is generated, not part of the reference.
            index = order.index(operation["record"]["identifier"])
            result["added_check"] = sorted({*result.get("added_check", []), index})
        elif operation["before_present"] != operation["after_present"]:
            # A recorded key order of the node must name the new key set.
            parent = ["eos_records", order.index(operation["record"])]
            parent += operation["path"][:-1]
            for ordering in result.get("object_orders", []):
                if ordering["path"] == parent:
                    key = operation["path"][-1]
                    if operation["after_present"]:
                        ordering["keys"].append(key)
                    else:
                        ordering["keys"].remove(key)
    ordered: dict[str, Any] = {}
    for key, value in result.items():
        if key == "catalog_revision_sha256":
            continue
        ordered[key] = value
        if key == "before_sha256":
            ordered["catalog_revision_sha256"] = _sha(reference)
    return ordered


def plan_revision(request: dict, root: Path = ROOT) -> Plan:
    """Compute and verify one revision in a staging copy; write nothing."""
    check_request(request)
    material = request["material"]
    kind = request["kind"]
    ledger = read_ledger(root)
    if any(r["id"] == request["id"] for r in ledger["revisions"]):
        raise Refused("duplicate_id", f"{request['id']} is already in {LEDGER}")
    document, sidecar = _read(root, material)
    receipt = json.loads((root / RECEIPT_PATH).read_bytes())
    binding = receipt["materials"].get(material)
    if binding is None:
        raise Refused("not_split", f"{material} has no audit-split binding")
    merged_before = merge_audit(document, sidecar)
    record_path, audit_path = _paths(material)
    files_before = ((root / record_path).read_bytes(), (root / audit_path).read_bytes())
    reference_before = pre_split_bytes(material, root=root)
    if kind == "add_record":
        forward = _add_record(
            material, document, sidecar, merged_before, request["record"], root
        )
    else:
        forward = _set_fields(document, request["changes"])
        after = {r["identifier"]: r for r in document["eos_records"]}
        for record in merged_before["eos_records"]:
            for name in SUPERSEDED_TEXT_FIELDS:
                if record_text(record, name) != record_text(
                    after[record["identifier"]], name
                ):
                    raise Refused(
                        "protected_field",
                        f"{record['identifier']}: {name} would change; use "
                        "scripts/apply_record_text.py",
                    )
    set_records = _set_records(forward)
    _refresh_checks(document, sidecar, set_records)
    files_after = (
        format_eosmat(document).encode("utf-8"),
        format_audit(sidecar).encode("utf-8"),
    )
    record_file = json.loads(files_after[0])
    try:
        validate_eosmat_document(record_file)
        merged_after = merge_audit(record_file, json.loads(files_after[1]))
        validate_eosmat_document(merged_after)
    except Exception as error:  # noqa: BLE001 - reported as a refusal
        raise Refused("invalid", str(error)) from None
    if kind == "add_record":
        identifier = forward[0]["record"]["identifier"]
        forward[0]["record"] = next(
            r for r in merged_after["eos_records"] if r["identifier"] == identifier
        )
        _executable(merged_after, identifier)
    reverse = reverse_operations(forward)
    # Exactly the declared change: the forward operations turn the complete
    # document into the result, the reverse operations turn it back.
    trial = copy.deepcopy(merged_before)
    apply_operations(trial["eos_records"], forward, "document")
    trial = _with_checks(trial, merged_after, set_records)
    if format_eosmat(trial) != format_eosmat(merged_after) or fingerprint(
        trial
    ) != fingerprint(merged_after):
        raise Refused("undeclared_change", "the files differ beyond the operations")
    back = copy.deepcopy(merged_after)
    apply_operations(back["eos_records"], reverse, "document")
    back = _with_checks(back, merged_before, set_records)
    if format_eosmat(back) != format_eosmat(merged_before):
        raise Refused(
            "not_reversible", "the reverse operations do not restore the file"
        )
    # The pre-split reference: forward operations on the current reference
    # and recovery from the new files through the revised binding must agree.
    reference_document = json.loads(reference_before)
    apply_operations(reference_document["eos_records"], forward, "reference")
    layout = LAYOUTS[binding.get("layout")]
    reference_after = format_eosmat(reference_document, layout=layout).encode("utf-8")
    new_binding = _revised_binding(
        binding, document, forward, files_after, reference_after
    )
    try:
        recover_pre_split(material, *files_after, new_binding)
    except ValueError as error:
        raise Refused(
            "not_bindable",
            f"the audit-split binding cannot restore the revised reference ({error})",
        ) from None
    new_receipt = copy.deepcopy(receipt)
    new_receipt["materials"][material] = new_binding
    new_receipt.setdefault(*SPLIT_RECEIPT_STATEMENT)
    previous = [r for r in ledger["revisions"] if r["material"] == material]
    with _stage(root) as stage:
        snapshot_before = scientific_snapshot(stage)
        (stage / record_path).write_bytes(files_after[0])
        (stage / audit_path).write_bytes(files_after[1])
        (stage / RECEIPT_PATH).write_text(receipt_text(new_receipt), encoding="utf-8")
        snapshot_after = scientific_snapshot(stage)
        records_before, rest_before = _snapshot_parts(snapshot_before, material)
        records_after, rest_after = _snapshot_parts(snapshot_after, material)
        if fingerprint(rest_before) != fingerprint(rest_after):
            raise Refused(
                "undeclared_change", "the scientific snapshot changed elsewhere"
            )
        trial_records = copy.deepcopy(records_before)
        apply_operations(trial_records, forward, "snapshot")
        if fingerprint(trial_records) != fingerprint(records_after):
            raise Refused(
                "undeclared_change",
                "the scientific snapshot differs from the operations",
            )
        policy = json.loads((stage / POLICY).read_bytes())
        pin = policy["material_sha256"].get(record_path)
        if pin != _sha(reference_before):
            raise Refused(
                "rights_review_stale",
                f"{POLICY} does not pin the current reference of {material}; run "
                "scripts/review_eos_structure_rights.py --check",
            )
        if policy["catalog_scientific_sha256"] != fingerprint(snapshot_before):
            raise Refused(
                "rights_review_stale",
                f"{POLICY} catalog_scientific_sha256 is not the current snapshot",
            )
        policy["material_sha256"][record_path] = _sha(reference_after)
        policy["catalog_scientific_sha256"] = fingerprint(snapshot_after)
        pinned = policy.setdefault(
            "catalog_revisions", {"ledger": LEDGER, "materials": {}}
        )
        pinned["materials"][material] = [
            *pinned["materials"].get(material, []),
            request["id"],
        ]
        pinned["materials"] = dict(sorted(pinned["materials"].items()))
        (stage / POLICY).write_bytes(_json_text(policy))
        derived = _review(stage)
        revision = {
            **{key: copy.deepcopy(request[key]) for key in COMMON_KEYS},
            "forward": forward,
            "reverse": reverse,
            "document_sha256": {
                "before": _sha(format_eosmat(merged_before).encode("utf-8")),
                "after": _sha(format_eosmat(merged_after).encode("utf-8")),
            },
            "reference_sha256": {
                "before": _sha(reference_before),
                "after": _sha(reference_after),
            },
            "scientific": {
                "records_before": fingerprint(records_before),
                "records_after": fingerprint(records_after),
                "added_records": {
                    r["identifier"]: fingerprint(r)
                    for r in records_after
                    if any(
                        o["op"] == "add_record"
                        and o["record"]["identifier"] == r["identifier"]
                        for o in forward
                    )
                },
                "catalog_before": fingerprint(snapshot_before),
                "catalog_after": fingerprint(snapshot_after),
            },
            "files": {
                record_path: {
                    "before": _sha(files_before[0]),
                    "after": _sha(files_after[0]),
                },
                audit_path: {
                    "before": _sha(files_before[1]),
                    "after": _sha(files_after[1]),
                },
            },
            "split_binding": {
                "before": copy.deepcopy(binding),
                "after_sha256": fingerprint(new_binding),
            },
            "policy": {
                "material_sha256": {
                    "before": _sha(reference_before),
                    "after": _sha(reference_after),
                },
                "catalog_scientific_sha256": {
                    "before": fingerprint(snapshot_before),
                    "after": fingerprint(snapshot_after),
                },
            },
            "derived_files": derived,
            "previous_sha256": fingerprint(previous[-1]) if previous else None,
        }
        new_ledger = copy.deepcopy(ledger)
        new_ledger["revisions"].append(revision)
        (stage / LEDGER).parent.mkdir(parents=True, exist_ok=True)
        (stage / LEDGER).write_text(ledger_text(new_ledger), encoding="utf-8")
        _verify_stage(stage, root, material, snapshot_before)
        plan = Plan(revision=revision, undo=False)
        plan.files = _changed_files(stage, root, material)
    plan.shown = _show(revision, merged_after, set_records)
    return plan


def _executable(document: dict, identifier: str) -> None:
    from peritheos.materials import Material

    try:
        material = Material.from_eosmat(document, require_primary_validation=False)
        record = next(r for r in material.eos_records if r.identifier == identifier)
        record.eos  # noqa: B018 - constructing the model is the check
    except Exception as error:  # noqa: BLE001 - reported as a refusal
        raise Refused("not_executable", f"{identifier}: {error}") from None


def _verify_stage(stage: Path, root: Path, material: str, before: dict) -> None:
    """Every recovery check on the staged files; raises on the first failure."""
    from scripts.migrate_eosmat_v5 import prior_material_bytes

    failures = verify_ledger(stage, [material])
    if failures:
        raise Refused("lineage", "; ".join(failures))
    if prior_material_bytes(material, root=stage) != prior_material_bytes(
        material, root=root
    ):
        raise Refused("lineage", "the v5-migration predecessor would change")
    # Reversing every revision gives the same snapshot as before this one.
    now = scientific_snapshot(stage)
    base = prior_revision_snapshot(before, root) or before
    if fingerprint(prior_revision_snapshot(now, stage) or now) != fingerprint(base):
        raise Refused("lineage", "the reversed scientific snapshot would change")
    outputs = build_review(stage)
    if any((stage / name).read_bytes() != content for name, content in outputs.items()):
        raise Refused("lineage", "the rights review outputs are stale")
    document, sidecar = _read(stage, material)
    merged = merge_audit(document, sidecar)
    for record, complete in zip(document["eos_records"], merged["eos_records"]):
        derived = check_block(
            {k: v for k, v in complete.items() if k != "check"},
            merged,
            audit=record["check"].get("audit"),
        )
        if record["check"] != derived:
            raise Refused("lineage", f"{record['identifier']}: check differs")


def _show(revision: dict, merged: dict, set_records: set[str]) -> list[str]:
    lines = []
    records = {r["identifier"]: r for r in merged["eos_records"]}
    for operation in revision["forward"]:
        if operation["op"] == "add_record":
            lines.append(f"Added record {operation['record']['identifier']}:")
            # Rendered in record order, as it will read in the complete document.
            shown = {"format_version": 5, "eos_records": [operation["record"]]}
            lines.append(format_eosmat(shown).rstrip())
        else:
            before = (
                json.dumps(operation["before"])
                if operation["before_present"]
                else "absent"
            )
            after = (
                json.dumps(operation["after"])
                if operation["after_present"]
                else "absent"
            )
            lines.append(
                f"Set {operation['record']}/{'.'.join(operation['path'])}: "
                f"{before} -> {after}"
            )
    for identifier in sorted(
        set_records
        | {
            o["record"]["identifier"]
            for o in revision["forward"]
            if o["op"] == "add_record"
        }
    ):
        lines.append(f"check of {identifier}:")
        lines.append(
            json.dumps(records[identifier]["check"], indent=2, ensure_ascii=False)
        )
    return lines


# Undo ---------------------------------------------------------------------------


def plan_undo(revision_id: str, root: Path = ROOT) -> Plan:
    """Compute the exact undo of a material's most recent revision; write nothing."""
    ledger = read_ledger(root)
    revision = next((r for r in ledger["revisions"] if r["id"] == revision_id), None)
    if revision is None:
        raise Refused("unknown_revision", f"{revision_id} is not in {LEDGER}")
    material = revision["material"]
    later = [
        r["id"]
        for r in ledger["revisions"][ledger["revisions"].index(revision) + 1 :]
        if r["material"] == material
    ]
    if later:
        raise Refused("out_of_order", f"undo the later revisions {later} first")
    record_path, audit_path = _paths(material)
    for relative in (record_path, audit_path):
        if _sha((root / relative).read_bytes()) != revision["files"][relative]["after"]:
            raise Refused(
                "changed_since",
                f"{relative} changed after {revision_id} (a later text revision?); "
                "undo that first",
            )
    receipt = json.loads((root / RECEIPT_PATH).read_bytes())
    if (
        fingerprint(receipt["materials"][material])
        != revision["split_binding"]["after_sha256"]
    ):
        raise Refused("changed_since", f"the audit-split binding of {material} changed")
    document, sidecar = _read(root, material)
    for operation in revision["reverse"]:
        if operation["op"] == "remove_record":
            sidecar["records"].pop(operation["record"], None)
    apply_operations(document["eos_records"], revision["reverse"], "document")
    _refresh_checks(document, sidecar, _set_records(revision["reverse"]))
    files = (
        format_eosmat(document).encode("utf-8"),
        format_audit(sidecar).encode("utf-8"),
    )
    for relative, content in zip((record_path, audit_path), files):
        if _sha(content) != revision["files"][relative]["before"]:
            raise Refused(
                "not_exact", f"{relative} would not return to its earlier bytes"
            )
    receipt["materials"][material] = copy.deepcopy(revision["split_binding"]["before"])
    if not any("catalog_revision_sha256" in b for b in receipt["materials"].values()):
        receipt.pop(SPLIT_RECEIPT_STATEMENT[0], None)
    new_ledger = copy.deepcopy(ledger)
    new_ledger["revisions"] = [r for r in ledger["revisions"] if r["id"] != revision_id]
    with _stage(root) as stage:
        (stage / record_path).write_bytes(files[0])
        (stage / audit_path).write_bytes(files[1])
        (stage / RECEIPT_PATH).write_text(receipt_text(receipt), encoding="utf-8")
        (stage / LEDGER).write_text(ledger_text(new_ledger), encoding="utf-8")
        policy = json.loads((stage / POLICY).read_bytes())
        policy["material_sha256"][record_path] = revision["policy"]["material_sha256"][
            "before"
        ]
        pinned = policy["catalog_revisions"]["materials"]
        pinned[material] = [i for i in pinned.get(material, []) if i != revision_id]
        if not pinned[material]:
            del pinned[material]
        if not pinned:
            del policy["catalog_revisions"]
        policy["catalog_scientific_sha256"] = fingerprint(scientific_snapshot(stage))
        (stage / POLICY).write_bytes(_json_text(policy))
        _review(stage)
        if (
            _sha(pre_split_bytes(material, root=stage))
            != revision["reference_sha256"]["before"]
        ):
            raise Refused("not_exact", "the pre-split reference would not return")
        records, _ = _snapshot_parts(scientific_snapshot(stage), material)
        if fingerprint(records) != revision["scientific"]["records_before"]:
            raise Refused("not_exact", "the scientific snapshot would not return")
        failures = verify_ledger(stage, [material])
        if failures:
            raise Refused("lineage", "; ".join(failures))
        plan = Plan(revision=revision, undo=True)
        plan.files = _changed_files(stage, root, material)
    plan.shown = [f"Undo {revision_id} ({revision['kind']}) of {material}"]
    return plan


# Writing --------------------------------------------------------------------------


def apply_plan(plan: Plan, root: Path = ROOT) -> None:
    """Write a verified plan; on any failed check of the written tree, put back."""
    originals = {
        relative: (root / relative).read_bytes() if (root / relative).exists() else None
        for relative in plan.files
    }
    material = plan.revision["material"]
    try:
        for relative, content in plan.files.items():
            (root / relative).parent.mkdir(parents=True, exist_ok=True)
            (root / relative).write_bytes(content)
        failures = verify_ledger(root, [material])
        if failures:
            raise ValueError("; ".join(failures))
        for name, content in build_review(root).items():
            if (root / name).read_bytes() != content:
                raise ValueError(f"{name} is stale after writing")
        pre_split_bytes(material, root=root)
    except Exception:
        for relative, content in originals.items():
            if content is None:
                (root / relative).unlink(missing_ok=True)
            else:
                (root / relative).write_bytes(content)
        raise


def verify(root: Path = ROOT) -> list[str]:
    """The whole ledger against the catalog (see ``verify_ledger``)."""
    return verify_ledger(root)


def _report(plan: Plan, check: bool) -> None:
    revision = plan.revision
    verb = "Would" if check else "Did"
    action = "undo" if plan.undo else "apply"
    print(
        f"{verb} {action} {revision['id']} ({revision['kind']}) for "
        f"{revision['material']}."
    )
    for line in plan.shown:
        print(line)
    print("Files:")
    for relative in plan.files:
        print(f"  {relative}")
    print("Then regenerate, in this order (each takes --check):")
    for number, (command, what) in enumerate(REGENERATE, start=1):
        print(f"  {number}. {command}  # {what}")
    if any(o["op"] == "add_record" for o in revision["forward"]):
        print("An added record must also be reflected by hand in:")
        for item in UPDATE_BY_HAND:
            print(f"  - {item}")
    print(
        "Then run python scripts/apply_catalog_revision.py --verify and the full "
        "checks in CONTRIBUTING.md."
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("request", nargs="?", type=Path, help="request JSON file")
    parser.add_argument("--check", action="store_true", help="dry run")
    parser.add_argument("--undo", metavar="ID", help="undo this revision")
    parser.add_argument("--verify", action="store_true", help="verify the ledger")
    arguments = parser.parse_args(argv)
    if arguments.verify:
        failures = verify()
        for failure in failures:
            print(f"FAILED {failure}")
        count = len(read_ledger()["revisions"])
        print(
            f"Verified {count} catalog revisions: "
            f"{'ok' if not failures else f'{len(failures)} failures'}."
        )
        return 1 if failures else 0
    if (arguments.request is None) == (arguments.undo is None):
        parser.error("give a request file or --undo ID")
    try:
        plan = (
            plan_undo(arguments.undo)
            if arguments.undo
            else plan_revision(read_request(arguments.request))
        )
        if not arguments.check:
            apply_plan(plan)
    except Refused as error:
        print(f"REFUSED {error}; nothing was written.")
        return 1
    _report(plan, arguments.check)
    return 0


if __name__ == "__main__":
    sys.exit(main())


__all__ = [
    "PROTECTED_FIELDS",
    "REGENERATE",
    "UPDATE_BY_HAND",
    "Plan",
    "Refused",
    "apply_plan",
    "check_request",
    "main",
    "plan_revision",
    "plan_undo",
    "read_request",
    "verify",
]
