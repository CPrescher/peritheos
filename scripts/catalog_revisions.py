"""Repository maintenance: the catalog revision ledger and its exact reversal.

``docs/data/catalog-revisions.json`` records every change of the scientific
content of a bundled material made after the audit split: an added EOS record
(``add_record``) or set or corrected structured record fields
(``set_fields``). ``scripts/apply_catalog_revision.py`` writes it; this module
reads it for the historical checks.

Each revision stores its forward and reverse operations on the complete inline
document, and the sha256 of three views before and after: the complete
canonical inline document (record file merged with its sidecar), the material's
*pre-split reference* and the scientific snapshot of its records. The pre-split
reference is what :func:`scripts.split_eosmat_audits.pre_split_bytes` recovers
from the current files: the bound text revisions undone, ``check`` blocks
removed, default flags and generic statements restored, rendered in the
binding's layout. Text revisions never change it, so the reference hashes of a
material's revisions form an unbroken chain whatever text revisions come
between them.

Invariant (binding): a catalog revision never moves a material's audit-split
binding to ``superseded``. Its ``before_sha256`` stays the original pre-split
file and is the first revision's reference ``before``; the binding's
``catalog_revision_sha256`` names the last revision's reference ``after``, and
restoring the current files with the binding gives exactly that reference. An
added record is appended, so the index-based binding entries of the other
records keep their meaning; its index joins ``added_check``. Text revisions
(``scripts/apply_record_text.py``, ``scripts/apply_check_text.py``) keep
updating the same binding, and ``scripts/split_eosmat_audits.py --add``
refuses a material that has catalog revisions.

:func:`prior_revision_bytes` walks a material's revisions backwards from its
current reference, verifying every recorded hash, and returns the original
pre-split file the older receipts expect. :func:`prior_revision_snapshot` does
the same for the scientific snapshot. With no revisions both return their input
unchanged (``None`` for the snapshot), so every existing check behaves as
before. Not an installed-package tool.
"""

from __future__ import annotations

import copy
import hashlib
import json
from collections.abc import Iterable, Mapping
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
LEDGER = "docs/data/catalog-revisions.json"
SCHEMA = "peritheos.catalog-revisions.v1"
#: Revision kinds this ledger can express. ``remove_record`` and material-level
#: fields (``default_records``) are not supported yet.
KINDS = ("add_record", "set_fields")
SCOPE = (
    "Changes of the scientific content of bundled materials made after the "
    "audit split, one revision per entry, in the order applied. Each revision "
    "names who decided it, why, the source it rests on and a rights "
    "statement, and stores explicit forward and reverse operations on the "
    "complete inline document with the sha256 of the complete canonical "
    "document, of the pre-split reference (scripts/split_eosmat_audits.py "
    "pre_split_bytes) and of the scientific snapshot of the material's records "
    "before and after. previous_sha256 chains the revisions of one material. "
    "Reversing every revision of a material recovers its original pre-split "
    "file (the audit-split binding's before_sha256) byte for byte. Written by "
    "scripts/apply_catalog_revision.py."
)
#: Statement added to the audit-split receipt while a binding holds
#: ``catalog_revision_sha256``.
SPLIT_RECEIPT_STATEMENT = (
    "catalog_revision",
    "A binding's catalog_revision_sha256 is the sha256 of the pre-split "
    "reference after the material's latest revision in "
    "docs/data/catalog-revisions.json; restoring the current files with the "
    "binding gives that reference, and reversing the material's revisions "
    "recovers before_sha256, which is unchanged.",
)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def fingerprint(value: Any) -> str:
    """sha256 of canonical JSON, as the scientific snapshot is hashed."""
    return _sha(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode("utf-8")
    )


def empty_ledger() -> dict:
    return {"schema": SCHEMA, "scope": SCOPE, "revisions": []}


def ledger_text(ledger: Mapping[str, Any]) -> str:
    """The ledger file's text."""
    return json.dumps(ledger, indent=2, ensure_ascii=False) + "\n"


def read_ledger(root: Path | None = None) -> dict:
    """The ledger at ``root`` (an empty ledger when the file is absent)."""
    path = (ROOT if root is None else root) / LEDGER
    if not path.exists():
        return empty_ledger()
    ledger = json.loads(path.read_bytes())
    if ledger.get("schema") != SCHEMA:
        raise ValueError(f"{LEDGER}: unknown schema {ledger.get('schema')!r}")
    return ledger


def material_revisions(material: str, root: Path | None = None) -> list[dict]:
    """A material's revisions, oldest first."""
    return [r for r in read_ledger(root)["revisions"] if r["material"] == material]


def added_record(
    material: str, identifier: str, root: Path | None = None
) -> dict | None:
    """The complete inline record as a revision added it, or ``None``."""
    for revision in material_revisions(material, root):
        for operation in revision["forward"]:
            if (
                operation["op"] == "add_record"
                and operation["record"]["identifier"] == identifier
            ):
                return copy.deepcopy(operation["record"])
    return None


# Operations ---------------------------------------------------------------

ABSENT = object()


def _same(left: Any, right: Any) -> bool:
    """Equal including JSON type (``1`` and ``1.0``, ``true`` and ``1`` differ)."""
    if left is ABSENT or right is ABSENT:
        return left is right
    return json.dumps(left, sort_keys=True) == json.dumps(right, sort_keys=True)


def _index(records: list, identifier: str) -> int:
    matches = [i for i, r in enumerate(records) if r.get("identifier") == identifier]
    if len(matches) != 1:
        raise ValueError(f"Record {identifier!r} is not in the document exactly once")
    return matches[0]


def _parent(record: dict, path: list[str]) -> dict:
    node = record
    for key in path[:-1]:
        node = node.get(key) if isinstance(node, dict) else None
        if not isinstance(node, dict):
            raise ValueError(
                f"Record {record.get('identifier')!r} has no object at {path[:-1]}"
            )
    return node


def read_path(record: dict, path: list[str]) -> Any:
    """The value at ``path`` in ``record``; ``ABSENT`` when the key is missing."""
    node = _parent(record, path)
    return copy.deepcopy(node[path[-1]]) if path[-1] in node else ABSENT


def _write(record: dict, path: list[str], present: bool, value: Any) -> None:
    node = _parent(record, path)
    if present:
        node[path[-1]] = copy.deepcopy(value)
    else:
        node.pop(path[-1], None)


def _view_record(record: Mapping[str, Any], view: str) -> dict:
    """An added record as it appears in one view."""
    if view == "document":
        return copy.deepcopy(dict(record))
    if view == "reference":
        # The pre-split reference holds no check blocks (``added_check``).
        return {k: copy.deepcopy(v) for k, v in record.items() if k != "check"}
    from scripts.scientific_source_snapshot import RECORD_FIELDS

    return {key: copy.deepcopy(record.get(key)) for key in RECORD_FIELDS}


def apply_operations(records: list, operations: Iterable[Mapping], view: str) -> None:
    """Apply operations to a material's ``eos_records`` list in place.

    ``view`` is ``document`` (complete inline records), ``reference`` (the
    pre-split reference: an added record without ``check``) or ``snapshot``
    (the scientific projection of ``scripts/scientific_source_snapshot.py``,
    where a missing top-level field reads ``null`` and a field outside the
    projection is skipped). A ``set`` operation first verifies the recorded
    ``before`` value; a removed record must be the last one.
    """
    if view not in ("document", "reference", "snapshot"):
        raise ValueError(f"Unknown view {view!r}")
    from scripts.scientific_source_snapshot import RECORD_FIELDS

    for operation in operations:
        kind = operation["op"]
        if kind == "add_record":
            record = _view_record(operation["record"], view)
            if any(r.get("identifier") == record["identifier"] for r in records):
                raise ValueError(f"Record {record['identifier']!r} already exists")
            records.append(record)
        elif kind == "remove_record":
            index = _index(records, operation["record"])
            if index != len(records) - 1:
                raise ValueError(
                    f"Record {operation['record']!r} is not the last record"
                )
            del records[index]
        elif kind == "set":
            path = list(operation["path"])
            if view == "snapshot" and path[0] not in RECORD_FIELDS:
                continue
            record = records[_index(records, operation["record"])]
            projected = view == "snapshot" and len(path) == 1
            before = operation["before"] if operation["before_present"] else ABSENT
            after_present = operation["after_present"]
            if projected:
                # The snapshot reads a missing field as null.
                before = None if before is ABSENT else before
                after_present = True
            current = read_path(record, path)
            if not _same(current, before):
                raise ValueError(
                    f"{operation['record']}/{'.'.join(path)} is not the recorded value"
                )
            _write(record, path, after_present, operation["after"])
        else:
            raise ValueError(f"Unknown operation {kind!r}")


def reverse_operations(forward: list[dict]) -> list[dict]:
    """The exact inverse of ``forward``, applied in reverse order."""
    result: list[dict] = []
    for operation in reversed(forward):
        if operation["op"] == "add_record":
            result.append(
                {"op": "remove_record", "record": operation["record"]["identifier"]}
            )
        elif operation["op"] == "set":
            result.append(
                {
                    "op": "set",
                    "record": operation["record"],
                    "path": list(operation["path"]),
                    "before_present": operation["after_present"],
                    "before": copy.deepcopy(operation["after"]),
                    "after_present": operation["before_present"],
                    "after": copy.deepcopy(operation["before"]),
                }
            )
        else:
            raise ValueError(f"No inverse for operation {operation['op']!r}")
    return result


# Reversal hooks -------------------------------------------------------------


def _binding(material: str, root: Path) -> dict:
    from scripts.split_eosmat_audits import RECEIPT
    from scripts.split_eosmat_audits import ROOT as SPLIT_ROOT

    path = root / RECEIPT.relative_to(SPLIT_ROOT)
    binding = json.loads(path.read_bytes())["materials"].get(material)
    if binding is None:
        raise ValueError(f"Catalog revisions of {material} have no audit-split binding")
    return binding


def prior_revision_bytes(data: bytes, material: str, root: Path | None = None) -> bytes:
    """Reverse a material's catalog revisions on its pre-split reference.

    ``data`` is the current reference (``pre_split_bytes``). Each revision,
    newest first, must find exactly its recorded ``after`` hash; its reverse
    operations are applied and the result must hash to its ``before``. The
    returned bytes are the original pre-split file (the binding's
    ``before_sha256``) that the v5-migration, garnet and rights receipts
    expect. A material without revisions returns ``data`` unchanged.
    """
    root = ROOT if root is None else root
    revisions = material_revisions(material, root)
    if not revisions:
        return data
    from peritheos._eosmat_format import format_eosmat
    from scripts.split_eosmat_audits import LAYOUTS

    binding = _binding(material, root)
    layout = LAYOUTS[binding.get("layout")]
    head = revisions[-1]["reference_sha256"]["after"]
    if binding.get("catalog_revision_sha256") != head:
        raise ValueError(
            f"Audit-split binding of {material} does not name its latest catalog "
            "revision"
        )
    if _sha(data) != head:
        raise ValueError(f"Catalog revision reference changed: {material}")
    for revision in reversed(revisions):
        if _sha(data) != revision["reference_sha256"]["after"]:
            raise ValueError(
                f"Catalog revision {revision['id']} does not follow its successor"
            )
        document = json.loads(data)
        apply_operations(document["eos_records"], revision["reverse"], "reference")
        data = format_eosmat(document, layout=layout).encode("utf-8")
        if _sha(data) != revision["reference_sha256"]["before"]:
            raise ValueError(
                f"Catalog revision {revision['id']} predecessor changed: {material}"
            )
    if _sha(data) != binding["before_sha256"]:
        raise ValueError(
            f"Catalog revisions do not reach the pre-split file: {material}"
        )
    return data


def prior_revision_snapshot(snapshot: dict, root: Path | None = None) -> dict | None:
    """Reverse every catalog revision in a scientific snapshot, verifying each step.

    A revised material's records must match its latest revision's
    ``records_after`` fingerprint (or already its first revision's
    ``records_before``); anything else raises. Each step checks the
    projection before and after, and each removed record's own fingerprint.
    Returns ``None`` when nothing was reversed.
    """
    ledger = read_ledger(root)
    if not ledger["revisions"]:
        return None
    by_material: dict[str, list[dict]] = {}
    for revision in ledger["revisions"]:
        by_material.setdefault(revision["material"], []).append(revision)
    result = copy.deepcopy(snapshot)
    changed = False
    for document in result["documents"]:
        revisions = by_material.get(document["identifier"])
        if not revisions:
            continue
        records = document["eos_records"]
        current = fingerprint(records)
        if current == revisions[-1]["scientific"]["records_after"]:
            pass
        elif current == revisions[0]["scientific"]["records_before"]:
            continue
        else:
            raise ValueError(
                f"Catalog revision snapshot changed: {document['identifier']}"
            )
        for revision in reversed(revisions):
            scientific = revision["scientific"]
            if fingerprint(records) != scientific["records_after"]:
                raise ValueError(
                    f"Catalog revision {revision['id']} snapshot does not follow "
                    "its successor"
                )
            for operation in revision["reverse"]:
                if operation["op"] == "remove_record":
                    record = records[_index(records, operation["record"])]
                    expected = scientific["added_records"].get(operation["record"])
                    if fingerprint(record) != expected:
                        raise ValueError(
                            f"Catalog revision {revision['id']} added record changed"
                        )
            apply_operations(records, revision["reverse"], "snapshot")
            if fingerprint(records) != scientific["records_before"]:
                raise ValueError(
                    f"Catalog revision {revision['id']} snapshot predecessor changed"
                )
        changed = changed or current != fingerprint(records)
    return result if changed else None


# Verification ---------------------------------------------------------------

REVISION_KEYS = (
    "id",
    "material",
    "kind",
    "date",
    "decided_by",
    "reason",
    "source",
    "rights",
    "forward",
    "reverse",
    "document_sha256",
    "reference_sha256",
    "scientific",
    "files",
    "split_binding",
    "policy",
    "derived_files",
    "previous_sha256",
)


def verify_ledger(
    root: Path | None = None, materials: Iterable[str] | None = None
) -> list[str]:
    """Check the ledger against the catalog; return the failures.

    For every revised material (or the named ones): the revision chain and
    the stored reverse operations, the audit-split binding (never superseded,
    ``catalog_revision_sha256`` naming the latest revision, ``before_sha256``
    the first revision's reference ``before``), the exact walk back to the
    original pre-split file, the scientific snapshot fingerprints and the
    rights-policy pins. Also checks that no binding names a revision the
    ledger lacks.
    """
    root = ROOT if root is None else root
    from scripts.review_eos_structure_rights import POLICY
    from scripts.split_eosmat_audits import RECEIPT, pre_split_bytes
    from scripts.split_eosmat_audits import ROOT as SPLIT_ROOT

    failures: list[str] = []
    ledger = read_ledger(root)
    receipt = json.loads((root / RECEIPT.relative_to(SPLIT_ROOT)).read_bytes())
    policy = json.loads((root / POLICY).read_bytes())
    ids = [r.get("id") for r in ledger["revisions"]]
    if len(ids) != len(set(ids)):
        failures.append("revision identifiers are not unique")
    by_material: dict[str, list[dict]] = {}
    for revision in ledger["revisions"]:
        if tuple(revision) != REVISION_KEYS:
            failures.append(f"{revision.get('id')}: keys differ from {REVISION_KEYS}")
            continue
        if revision["kind"] not in KINDS:
            failures.append(f"{revision['id']}: unknown kind {revision['kind']!r}")
        if revision["reverse"] != reverse_operations(revision["forward"]):
            failures.append(f"{revision['id']}: reverse is not the inverse of forward")
        by_material.setdefault(revision["material"], []).append(revision)
    bound = {
        m for m, b in receipt["materials"].items() if "catalog_revision_sha256" in b
    }
    for material in sorted(bound - set(by_material)):
        failures.append(f"{material}: binding names a revision the ledger lacks")
    statement, _ = SPLIT_RECEIPT_STATEMENT
    if bool(bound) != (statement in receipt):
        failures.append(f"audit-split receipt statement {statement!r} out of step")
    pinned = policy.get("catalog_revisions", {}).get("materials", {})
    for material in sorted(set(pinned) - set(by_material)):
        failures.append(f"{material}: rights policy names revisions the ledger lacks")
    selected = set(by_material) if materials is None else set(materials)
    for material in sorted(selected & set(by_material)):
        revisions = by_material[material]
        binding = receipt["materials"].get(material)
        if binding is None:
            failures.append(f"{material}: revised material has no audit-split binding")
            continue
        previous: dict | None = None
        for revision in revisions:
            expected_previous = None if previous is None else fingerprint(previous)
            if revision["previous_sha256"] != expected_previous:
                failures.append(f"{revision['id']}: previous_sha256 breaks the chain")
            start = (
                binding["before_sha256"]
                if previous is None
                else previous["reference_sha256"]["after"]
            )
            if revision["reference_sha256"]["before"] != start:
                failures.append(f"{revision['id']}: reference chain is broken")
            if previous is not None and (
                revision["scientific"]["records_before"]
                != previous["scientific"]["records_after"]
            ):
                failures.append(f"{revision['id']}: snapshot chain is broken")
            previous = revision
        try:
            prior_revision_bytes(pre_split_bytes(material, root=root), material, root)
        except (ValueError, KeyError, IndexError) as error:
            failures.append(f"{material}: {error}")
        relative = f"peritheos/data/materials/{material}.eosmat"
        last = revisions[-1]
        try:
            reference = _sha(pre_split_bytes(material, root=root))
        except ValueError as error:
            failures.append(f"{material}: {error}")
            reference = None
        if policy["material_sha256"].get(relative) != reference or (
            reference != last["policy"]["material_sha256"]["after"]
        ):
            failures.append(f"{material}: rights policy pin is not the latest revision")
        if pinned.get(material) != [r["id"] for r in revisions]:
            failures.append(f"{material}: rights policy does not list its revisions")
    if any(m in selected for m in by_material):
        from scripts.scientific_source_snapshot import scientific_snapshot

        snapshot = scientific_snapshot(root)
        documents = {d["identifier"]: d for d in snapshot["documents"]}
        for material in sorted(selected & set(by_material)):
            document = documents.get(material)
            last = by_material[material][-1]
            if document is None or (
                fingerprint(document["eos_records"])
                != last["scientific"]["records_after"]
            ):
                failures.append(
                    f"{material}: scientific snapshot is not the latest revision"
                )
        try:
            prior_revision_snapshot(snapshot, root)
        except (ValueError, KeyError, IndexError) as error:
            failures.append(str(error))
    return failures


__all__ = [
    "KINDS",
    "LEDGER",
    "added_record",
    "apply_operations",
    "fingerprint",
    "material_revisions",
    "prior_revision_bytes",
    "prior_revision_snapshot",
    "read_ledger",
    "reverse_operations",
    "verify_ledger",
]
