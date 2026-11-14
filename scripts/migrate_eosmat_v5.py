"""Repository maintenance: migrate the catalog to v5 with reversible receipts.

Standalone users can use ``python -m peritheos convert INPUT OUTPUT``; this
catalog receipt is not required to author, load or exchange an EOSMAT file.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from peritheos import convert_eosmat, validate_eosmat_document
from peritheos._eosmat_format import format_eosmat

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "docs/data/eosmat-v5-migration.json"


def fingerprint(value: dict) -> str:
    return hashlib.sha256(
        json.dumps(
            value, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()


def _changes(before, after, path=()):
    if before == after:
        return []
    if isinstance(before, dict) and isinstance(after, dict):
        result = []
        for key in dict.fromkeys([*before, *after]):
            if key in before and key in after:
                result.extend(_changes(before[key], after[key], (*path, key)))
            else:
                result.append(
                    {
                        "path": [*path, key],
                        "before_present": key in before,
                        "before": before.get(key),
                    }
                )
        return result
    if path == ("eos_records",):
        assert len(before) == len(after)
        return [
            change
            for i, (a, b) in enumerate(zip(before, after))
            for change in _changes(a, b, (*path, i))
        ]
    return [{"path": list(path), "before_present": True, "before": before}]


def restore_document(document: dict, binding: dict) -> dict:
    result = copy.deepcopy(document)
    for change in binding["changes"]:
        node = result
        for key in change["path"][:-1]:
            node = node[key]
        key = change["path"][-1]
        if change["before_present"]:
            node[key] = copy.deepcopy(change["before"])
        else:
            node.pop(key, None)
    for ordering in binding.get("object_orders", []):
        node = result
        for key in ordering["path"]:
            node = node[key]
        assert set(node) == set(ordering["keys"])
        values = {key: node[key] for key in ordering["keys"]}
        node.clear()
        node.update(values)
    return result


def prior_material_bytes(material: str, *, root: Path | None = None) -> bytes:
    """Recover the exact v4 predecessor for frozen historical audits."""
    relative = f"peritheos/data/materials/{material}.eosmat"
    data = ((ROOT if root is None else root) / relative).read_bytes()
    if not RECEIPT.exists():
        return data
    binding = json.loads(RECEIPT.read_bytes())["materials"].get(relative)
    # Later, separately reviewed additions retain their own exact v4/v5 pair.
    # Keep the original migration receipt immutable and verify that removing
    # only the named additive records still recovers its migrated predecessor.
    update_path = ROOT / "docs/data/garnet-angel-2022/v5-integration.json"
    if update_path.exists():
        update = json.loads(update_path.read_bytes())["materials"].get(relative)
        if update is not None:
            current = json.loads(data)
            if binding is not None:
                predecessor = copy.deepcopy(current)
                added = set(update["added_records"])
                present = {r["identifier"] for r in predecessor["eos_records"]}
                if not added <= present:
                    raise ValueError("Missing reviewed post-migration addition")
                predecessor["eos_records"] = [
                    r
                    for r in predecessor["eos_records"]
                    if r["identifier"] not in added
                ]
                if (
                    hashlib.sha256(format_eosmat(predecessor).encode()).hexdigest()
                    != binding["after_sha256"]
                ):
                    raise ValueError(
                        "Post-migration addition changed its historical predecessor"
                    )
            binding = update
    if binding is None:
        return data
    if hashlib.sha256(data).hexdigest() != binding["after_sha256"]:
        raise ValueError(f"EOSMAT v5 migration source changed: {relative}")
    old = format_eosmat(restore_document(json.loads(data), binding)).encode()
    if hashlib.sha256(old).hexdigest() != binding["before_sha256"]:
        raise ValueError(f"EOSMAT v5 migration predecessor changed: {relative}")
    return old


def migrate_catalog() -> None:
    from scripts.scientific_source_snapshot import RECORD_FIELDS, scientific_sha256

    if RECEIPT.exists():
        raise ValueError(
            "The v5 migration already has a receipt; do not overwrite history"
        )
    prepared = []
    proof = {
        "schema": "peritheos.eosmat-v5-migration.v1",
        "date": "2026-11-12",
        "scope": "Parameter-unit uncertainties only; coefficients, observations, confidence semantics, evidence and scientific classifications are unchanged.",
        "before_scientific_sha256": scientific_sha256(),
        "materials": {},
        "scientific_record_changes": {},
    }
    for path in sorted((ROOT / "peritheos/data/materials").glob("*.eosmat")):
        original = path.read_bytes()
        before = json.loads(original)
        if before["format_version"] != 4:
            raise ValueError(f"Expected the reviewed v4 catalog: {path}")
        after = convert_eosmat(before, version=5)
        changes = _changes(before, after)
        for change in changes:
            p = change["path"]
            allowed = (
                p == ["format_version"]
                or len(p) >= 3
                and p[0] == "eos_records"
                and (
                    p[2:] == ["parameter_error_volume_basis"]
                    or p[2:] == ["parameters", "V0", "error"]
                    or p[2:]
                    in (
                        ["parameter_covariance", "matrix"],
                        ["parameter_covariance", "volume_basis"],
                    )
                )
            )
            if not allowed:
                raise ValueError(f"Unexpected migration change: {path}: {p}")
        rendered = format_eosmat(after).encode()
        binding = {
            "before_sha256": hashlib.sha256(original).hexdigest(),
            "after_sha256": hashlib.sha256(rendered).hexdigest(),
            "changes": changes,
        }
        binding["object_orders"] = []
        for change in changes:
            if change["before_present"] and change["path"][-1] in {
                "parameter_error_volume_basis",
                "volume_basis",
            }:
                node = before
                for key in change["path"][:-1]:
                    node = node[key]
                binding["object_orders"].append(
                    {"path": change["path"][:-1], "keys": list(node)}
                )
        assert restore_document(after, binding) == before
        assert format_eosmat(before).encode() == original
        proof["materials"][path.relative_to(ROOT).as_posix()] = binding
        for old, new in zip(before["eos_records"], after["eos_records"]):
            a, b = ({key: r.get(key) for key in RECORD_FIELDS} for r in (old, new))
            if a != b:
                proof["scientific_record_changes"][old["identifier"]] = {
                    "before": a,
                    "after": b,
                }
        prepared.append((path, original, rendered))
    if any(path.read_bytes() != before for path, before, _ in prepared):
        raise ValueError("A source file changed while preparing the v5 migration")
    for path, _, rendered in prepared:
        path.write_bytes(rendered)
    proof["after_scientific_sha256"] = scientific_sha256()
    RECEIPT.write_text(json.dumps(proof, indent=2, ensure_ascii=False) + "\n")
    for path, _, _ in prepared:
        prior_material_bytes(path.stem)
        validate_eosmat_document(json.loads(path.read_bytes()))
    print(
        f"Migrated {len(prepared)} materials; all v4 predecessors recover byte for byte."
    )


if __name__ == "__main__":
    migrate_catalog()
