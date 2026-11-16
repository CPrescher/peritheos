"""Exact pre-dataset projections for frozen garnet audits (maintenance only)."""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from peritheos._eosmat_format import V4_LAYOUT, format_eosmat
from peritheos_build import digest

RECEIPT = "docs/data/garnet-angel-2022/dataset-integration-2026-11-15.json"


def prior_dataset_bytes(data: bytes, material: str, root: Path) -> bytes:
    """Reverse only the explicitly bound additions and prose updates."""
    path = root / RECEIPT
    if not path.exists():
        return data
    binding = json.loads(path.read_bytes())["materials"].get(material)
    if binding is None:
        return data
    actual = hashlib.sha256(data).hexdigest()
    if actual == binding["before_sha256"]:
        return data
    if actual != binding["after_sha256"]:
        raise ValueError(f"Garnet dataset integration source changed: {material}")
    document = json.loads(data)
    additions = binding["appended_datasets"]
    found = {
        d["identifier"]: d for d in document["datasets"] if d["identifier"] in additions
    }
    if set(found) != set(additions) or any(
        digest(found[k]) != v for k, v in additions.items()
    ):
        raise ValueError(f"Garnet dataset declaration changed: {material}")
    document["datasets"] = [
        d for d in document["datasets"] if d["identifier"] not in additions
    ]
    for field in binding["restored_fields"]:
        node = document
        for key in field["path"][:-1]:
            node = node[key]
        if field.get("before_present", True):
            node[field["path"][-1]] = copy.deepcopy(field["before"])
        else:
            node.pop(field["path"][-1], None)
    # The predecessor predates the v5 reading-order revision.
    previous = format_eosmat(document, layout=V4_LAYOUT).encode()
    if hashlib.sha256(previous).hexdigest() != binding["before_sha256"]:
        raise ValueError(f"Garnet dataset predecessor changed: {material}")
    return previous


def prior_dataset_snapshot(snapshot: dict, root: Path) -> dict | None:
    """Remove exact new dataset projections, retaining every old observation."""
    path = root / RECEIPT
    if not path.exists():
        return None
    receipt = json.loads(path.read_bytes())
    result = copy.deepcopy(snapshot)
    changed = False

    def fingerprint(value):
        return hashlib.sha256(
            json.dumps(
                value, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode()
        ).hexdigest()

    for document in result["documents"]:
        binding = receipt["materials"].get(document["identifier"])
        if binding is None:
            continue
        additions = binding["appended_datasets"]
        if not any(d["identifier"] in additions for d in document["datasets"]):
            continue
        if (
            fingerprint(document["datasets"])
            != binding["scientific_datasets_after_sha256"]
        ):
            raise ValueError("Garnet observation projection changed")
        document["datasets"] = [
            d for d in document["datasets"] if d["identifier"] not in additions
        ]
        if (
            fingerprint(document["datasets"])
            != binding["scientific_datasets_before_sha256"]
        ):
            raise ValueError("Garnet observation predecessor changed")
        changed = True
    return result if changed else None
