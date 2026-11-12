"""Recover hash-bound pre-annotation documents for historical audits only.

Runtime readers always retain the explicit unit annotations. This projection
does not erase any coefficient, observation, uncertainty or source evidence.
"""

from __future__ import annotations

import copy
import hashlib
import json
from pathlib import Path

from peritheos._eosmat_format import format_eosmat
from peritheos_build import digest
from scripts.migrate_eosmat_v5 import prior_material_bytes as before_v5_bytes

ROOT = Path(__file__).resolve().parents[1]
RECEIPT = ROOT / "docs/data/uncertainty-volume-basis-review.json"


def _strip(record: dict) -> None:
    record.pop("parameter_error_volume_basis", None)
    covariance = record.get("parameter_covariance")
    if isinstance(covariance, dict):
        covariance.pop("volume_basis", None)


def prior_material_bytes(material: str) -> bytes:
    """Restore exact prior bytes, checking both the current and prior hashes."""
    relative = f"peritheos/data/materials/{material}.eosmat"
    data = before_v5_bytes(material, root=ROOT)
    binding = json.loads(RECEIPT.read_bytes())["materials"].get(relative)
    if binding is None:
        return data
    if hashlib.sha256(data).hexdigest() != binding["after_sha256"]:
        raise ValueError(f"Uncertainty-unit review source changed: {relative}")
    document = json.loads(data)
    for record in document["eos_records"]:
        _strip(record)
    if digest(document) != binding["preserved_document_sha256"]:
        raise ValueError(
            f"Uncertainty-unit review scientific content changed: {relative}"
        )
    original = format_eosmat(document).encode("utf-8")
    if hashlib.sha256(original).hexdigest() != binding["before_sha256"]:
        raise ValueError(f"Uncertainty-unit review predecessor changed: {relative}")
    return original


def prior_material_document(material: str) -> dict:
    """Load the hash-verified predecessor for a frozen historical comparison."""
    return json.loads(prior_material_bytes(material))


def historical_unit_view(document: dict) -> dict:
    """Project only exact reviewed legacy records to their earlier fingerprints."""
    reviewed = {r["identifier"]: r for r in json.loads(RECEIPT.read_bytes())["records"]}
    result = copy.deepcopy(document)
    if document.get("format_version") == 3 and document.get("identifier"):
        from peritheos import convert_eosmat

        material = document["identifier"]
        path = ROOT / f"peritheos/data/materials/{material}.eosmat"
        if path.is_file() and json.loads(path.read_bytes()).get("format_version") == 5:
            current = convert_eosmat(
                json.loads(path.read_bytes()), version=3, legacy_uncertainty_units=False
            )
            previous = convert_eosmat(
                json.loads(before_v5_bytes(material)),
                version=3,
                legacy_uncertainty_units=False,
            )
            old_records = {r["identifier"]: r for r in previous["eos_records"]}
            current_records = {r["identifier"]: r for r in current["eos_records"]}
            for i, record in enumerate(result.get("eos_records", [])):
                if record == current_records.get(record.get("identifier")):
                    result["eos_records"][i] = old_records[record["identifier"]]
    for record in result.get("eos_records", []):
        binding = reviewed.get(record.get("identifier"))
        if binding is None:
            continue
        original = copy.deepcopy(record)
        _strip(original)
        if digest(original) == binding["preserved_legacy_record_sha256"]:
            record.clear()
            record.update(original)
    return result
