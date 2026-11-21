"""The catalog, its manifest and the ledger-anchored expectation agree."""

import json
import shutil
from pathlib import Path

from tests.catalog_state import ROOT, catalog_counts, expected, manifest

KEYS = (
    "materials",
    "records",
    "validation",
    "reproduction",
    "compared",
    "source_check",
    "with_warning",
    "calibration_status",
    "recalculation",
)


def mismatches(root: Path) -> list[str]:
    """Totals of the catalog or manifest under ``root`` that differ from BASE
    plus the ledger's effect."""
    actual = catalog_counts(root)
    want = expected()
    data = manifest(root)
    found = [key for key in KEYS if actual[key] != want[key]]
    if actual["identifiers"] != actual["records"]:
        found.append("duplicate identifiers")
    if data["materials"] != want["materials"]:
        found.append("manifest materials")
    if data["eos_records"] != want["records"]:
        found.append("manifest eos_records")
    if data["scientific_validation"]["counts"] != dict(
        sorted(want["validation"].items())
    ):
        found.append("manifest validation counts")
    calibration = data["pressure_calibration"]
    if calibration["status_counts"] != {
        key: want["calibration_status"].get(key, 0)
        for key in calibration["status_counts"]
    }:
        found.append("manifest calibration counts")
    if calibration["recalculation_counts"] != dict(
        sorted(want["recalculation"].items())
    ):
        found.append("manifest recalculation counts")
    return found


def test_catalog_manifest_and_ledger_anchored_totals_agree():
    assert mismatches(ROOT) == []


def test_a_record_lost_or_duplicated_without_a_revision_is_caught(tmp_path):
    for relative in ("peritheos/data/materials", "peritheos/data/audits"):
        shutil.copytree(ROOT / relative, tmp_path / relative)
    assert mismatches(tmp_path) == []
    path = tmp_path / "peritheos/data/materials/qandilite.eosmat"
    document = json.loads(path.read_bytes())
    manifest_path = tmp_path / "peritheos/data/materials/manifest.json"
    # Lost, with the manifest rewritten to count what is left.
    path.write_text(json.dumps(dict(document, eos_records=document["eos_records"][1:])))
    rewritten = manifest(tmp_path)
    rewritten["eos_records"] = catalog_counts(tmp_path)["records"]
    manifest_path.write_text(json.dumps(rewritten))
    assert "records" in mismatches(tmp_path)
    # Duplicated: the same record twice.
    records = [*document["eos_records"], document["eos_records"][0]]
    path.write_text(json.dumps(dict(document, eos_records=records)))
    found = mismatches(tmp_path)
    assert "records" in found and "duplicate identifiers" in found
