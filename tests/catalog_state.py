"""Catalog totals for tests, from one place.

Pure catalog counts come from ``peritheos/data/materials/manifest.json``;
``tests/test_catalog_manifest.py`` checks the manifest against the catalog,
so a lost or extra record still fails there. Records added later through
``scripts/apply_catalog_revision.py`` are named by its ledger, so tests that
describe the catalog as curated before the ledger can leave them out.
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]


def manifest() -> dict[str, Any]:
    """The bundled material manifest."""
    return json.loads(
        (ROOT / "peritheos/data/materials/manifest.json").read_text(encoding="utf-8")
    )


def record_count() -> int:
    return int(manifest()["eos_records"])


def material_count() -> int:
    return int(manifest()["materials"])


def validation_count(status: str) -> int:
    return int(manifest()["scientific_validation"]["counts"].get(status, 0))


def _revisions() -> list[dict[str, Any]]:
    path = ROOT / "docs/data/catalog-revisions.json"
    if not path.exists():
        return []
    return list(json.loads(path.read_text(encoding="utf-8"))["revisions"])


def added_records() -> set[str]:
    """Records added by catalog revisions (absent from the pre-ledger catalog)."""
    return {
        operation["record"]["identifier"]
        for revision in _revisions()
        for operation in revision["forward"]
        if operation["op"] == "add_record"
    }


def revised_records() -> set[str]:
    """Existing records whose fields a catalog revision set."""
    return {
        operation["record"]
        for revision in _revisions()
        for operation in revision["forward"]
        if operation["op"] == "set"
    }
