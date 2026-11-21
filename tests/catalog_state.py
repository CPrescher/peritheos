"""Catalog totals for tests, anchored independently of the files they check.

``BASE`` holds the totals of the bundled catalog as curated before the
catalog revision ledger (``docs/data/catalog-revisions.json``), as literals.
The expected current totals are ``BASE`` plus the effect of the ledger's
revisions, computed from the ledger alone: each added record's stored content
(headline, statuses, calibration) and each acknowledged headline change or
refit-row change. ``tests/test_catalog_manifest.py`` asserts that the actual
catalog, the manifest and these expectations agree, so a record lost or
duplicated without a revision fails there even when the manifest was
rewritten to match.
"""

from __future__ import annotations

import copy
import json
from collections import Counter
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]

#: Totals before the revision ledger (commit 5dab6977).
BASE: dict[str, Any] = {
    "materials": 235,
    "records": 687,
    "validation": {
        "primary_source_validated": 666,
        "not_reproduced": 18,
        "deferred": 3,
    },
    "source_check": {"verified": 666, "unresolved": 18, "deferred": 3},
    "reproduction": {"reproduced": 358, "close": 72, "disagrees": 61, "untested": 196},
    "compared": {
        "published_fit": 474,
        "peritheos_refit": 16,
        "source_equation": 19,
        "none": 178,
    },
    "with_warning": 34,
    "refit_status": {
        "not_refittable": 178,
        "parity": 358,
        "parity_not_achieved": 61,
        "similar": 72,
        "source_reconstruction": 18,
    },
    "calibration_status": {
        "not_applicable": 164,
        "partially_resolved": 193,
        "resolved": 242,
        "unresolved": 88,
    },
    "recalculation": {
        "missing_calibrant_observations": 200,
        "not_applicable": 161,
        "not_possible": 127,
        "ready": 136,
        "reference_eos_not_bundled": 56,
        "reference_model_not_supported": 7,
    },
    "eos_source_groups": 273,
    "numerical_parameter_assessment": 680,
}
#: Source-check codes of scientific_validation.status (EOS labels).
SOURCE_CHECK = {
    "primary_source_validated": "verified",
    "not_reproduced": "unresolved",
    "deferred": "deferred",
}


def manifest(root: Path = ROOT) -> dict[str, Any]:
    """The bundled material manifest."""
    return json.loads(
        (root / "peritheos/data/materials/manifest.json").read_text(encoding="utf-8")
    )


def _revisions(root: Path = ROOT) -> list[dict[str, Any]]:
    path = root / "docs/data/catalog-revisions.json"
    if not path.exists():
        return []
    return list(json.loads(path.read_text(encoding="utf-8"))["revisions"])


def added_records(root: Path = ROOT) -> set[str]:
    """Records added by catalog revisions (absent from the pre-ledger catalog)."""
    return set(_added(root))


def _added(root: Path = ROOT) -> dict[str, dict[str, Any]]:
    return {
        operation["record"]["identifier"]: operation["record"]
        for revision in _revisions(root)
        for operation in revision["forward"]
        if operation["op"] == "add_record"
    }


def revised_records(root: Path = ROOT) -> set[str]:
    """Existing records whose fields a catalog revision set."""
    return {
        operation["record"]
        for revision in _revisions(root)
        for operation in revision["forward"]
        if operation["op"] == "set"
    }


def expected(root: Path = ROOT) -> dict[str, Any]:
    """``BASE`` plus the effect of every revision in the ledger."""
    result = copy.deepcopy(BASE)

    def bump(key: str, value: Any, step: int = 1) -> None:
        if value is None:
            return
        counts = result[key]
        counts[value] = counts.get(value, 0) + step
        if counts[value] == 0:
            del counts[value]

    thermodynamic = json.loads(
        (root / "docs/data/eos-structure-rights-policy.json").read_text(
            encoding="utf-8"
        )
    )["thermodynamic_collection_doi"]
    for record in _added(root).values():
        result["records"] += 1
        status = record["evidence"]["scientific_validation"]["status"]
        bump("validation", status)
        bump("source_check", SOURCE_CHECK.get(status))
        bump("reproduction", record["check"]["reproduction"])
        bump("compared", record["check"]["compared"])
        result["with_warning"] += record["check"].get("warning") is not None
        calibration = record["pressure_calibration"]
        bump("calibration_status", calibration["status"])
        bump("recalculation", calibration["recalculation"]["status"])
        if (record.get("reference") or {}).get("doi") != thermodynamic:
            result["numerical_parameter_assessment"] += 1
    for revision in _revisions(root):
        for fields in revision.get("acknowledge_check_change", {}).values():
            for field, change in fields.items():
                if field == "warning":
                    result["with_warning"] += (change["to"] is not None) - (
                        change["from"] is not None
                    )
                elif field in ("reproduction", "compared"):
                    bump(field, change["from"], -1)
                    bump(field, change["to"])
                elif field == "source_check":
                    bump("source_check", change["from"], -1)
                    bump("source_check", change["to"])
        for pair in revision.get("refit_rows", {}).values():
            if pair["before"] is not None:
                bump("refit_status", pair["before"]["status"], -1)
            if pair["after"] is not None:
                bump("refit_status", pair["after"]["status"])
    return result


def record_count(root: Path = ROOT) -> int:
    return int(expected(root)["records"])


def material_count(root: Path = ROOT) -> int:
    return int(expected(root)["materials"])


def validation_count(status: str, root: Path = ROOT) -> int:
    return int(expected(root)["validation"].get(status, 0))


def catalog_counts(root: Path = ROOT) -> dict[str, Any]:
    """Totals read directly from the record files and sidecars under ``root``."""
    materials = 0
    records: list[tuple[dict, dict]] = []
    for path in sorted((root / "peritheos/data/materials").glob("*.eosmat")):
        materials += 1
        document = json.loads(path.read_text(encoding="utf-8"))
        sidecar = json.loads(
            (root / f"peritheos/data/audits/{path.stem}.audit.json").read_text(
                encoding="utf-8"
            )
        )["records"]
        for record in document["eos_records"]:
            records.append((record, sidecar.get(record["identifier"], {})))
    return {
        "materials": materials,
        "records": len(records),
        "identifiers": len({r["identifier"] for r, _ in records}),
        "validation": dict(
            Counter(
                e["evidence"]["scientific_validation"]["status"] for _, e in records
            )
        ),
        "reproduction": dict(Counter(r["check"]["reproduction"] for r, _ in records)),
        "compared": dict(Counter(r["check"]["compared"] for r, _ in records)),
        "source_check": dict(
            Counter(r["check"].get("source_check") for r, _ in records)
        ),
        "with_warning": sum(r["check"].get("warning") is not None for r, _ in records),
        "calibration_status": dict(
            Counter(e["pressure_calibration"]["status"] for _, e in records)
        ),
        "recalculation": dict(
            Counter(
                e["pressure_calibration"]["recalculation"]["status"] for _, e in records
            )
        ),
    }
