#!/usr/bin/env python3
"""Audit every stored EOS and its observation views without altering the catalog.

Run: python scripts/audit_observation_coverage.py --output /tmp/coverage.json
The output reports current checkout contents, not an assumed release revision.
"""

from __future__ import annotations

import argparse
import json
import sys
from collections import Counter
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

import numpy as np  # noqa: E402

from peritheos.datasets import Dataset  # noqa: E402
from peritheos.observations import (  # noqa: E402
    observation_registry,
    plan_dataset_observations,
    resolve_dataset_observations,
)


def audit_observation_coverage(materials: Path, *, resolve: bool = True) -> dict:
    registry = observation_registry()
    records = []
    material_paths = sorted(materials.glob("*.eosmat"))
    for path in material_paths:
        document = json.loads(path.read_text())
        for record in document["eos_records"]:
            identifier = record["identifier"]
            audit = record.get("scientific_validation", {}).get(
                "primary_data_check", {}
            )
            audit_ids = set(audit.get("dataset_identifiers", [])) | set(
                audit.get("digitized_dataset_identifiers", [])
            )
            declared = {
                e["dataset_id"]
                for e in registry["datasets"]
                if e["material_id"] == document["identifier"]
                and any(identifier in v["eos_records"] for v in e["views"])
            }
            linked = [
                d
                for d in document.get("datasets", [])
                if identifier in d.get("used_by_eos_records", [])
                or identifier
                in {r for v in d.get("observation_views", []) for r in v["eos_records"]}
                or d["identifier"] in audit_ids | declared
            ]
            tables = []
            for dataset in linked:
                entry = {
                    "dataset_id": dataset["identifier"],
                    "source_location": dataset["source_location"],
                }
                try:
                    raw = Dataset.from_mapping(dataset)
                    native = []
                    for column in raw.columns:
                        if column.role != "value":
                            continue
                        count = 0
                        for value in raw[column.name]:
                            try:
                                count += int(np.isfinite(float(value)))
                            except (ValueError, TypeError):
                                pass
                        if count:
                            native.append(
                                {
                                    "column": column.name,
                                    "unit": column.unit,
                                    "quantity": column.quantity,
                                    "numeric_values": count,
                                }
                            )
                    entry["native_observables"] = native
                except Exception as error:
                    entry["source_access_error"] = str(error)
                try:
                    plan = plan_dataset_observations(
                        document, identifier, dataset["identifier"], registry=registry
                    )
                    entry.update(plan)
                    if plan["status"] == "ready" and resolve:
                        data = resolve_dataset_observations(
                            document,
                            identifier,
                            dataset["identifier"],
                            registry=registry,
                        )
                        entry.update(
                            points=len(data.row_indices),
                            omitted_rows=len(data.omissions),
                            selected_rows=len(data.row_indices) + len(data.omissions),
                            total_source_rows=len(data.dataset),
                        )
                        if len(data.row_indices) == 0:
                            entry.update(
                                status="requires_review",
                                reason="no_finite_pressure_volume_pairs",
                            )
                except Exception as error:
                    entry.update(
                        status="requires_review",
                        reason="invalid_or_stale_view",
                        detail=str(error),
                    )
                tables.append(entry)
            measured = [
                t
                for t in tables
                if t["status"] == "ready"
                and t["view"]["evidence_kind"]
                in {"experimental", "digitized", "experimental_reduction", "shock"}
            ]
            records.append(
                {
                    "material_id": document["identifier"],
                    "eos_record": identifier,
                    "determination_method": record.get(
                        "determination_method", "unknown"
                    ),
                    "experimental_data_status": "plottable"
                    if measured
                    else ("no_linked_data" if not tables else "requires_review"),
                    "primary_data_status": audit.get("status"),
                    "primary_data_finding": audit.get("finding"),
                    "missing_referenced_datasets": sorted(
                        audit_ids
                        - {d["identifier"] for d in document.get("datasets", [])}
                    ),
                    "tables": tables,
                }
            )
    experimental = [
        r for r in records if r["determination_method"] in {"experimental", "hybrid"}
    ]
    return {
        "format": "peritheos.observation-coverage",
        "format_version": 1,
        "scope": "Current source checkout; display-coordinate resolution, not publication-fit reproduction.",
        "counts": {
            "materials": len(material_paths),
            "records": len(records),
            "experimental_or_hybrid": len(experimental),
            "experimental_data_status": dict(
                Counter(r["experimental_data_status"] for r in experimental)
            ),
            "unresolved_reasons": dict(
                Counter(
                    t.get("reason")
                    for r in experimental
                    for t in r["tables"]
                    if t["status"] != "ready"
                )
            ),
        },
        "records": records,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--materials",
        type=Path,
        default=Path(__file__).resolve().parents[1] / "peritheos/data/materials",
    )
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--plans-only", action="store_true")
    args = parser.parse_args()
    report = audit_observation_coverage(args.materials, resolve=not args.plans_only)
    args.output.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n")
    print(json.dumps(report["counts"], indent=2))


if __name__ == "__main__":
    main()
