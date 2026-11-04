"""Repository-only snapshot of scientific catalog values and source arrays.

This deliberately excludes names, aliases and reconstruction narratives. A
reviewed fingerprint protects source values independently of the metadata-only
compatibility exceptions in merge_source_bindings.py.
"""

from __future__ import annotations

import hashlib
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
RECORD_FIELDS = (
    "identifier",
    "model",
    "parameters",
    "thermal",
    "temperature_ref",
    "parameter_error_confidence",
    "experimental_pressure_range_gpa",
    "experimental_temperature_range_k",
    "record_kind",
    "initial_state",
    "loading_path",
    "branch_kind",
    "pressure_scale_variant_of",
)


def scientific_snapshot(root: Path = ROOT) -> dict:
    """Retain literal scientific fields and the actual bytes of every resource."""
    documents = []
    for path in sorted((root / "peritheos/data/materials").glob("*.eosmat")):
        document = json.loads(path.read_text(encoding="utf-8"))
        datasets = []
        for dataset in document.get("datasets", []):
            entry = {
                key: dataset.get(key)
                for key in (
                    "identifier",
                    "columns",
                    "rows",
                    "resource",
                    "observation_views",
                    "uncertainty",
                    "uncertainty_convention",
                    "quantity_origins",
                    "provenance_only",
                    "usage_recommendation",
                    "used_by_eos_records",
                )
            }
            resource = dataset.get("resource")
            if resource:
                source = root / "peritheos/data" / resource["path"]
                entry["resource_bytes_sha256"] = hashlib.sha256(
                    source.read_bytes()
                ).hexdigest()
            datasets.append(entry)
        documents.append(
            {
                **{
                    key: document.get(key)
                    for key in (
                        "identifier",
                        "formula",
                        "units",
                        "formula_units_per_cell",
                        "phase",
                        "symmetry",
                        "lattice",
                        "space_group",
                        "space_group_number",
                        "atom_sites",
                        "peaks",
                    )
                },
                "eos_records": [
                    {key: record.get(key) for key in RECORD_FIELDS}
                    for record in document["eos_records"]
                ],
                "datasets": datasets,
            }
        )
    return {"documents": documents}


def scientific_sha256(root: Path = ROOT) -> str:
    payload = json.dumps(
        scientific_snapshot(root),
        sort_keys=True,
        separators=(",", ":"),
        allow_nan=False,
    )
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


if __name__ == "__main__":
    print(scientific_sha256())
