"""Repository-only snapshot of scientific catalog values and source arrays.

This deliberately excludes names, aliases and reconstruction narratives. A
reviewed fingerprint protects source values independently of the metadata-only
compatibility exceptions in merge_source_bindings.py.

``thermal`` is retained whole, including its free-text ``notes``. A later text
revision keeps the replaced original in the record's ``superseded_record_text``
sidecar block; the snapshot reads ``thermal.notes`` at that stored original, so
a text-only revision does not change the fingerprint.
"""

from __future__ import annotations

import copy
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


def _record_fields(record: dict, entry: dict) -> dict:
    """Snapshot fields of one record; ``thermal.notes`` at its stored original."""
    fields = {key: record.get(key) for key in RECORD_FIELDS}
    original = entry.get("superseded_record_text", {})
    if "thermal_notes" in original and isinstance(fields["thermal"], dict):
        thermal = copy.deepcopy(fields["thermal"])
        thermal.pop("notes", None)
        if original["thermal_notes"] not in (None, ""):
            thermal["notes"] = original["thermal_notes"]
        fields["thermal"] = thermal
    return fields


def scientific_snapshot(root: Path = ROOT) -> dict:
    """Retain literal scientific fields and the actual bytes of every resource."""
    documents = []
    for path in sorted((root / "peritheos/data/materials").glob("*.eosmat")):
        document = json.loads(path.read_text(encoding="utf-8"))
        audit = root / "peritheos/data/audits" / f"{path.stem}.audit.json"
        entries = (
            json.loads(audit.read_text(encoding="utf-8"))["records"]
            if audit.exists()
            else {}
        )
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
                    _record_fields(record, entries.get(record["identifier"], {}))
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


def assert_snapshot_lineage(snapshot: dict, historical_sha256: str) -> None:
    """Verify history across the explicitly authorized powder-pattern replacement.

    Original peak arrays remain private. Public receipts bind their digests,
    the exact new arrays, and unchanged scientific projections at each earlier
    stage. This never reports a historical digest as the current full digest.
    """
    import copy

    def fingerprint(value):
        return hashlib.sha256(
            json.dumps(
                value, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode()
        ).hexdigest()

    if fingerprint(snapshot) == historical_sha256:
        return
    from scripts.garnet_dataset_bindings import prior_dataset_snapshot

    previous = prior_dataset_snapshot(snapshot, ROOT)
    if previous is not None:
        assert_snapshot_lineage(previous, historical_sha256)
        return
    migration_path = ROOT / "docs/data/eosmat-v5-migration.json"
    if migration_path.exists():
        migration = json.loads(migration_path.read_bytes())
        changes = migration["scientific_record_changes"]
        projected = copy.deepcopy(snapshot)
        seen = set()
        for document in projected["documents"]:
            for index, record in enumerate(document["eos_records"]):
                identifier = record["identifier"]
                if identifier in changes and record == changes[identifier]["after"]:
                    document["eos_records"][index] = changes[identifier]["before"]
                    seen.add(identifier)
        if seen:
            assert seen == set(changes)
            assert_snapshot_lineage(projected, historical_sha256)
            return
    # A later curation may append exact, separately audited records/materials.
    # Verify their complete scientific projections before removing them to test
    # the historical invariant. Never replace a historical receipt with a new hash.
    addition_path = ROOT / "docs/data/garnet-angel-2022/integration-preservation.json"
    if addition_path.exists():
        additions = json.loads(addition_path.read_bytes())
        new_materials = additions["new_materials"]
        new_records = additions["new_records"]
        if any(
            d["identifier"] in new_materials
            or any(r["identifier"] in new_records for r in d["eos_records"])
            for d in snapshot["documents"]
        ):
            projected = copy.deepcopy(snapshot)
            kept = []
            seen_materials, seen_records = set(), set()
            for document in projected["documents"]:
                material = document["identifier"]
                if material in new_materials:
                    assert fingerprint(document) == new_materials[material]
                    seen_materials.add(material)
                    continue
                records = []
                for record in document["eos_records"]:
                    identifier = record["identifier"]
                    if identifier in new_records:
                        assert fingerprint(record) == new_records[identifier]
                        seen_records.add(identifier)
                    else:
                        records.append(record)
                document["eos_records"] = records
                kept.append(document)
            projected["documents"] = kept
            assert seen_materials == set(new_materials)
            assert seen_records == set(new_records)
            # Some callers have already removed later dataset additions to
            # verify an earlier historical stage. Its original receipt, below,
            # must validate that projected state instead of the latest baseline.
            if historical_sha256 == additions["before_scientific_sha256"]:
                assert fingerprint(projected) == additions["before_scientific_sha256"]
            assert_snapshot_lineage(projected, historical_sha256)
            return
    proof = json.loads(
        (ROOT / "docs/data/calculated-pattern-replacement.json").read_bytes()
    )
    assert historical_sha256 in proof["historical_snapshot_projections"]
    projected = copy.deepcopy(snapshot)
    seen = set()
    for document in projected["documents"]:
        material = document["identifier"]
        if material in proof["material_bindings"]:
            binding = proof["material_bindings"][material]
            # build.digest uses UTF-8 text; the scientific snapshot uses ASCII.
            from peritheos_build import digest

            assert digest(document["peaks"]) == binding["calculated_peaks_sha256"]
            document.pop("peaks")
            if material == "fe3s":
                assert (
                    digest(document["atom_sites"])
                    == binding["calculated_atom_sites_sha256"]
                )
                document.pop("atom_sites")
            seen.add(material)
    assert seen == set(proof["material_bindings"])
    assert (
        fingerprint(projected)
        == proof["historical_snapshot_projections"][historical_sha256]
    )
