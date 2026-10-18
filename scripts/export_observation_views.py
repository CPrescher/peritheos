#!/usr/bin/env python3
"""Export source-bound observation metadata for an immutable material revision.

The views may come from a local unreleased change. The manifest records their
content hash separately from the scientific source revision; it never repins EOS.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import subprocess
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

from peritheos.eosmat import convert_eosmat  # noqa: E402
from peritheos.observations import (  # noqa: E402
    observation_material_fingerprint,
    observation_source_fingerprint,
    validate_observation_view,
)


def export_observation_views(revision: str) -> dict:
    def git(*args: str) -> bytes:
        return subprocess.check_output(["git", "-C", str(ROOT), *args])

    commit = git("rev-parse", "--verify", revision + "^{commit}").decode().strip()
    registry_path = "peritheos/data/observation-views.json"
    archived = (
        subprocess.run(
            ["git", "-C", str(ROOT), "cat-file", "-e", commit + ":" + registry_path],
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        ).returncode
        == 0
    )
    registry_bytes = (
        git("show", commit + ":" + registry_path)
        if archived
        else (ROOT / registry_path).read_bytes()
    )
    registry = json.loads(registry_bytes)
    entries = []
    rejected = []
    for entry in registry["datasets"]:
        path = "peritheos/data/materials/" + entry["material_id"] + ".eosmat"
        try:
            document = convert_eosmat(
                json.loads(git("show", commit + ":" + path)), version=3
            )
        except subprocess.CalledProcessError:
            rejected.append(
                {
                    "material_id": entry["material_id"],
                    "dataset_id": entry["dataset_id"],
                    "reason": "material_absent",
                }
            )
            continue
        dataset = next(
            (
                d
                for d in document.get("datasets", [])
                if d["identifier"] == entry["dataset_id"]
            ),
            None,
        )
        if (
            dataset is None
            or observation_source_fingerprint(dataset) != entry["declaration_sha256"]
            or observation_material_fingerprint(document)
            != entry["material_identity_sha256"]
        ):
            rejected.append(
                {
                    "material_id": entry["material_id"],
                    "dataset_id": entry["dataset_id"],
                    "reason": "source_binding_differs",
                }
            )
            continue
        records = {r["identifier"] for r in document["eos_records"]}
        calibrations = {
            r["identifier"]: r.get("pressure_calibration")
            for r in document["eos_records"]
        }
        if any(
            entry["record_pressure_calibration_sha256"].get(rid)
            != observation_source_fingerprint(
                {"pressure_calibration": calibrations.get(rid)}
            )
            for v in entry["views"]
            for rid in v["eos_records"]
        ):
            rejected.append(
                {
                    "material_id": entry["material_id"],
                    "dataset_id": entry["dataset_id"],
                    "reason": "record_pressure_calibration_differs",
                }
            )
            continue
        entry["record_pressure_calibrations"] = {
            rid: calibrations[rid]
            for view in entry["views"]
            for rid in view["eos_records"]
        }
        # An upstream record deletion must not remove a still-valid archived view.
        for view in entry["views"]:
            validate_observation_view(view, dataset, records)
        entries.append(entry)
    # New contributions can declare views directly on datasets; a sidecar entry
    # must not be required for new materials to reach downstream applications.
    paths = (
        git("ls-tree", "-r", "--name-only", commit, "peritheos/data/materials")
        .decode()
        .splitlines()
    )
    for path in paths:
        if not path.endswith(".eosmat"):
            continue
        document = convert_eosmat(
            json.loads(git("show", commit + ":" + path)), version=3
        )
        records = {r["identifier"]: r for r in document.get("eos_records", [])}
        for dataset in document.get("datasets", []):
            views = dataset.get("observation_views", [])
            if not views:
                continue
            for view in views:
                validate_observation_view(view, dataset, set(records))
            existing = next(
                (
                    e
                    for e in entries
                    if e["material_id"] == document["identifier"]
                    and e["dataset_id"] == dataset["identifier"]
                ),
                None,
            )
            if existing and {
                r for v in existing["views"] for r in v["eos_records"]
            }.intersection(r for v in views for r in v["eos_records"]):
                raise ValueError("Duplicate embedded and registry observation views")
            if existing is None:
                existing = {
                    "material_id": document["identifier"],
                    "dataset_id": dataset["identifier"],
                    "declaration_sha256": observation_source_fingerprint(dataset),
                    "material_identity_sha256": observation_material_fingerprint(
                        document
                    ),
                    "views": [],
                    "record_pressure_calibration_sha256": {},
                    "record_pressure_calibrations": {},
                }
                entries.append(existing)
            existing["views"].extend(views)
            for view in views:
                for identifier in view["eos_records"]:
                    calibration = records[identifier].get("pressure_calibration")
                    existing["record_pressure_calibrations"][identifier] = calibration
                    existing["record_pressure_calibration_sha256"][identifier] = (
                        observation_source_fingerprint(
                            {"pressure_calibration": calibration}
                        )
                    )
    return {
        "format": "peritheos.observation-view-export",
        "format_version": 1,
        "scientific_revision": commit,
        "view_registry_sha256": hashlib.sha256(registry_bytes).hexdigest(),
        "view_registry_origin": (
            "immutable Peritheos commit " + commit
            if archived
            else "explicit local Peritheos metadata export; not a scientific release"
        ),
        "datasets": entries,
        "incompatible_datasets": rejected,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--revision", required=True)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = (
        json.dumps(
            export_observation_views(args.revision), ensure_ascii=False, indent=2
        )
        + "\n"
    )
    if args.output:
        args.output.write_text(result, encoding="utf-8")
    else:
        sys.stdout.write(result)


if __name__ == "__main__":
    main()
