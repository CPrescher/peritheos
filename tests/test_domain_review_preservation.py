"""The domain review freezes scientific values independently of metadata pins."""

import copy
import hashlib
import json
from pathlib import Path

from scripts.merge_source_bindings import accepts_merge_binding
from scripts.scientific_source_snapshot import scientific_sha256, scientific_snapshot

ROOT = Path(__file__).resolve().parents[1]


def test_reviewed_coefficients_units_and_source_arrays_remain_exact():
    review = json.loads(
        (ROOT / "docs/data/pr27-scientific-domain-review.json").read_text(
            encoding="utf-8"
        )
    )
    replacement = json.loads(
        (ROOT / "docs/data/brucite-structure-replacement.json").read_bytes()
    )
    assert (
        review["catalog_scientific_sha256"]
        == replacement["before_catalog_scientific_sha256"]
    )
    followup = json.loads(
        (ROOT / "docs/data/morard-campbell-redistribution-followup.json").read_bytes()
    )
    assert (
        followup["before_catalog_scientific_sha256"]
        == replacement["after_catalog_scientific_sha256"]
    )
    assert scientific_sha256() == followup["after_catalog_scientific_sha256"]
    snapshot = copy.deepcopy(scientific_snapshot())
    added = set(followup["added_dataset_declarations"])
    removed = set()
    for document in snapshot["documents"]:
        retained = []
        for dataset in document["datasets"]:
            key = document["identifier"] + "/" + dataset["identifier"]
            if key in added:
                removed.add(key)
            else:
                retained.append(dataset)
        document["datasets"] = retained
    assert removed == added
    assert (
        hashlib.sha256(
            json.dumps(
                snapshot, sort_keys=True, separators=(",", ":"), allow_nan=False
            ).encode()
        ).hexdigest()
        == replacement["after_catalog_scientific_sha256"]
    )
    for document in snapshot["documents"]:
        if document["identifier"] == "brucite":
            for key in ("lattice", "atom_sites"):
                document.pop(key)
    fingerprint = hashlib.sha256(
        json.dumps(
            snapshot, sort_keys=True, separators=(",", ":"), allow_nan=False
        ).encode()
    ).hexdigest()
    assert (
        fingerprint
        == replacement[
            "unchanged_catalog_excluding_brucite_lattice_and_atom_sites_sha256"
        ]
    )


def test_structure_replacement_accepts_only_the_exact_reviewed_pair():
    replacement = json.loads(
        (ROOT / "docs/data/brucite-structure-replacement.json").read_bytes()
    )
    context = "release_document"
    name = "peritheos/data/materials/brucite.eosmat"
    binding = replacement["bindings"][context][name]
    old, new = binding["historical_sha256"], binding["replacement_sha256"]
    assert accepts_merge_binding(context, name, old, new)
    assert not accepts_merge_binding(context, name, "0" * 64, new)
    assert not accepts_merge_binding(context, name, old, "0" * 64)
    assert not accepts_merge_binding("unreviewed", name, old, new)
