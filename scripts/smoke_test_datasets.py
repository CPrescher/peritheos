"""Verify the installed optional package outside the source checkout."""

import hashlib
import json
from importlib import resources
from importlib.metadata import version

from peritheos import (
    Dataset,
    get_material,
    get_material_document,
    get_study_dataset,
    get_study_document,
    observation_registry,
)
from peritheos.observations import (
    observation_material_fingerprint,
    observation_source_fingerprint,
)


def main():
    package = resources.files("peritheos_datasets")
    assert package.joinpath("eos-structure-rights-review.json").read_bytes() == (
        resources.files("peritheos")
        .joinpath("data", "eos-structure-rights-review.json")
        .read_bytes()
    )
    for notice in (
        "dioptas-materials.LICENSE.md",
        "eos-structure-compilations.LICENSE.md",
    ):
        assert package.joinpath("notices", notice).is_file()
    manifest = json.loads(package.joinpath("manifest.json").read_text(encoding="utf-8"))
    assert version("peritheos") == version("peritheos-datasets") == manifest["version"]
    review_bytes = package.joinpath("redistribution-review.json").read_bytes()
    assert (
        hashlib.sha256(review_bytes).hexdigest()
        == manifest["policy"]["redistribution_review"]["sha256"]
    )
    review = json.loads(review_bytes)
    reviewed = {
        e["review_key"]: e
        for e in review["entries"] + review.get("supporting_studies", [])
        if "review_key" in e
    }
    held_files = {
        e["resource"]["path"]
        for e in reviewed.values()
        if e["decision"] == "hold" and e["resource"]
    }
    assert held_files.isdisjoint(manifest["resources"])
    for relative, checksum in manifest["resources"].items():
        assert (
            hashlib.sha256(
                package.joinpath("data", *relative.split("/")).read_bytes()
            ).hexdigest()
            == checksum
        )
    for entry in manifest["policy"]["datasets"] + manifest["policy"].get(
        "study_datasets", []
    ):
        document = (
            get_study_document(entry["study"])
            if "study" in entry
            else get_material_document(entry["material"])
        )
        declaration = next(
            d for d in document["datasets"] if d["identifier"] == entry["dataset"]
        )
        if entry.get("review_key"):
            assessment = reviewed[entry["review_key"]]
            assert assessment["decision"] == "admit"
            assert assessment["basis"] == entry["redistribution_basis"]
            assert assessment["source_status"] == entry["source_status"]
            assert assessment["declaration_sha256"] == entry["declaration_sha256"]
            assert assessment["resource"] == declaration.get("resource")
        assert len(Dataset.from_mapping(declaration)) > 0
        # Package policy binds the complete declaration, including inline views.
        # Observation registry bindings below deliberately exclude those views.
        declaration_bytes = json.dumps(
            declaration,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode()
        assert (
            hashlib.sha256(declaration_bytes).hexdigest() == entry["declaration_sha256"]
        )
    for entry in manifest["policy"].get("study_datasets", []):
        assert len(get_study_dataset(entry["study"], entry["dataset"])) > 0
    assert len(get_material("mgo").get_dataset("zhu_2025_mgo_pvt")) > 0
    for entry in observation_registry()["datasets"]:
        document = get_material_document(entry["material_id"])
        dataset = next(
            d for d in document["datasets"] if d["identifier"] == entry["dataset_id"]
        )
        assert observation_source_fingerprint(dataset) == entry["declaration_sha256"]
        assert (
            observation_material_fingerprint(document)
            == entry["material_identity_sha256"]
        )
    print(
        f"Optional data smoke passed: {manifest['dataset_count']} material and {manifest.get('study_dataset_count', 0)} study declarations, {len(manifest['resources'])} resource files."
    )


if __name__ == "__main__":
    main()
