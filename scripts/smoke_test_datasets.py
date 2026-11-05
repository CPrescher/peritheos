"""Verify the installed optional package outside the source checkout."""

import hashlib
import json
from importlib import resources
from importlib.metadata import version

from peritheos import Dataset, get_material, get_material_document, observation_registry
from peritheos.observations import (
    observation_material_fingerprint,
    observation_source_fingerprint,
)


def main():
    package = resources.files("peritheos_datasets")
    manifest = json.loads(package.joinpath("manifest.json").read_text(encoding="utf-8"))
    assert version("peritheos") == version("peritheos-datasets") == manifest["version"]
    for relative, checksum in manifest["resources"].items():
        assert (
            hashlib.sha256(
                package.joinpath("data", *relative.split("/")).read_bytes()
            ).hexdigest()
            == checksum
        )
    for entry in manifest["policy"]["datasets"]:
        document = get_material_document(entry["material"])
        declaration = next(
            d for d in document["datasets"] if d["identifier"] == entry["dataset"]
        )
        assert len(Dataset.from_mapping(declaration)) > 0
        assert (
            observation_source_fingerprint(declaration) == entry["declaration_sha256"]
        )
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
        f"Optional data smoke passed: {manifest['dataset_count']} exact dataset declarations, {len(manifest['resources'])} resource files."
    )


if __name__ == "__main__":
    main()
