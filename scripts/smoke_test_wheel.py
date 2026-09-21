#!/usr/bin/env python3
"""Check an installed wheel from a directory outside the source checkout."""

import json
from collections import Counter
from importlib import resources

from peritheos import (
    get_material_document,
    group_materials,
    list_eos_records,
    list_material_documents,
    list_material_families,
    list_materials,
)
from peritheos.eos.rt import BM3, Baonza


def main():
    package = resources.files("peritheos")
    manifest = json.loads(
        package.joinpath("data", "materials", "manifest.json").read_text()
    )
    documents = [
        get_material_document(identifier) for identifier in list_material_documents()
    ]
    assert package.joinpath("py.typed").is_file()
    assert package.joinpath("data", "material-families.json").is_file()
    assert list_material_families()
    assert package.joinpath(
        "data", "datasets", "experimental-metal-eos-batch-47.csv"
    ).is_file()
    assert BM3(10, 120, 4).pressure(9) > 0
    assert Baonza(10, 120, 4).pressure(9) > 0

    # Source documents include deferred records and cards with no executable EOS.
    records = [record for document in documents for record in document["eos_records"]]
    assert len(documents) == manifest["materials"]
    assert len(records) == manifest["eos_records"]
    assert (
        Counter(record["scientific_validation"]["status"] for record in records)
        == (manifest["scientific_validation"]["counts"])
    )

    expected = {}
    for document in documents:
        identifiers = {
            record["identifier"]
            for record in document["eos_records"]
            if record["scientific_validation"]["status"] == "primary_source_validated"
        }
        if identifiers:
            expected[document["identifier"]] = identifiers

    # Check identities per material so a deferred record cannot replace an
    # accepted one, or cause accepted records in the same material to disappear.
    materials = list_materials()
    assert {
        material.identifier: {record.identifier for record in material.eos_records}
        for material in materials
    } == expected
    executable_ids = {identifier for ids in expected.values() for identifier in ids}
    assert {record.identifier for record in list_eos_records()} == executable_ids
    assert (
        sum(len(material.eos_records) for material in materials)
        == (manifest["scientific_validation"]["counts"]["primary_source_validated"])
    )
    assert sum(len(group.materials) for group in group_materials(materials)) == len(
        materials
    )
    print(
        f"Wheel smoke passed: {len(documents)} documents, {len(records)} source "
        f"records, {len(materials)} executable materials, "
        f"{len(executable_ids)} executable records."
    )


if __name__ == "__main__":
    main()
