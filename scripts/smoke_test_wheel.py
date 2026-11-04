#!/usr/bin/env python3
"""Check an installed wheel from a directory outside the source checkout."""

import json
from collections import Counter
from importlib import resources
from importlib.metadata import distribution

from peritheos import (
    get_material_document,
    group_materials,
    list_eos_records,
    list_material_documents,
    list_material_families,
    list_materials,
)
from peritheos.eos.rt import BM3, Baonza, BM2Reference


def main():
    package = resources.files("peritheos")
    manifest = json.loads(
        package.joinpath("data", "materials", "manifest.json").read_text(
            encoding="utf-8"
        )
    )
    documents = [
        get_material_document(identifier) for identifier in list_material_documents()
    ]
    assert package.joinpath("py.typed").is_file()
    # License metadata must retain the installed data notices, rather than
    # presenting the mixed code/data distribution as uniformly MIT.
    installed = distribution("peritheos")
    assert (
        installed.metadata["License-Expression"] == "LicenseRef-Peritheos-Distribution"
    )
    license_files = set(installed.metadata.get_all("License-File", []))
    assert {
        "LICENSE",
        "LICENSES/CC-BY-4.0.txt",
        "LICENSES/CC0-1.0.txt",
        "LICENSES/LicenseRef-Peritheos-Distribution.txt",
        "peritheos/data/LICENSING.md",
    } <= license_files
    installed_licenses = {
        str(path).split(".dist-info/licenses/", 1)[1]: installed.locate_file(path)
        for path in installed.files or []
        if ".dist-info/licenses/" in str(path)
    }
    assert license_files == set(installed_licenses)
    for path in installed_licenses.values():
        assert path.is_file() and path.stat().st_size > 0
    data = package.joinpath("data")
    assert not data.joinpath(
        "datasets", "ma_2024_sources", "supporting-information.docx"
    ).is_file(), "The Ma publisher supplement must not be distributed"
    for source in (
        ("kcl_variant_sources", "6779TableS1 revised.xlsx"),
        ("chidester_2018_sources", "AM-18-56212.zip"),
    ):
        assert not data.joinpath("datasets", *source).is_file(), source
    # Removal of containers must not remove the curated numerical resources.
    for resource in (
        "kcl-tateno-2019-official-table-s1.csv",
        "chidester_2018_sources/6212_supp1.csv",
        "chidester_2018_sources/6212_supp2.csv",
    ):
        assert data.joinpath("datasets", resource).is_file(), resource
    guide = data.joinpath("LICENSING.md").read_text(encoding="utf-8")
    datasets = data.joinpath("datasets")
    notices = [path for path in datasets.iterdir() if path.name.endswith(".LICENSE.md")]
    notices.append(datasets.joinpath("ma_2024_sources", "LICENSE.md"))
    for notice in notices:
        relative = (
            "datasets/ma_2024_sources/LICENSE.md"
            if notice.name == "LICENSE.md"
            else f"datasets/{notice.name}"
        )
        assert f"]({relative})" in guide
        assert (
            notice.read_bytes()
            == installed_licenses[f"peritheos/data/{relative}"].read_bytes()
        )
    assert package.joinpath("data", "material-families.json").is_file()
    assert list_material_families()
    assert package.joinpath(
        "data", "datasets", "experimental-metal-eos-batch-47.csv"
    ).is_file()
    reference = BM2Reference(375.04, 254, 32.7)
    assert reference.pressure(375.04) == 32.7
    assert abs(reference.pressure(375.04 / 8) - 36608.7) < 1e-8
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
            or (
                record["scientific_validation"]["status"] == "not_reproduced"
                and record.get("catalog_access") == "explicit_selection"
                and record.get("default") is not True
                and record.get("default_for") is None
            )
        }
        if identifiers:
            expected[document["identifier"]] = identifiers

    # Check identities per material so a deferred record cannot replace an
    # eligible one, or cause eligible records in the same material to disappear.
    materials = list_materials()
    assert {
        material.identifier: {record.identifier for record in material.eos_records}
        for material in materials
    } == expected
    executable_ids = {identifier for ids in expected.values() for identifier in ids}
    assert {
        record.identifier
        for record in list_eos_records(include_reference_variants=True)
    } == executable_ids
    assert sum(len(material.eos_records) for material in materials) == len(
        executable_ids
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
