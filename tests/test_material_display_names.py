"""Family-first labels retain exact source compositions and searchable old names."""

import json
from pathlib import Path

import pytest

from peritheos import (
    Material,
    get_material_document,
    get_material_family,
    group_materials,
    list_material_documents,
    list_materials,
    search_eos_records,
    search_materials,
)

CASES = json.loads(
    (Path(__file__).parent / "fixtures/material-display-names.json").read_text(
        encoding="utf-8"
    )
)
MIXED_OXIDES = {
    case["identifier"] for case in CASES if case["family_id"] == "mg_fe_monoxide"
}


def test_ferropericlase_folder_label_uses_the_existing_family_identifier():
    family = get_material_family("mg_fe_monoxide")
    assert family.name == "Ferropericlase"
    groups = group_materials(list_materials(family_id=family.identifier))
    assert len(groups) == 1
    assert groups[0].family == family
    assert MIXED_OXIDES <= {material.identifier for material in groups[0].materials}


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["identifier"])
def test_family_first_names_preserve_composition_and_previous_searches(case):
    identifier = case["identifier"]
    for version in (3, 4):
        document = get_material_document(identifier, format_version=version)
        assert document["name"] == case["name"]
        assert document["formula"] == case["formula"]
        assert document.get("phase") == case["phase"]
        assert set(case["previous_aliases"]) <= set(document.get("aliases", []))
        if case["previous_name"] != case["name"]:
            assert case["previous_name"] in document["aliases"]

    # Consumers such as Dioptas append the formula to the name themselves.
    assert case["formula"] not in document["name"]
    label = f"{document['name']} ({document['formula']})"
    assert label.count(case["formula"]) == 1
    assert label.startswith(
        "Bridgmanite" if case["family_id"] == "bridgmanite" else "Ferropericlase"
    )
    assert identifier in {
        material.identifier for material in search_materials(case["previous_name"])
    }

    material = Material.from_eosmat(document, require_primary_validation=False)
    exported = material.to_eosmat()
    assert exported["name"] == case["name"]
    assert exported["formula"] == case["formula"]
    assert exported["aliases"] == document.get("aliases", [])


def test_name_audit_covers_every_family_member_and_keeps_endmembers_distinct():
    documents = [get_material_document(key) for key in list_material_documents()]
    bridgmanites = {
        d["identifier"] for d in documents if d.get("family_id") == "bridgmanite"
    }
    assert bridgmanites == {
        case["identifier"] for case in CASES if case["family_id"] == "bridgmanite"
    }
    assert MIXED_OXIDES == {
        d["identifier"]
        for d in documents
        if d.get("family_id") == "mg_fe_monoxide"
        and "Mg" in d["formula"]
        and "Fe" in d["formula"]
    }
    for identifier, name in (
        ("mgo", "Periclase B1"),
        ("mgo_b2", "B2 MgO"),
        ("feo", "Wüstite B1"),
        ("feo_b8_2", "FeO (B8)"),
        ("fe093o_rhombohedral", "Fe0.93O (rhombohedrally distorted B1)"),
    ):
        assert get_material_document(identifier)["name"] == name

    # Sorting consumer labels puts each family's entries in a contiguous block.
    names = sorted(d["name"] for d in documents)
    for prefix, count in (
        ("Bridgmanite", len(bridgmanites)),
        ("Ferropericlase", len(MIXED_OXIDES)),
    ):
        positions = [i for i, name in enumerate(names) if name.startswith(prefix)]
        assert len(positions) == count
        assert positions == list(range(positions[0], positions[0] + count))


@pytest.mark.parametrize(
    "query", ["ferropericlase", "magnesiowustite", "magnesiowüstite"]
)
def test_mixed_oxide_synonyms_find_all_compositions_and_their_eos_records(query):
    for keyword in ({"text": query}, {"alias": query}):
        assert MIXED_OXIDES <= {
            material.identifier for material in search_materials(**keyword)
        }
        records = search_eos_records(**keyword)
        matched_ids = {record.identifier for record in records}
        for identifier in MIXED_OXIDES:
            expected = {
                record.identifier
                for record in search_eos_records(
                    formula=get_material_document(identifier)["formula"],
                )
            }
            assert expected
            assert expected <= matched_ids
