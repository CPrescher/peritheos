"""Compatibility and consumer-facing contracts for material family browsing."""

from collections import Counter
from dataclasses import replace

import jsonschema
import pytest

from peritheos import (
    EosmatError,
    Material,
    MaterialError,
    MaterialFamily,
    MaterialLookupError,
    eosmat_schema,
    get_material,
    get_material_document,
    get_material_family,
    group_materials,
    list_material_families,
    list_materials,
    search_eos_records,
    search_materials,
    validate_eosmat_document,
)
from peritheos.catalog_families import _parse_family_registry
from peritheos.materials import list_materials as compatibility_list_materials


def test_families_have_separate_lookup_namespace_and_curated_scope():
    assert [f.identifier for f in list_material_families()] == [
        "bridgmanite",
        "garnet",
        "magnesite",
        "mg_fe_monoxide",
        "naalsio4_calcium_ferrite",
        "olivine",
        "phase_d",
        "post_perovskite",
    ]
    composition_families = {
        "phase_d": ["phase_d", "al_bearing_phase_d", "fe_bearing_al_phase_d"],
        "naalsio4_calcium_ferrite": [
            "naalsio4_calcium_ferrite",
            "na093al102si100o4_calcium_ferrite",
            "na088al099fe013si094o4_calcium_ferrite",
        ],
        "magnesite": ["magnesite", "mg0991fe0008mn0001co3_magnesite"],
        "olivine": ["forsterite", "fayalite"],
        "garnet": ["pyrope", "almandine", "majorite"],
    }
    for family_id, member_ids in composition_families.items():
        assert {m.identifier for m in list_materials(family_id=family_id)} == set(
            member_ids
        )
    assert isinstance(get_material_family("bridgmanite"), MaterialFamily)
    assert get_material("bridgmanite").formula == "MgSiO3"
    assert get_material_family("bridgmanite").formula != "MgSiO3"
    for identifier in (
        "bridgmanite",
        "klb1_mg_perovskite",
        "al2o3_perovskite",
        "perovskite_orthorhombic",
        "fe088sio3_bridgmanite",
    ):
        assert get_material(identifier).family_id == "bridgmanite"
    for identifier in (
        "mgo",
        "feo",
        "mgfe94o_rhombohedral",
        "mgo_b2",
        "feo_b8_2",
        "klb1_ferropericlase",
    ):
        assert get_material(identifier).family_id == "mg_fe_monoxide"
    for identifier in (
        "gold",
        "mgo_liquid",
        "fe2o3",
        "ca_perovskite",
        "ca_perovskite_tetragonal",
        "casio3_perovskite_tetragonal",
        "klb1_ca_perovskite",
        "wadsleyite",
        "ringwoodite",
        "mgal2o4_cafe2o4",
        "namg2al5sio12_cf",
        "namg2al5sio12_nal",
        "phase_h",
        "cairo3_perovskite",
        "mggeo3_post_perovskite",
        "fesio3_post_perovskite_ii",
    ):
        assert get_material(identifier).family_id is None
    assert get_material("mgsio3_post_perovskite").family_id == "post_perovskite"


def test_grouped_catalog_is_a_lossless_deterministic_partition():
    materials = list_materials()
    groups = group_materials(materials)
    assert groups == group_materials(reversed(materials))
    assert Counter(
        m.identifier for group in groups for m in group.materials
    ) == Counter(m.identifier for m in materials)
    assert len(groups) < len(materials)
    for group in groups:
        assert tuple(m.identifier for m in group.materials) == tuple(
            sorted(m.identifier for m in group.materials)
        )
        if group.family is None:
            assert len(group.materials) == 1
        else:
            assert all(m.family_id == group.family.identifier for m in group.materials)
        assert all(m is get_material(m.identifier) for m in group.materials)
    assert group_materials([]) == ()
    with pytest.raises(MaterialError, match="Duplicate material identifier"):
        group_materials([materials[0], materials[0]])


def test_family_filters_combine_with_existing_material_and_record_filters():
    members = list_materials(family_id="bridgmanite")
    assert members == tuple(m for m in list_materials() if m.family_id == "bridgmanite")
    assert members == compatibility_list_materials(family_id="bridgmanite")
    assert list_materials(family_id="bridgmanite", formula="MgSiO3") == (
        get_material("bridgmanite"),
    )
    assert search_materials(family_id="bridgmanite", formula="MgSiO3") == (
        get_material("bridgmanite"),
    )
    assert search_materials(family_id="bridgmanite", formula="Au") == ()
    expected = {r.identifier for m in members for r in m.eos_records}
    assert {
        r.identifier for r in search_eos_records(family_id="bridgmanite")
    } == expected
    filtered = search_eos_records(family_id="bridgmanite", thermal=True, author="Tange")
    assert filtered
    assert filtered == tuple(
        r
        for r in search_eos_records(thermal=True, author="Tange")
        if r.identifier in expected
    )
    with pytest.raises(MaterialLookupError, match="Unknown material family"):
        get_material_family("missing_family")
    for operation in (list_materials, search_materials, search_eos_records):
        with pytest.raises(MaterialLookupError, match="Unknown material family"):
            operation(family_id="missing_family")
        with pytest.raises(MaterialError, match="Family identifier"):
            operation(family_id=42)


def test_grouping_filtered_results_does_not_expand_the_family():
    matches = search_materials(family_id="bridgmanite", formula="MgSiO3")
    groups = group_materials(iter(matches))
    assert len(groups) == 1
    assert groups[0].family == get_material_family("bridgmanite")
    assert groups[0].materials == matches
    # Material and family names may collide without merging standalone entries.
    standalone = replace(get_material("gold"), identifier="bridgmanite")
    groups = group_materials([standalone, get_material("mg09fe01sio3_bridgmanite")])
    assert len(groups) == 2
    assert sum(group.family is None for group in groups) == 1


def test_family_round_trip_old_documents_and_external_family_fallback():
    document = get_material_document("gold")
    assert "family_id" not in document
    material = Material.from_eosmat(document)
    assert material.family_id is None
    assert "family_id" not in material.to_eosmat()
    document["family_id"] = "external_family"
    validate_eosmat_document(document)
    jsonschema.validate(document, eosmat_schema())
    imported = Material.from_eosmat(document)
    assert imported.family_id == "external_family"
    assert imported.to_eosmat()["family_id"] == "external_family"
    assert Material.from_eosmat(imported.to_eosmat()).family_id == "external_family"
    groups = group_materials([imported])
    assert groups[0].family is None
    assert groups[0].materials == (imported,)
    assert "family_id" not in replace(imported, family_id=None).to_eosmat()
    assert imported.default_record().pressure(
        60.0
    ) == material.default_record().pressure(60.0)
    bundled = get_material("bridgmanite")
    assert Material.from_eosmat(bundled.to_eosmat()).family_id == "bridgmanite"


@pytest.mark.parametrize(
    "value", [None, 3, [], "", "Upper", "two words", "a__b", "a_", "é"]
)
def test_invalid_family_metadata_rejected_by_schema_and_python(value):
    document = get_material_document("gold")
    document["family_id"] = value
    with pytest.raises(EosmatError, match="family_id"):
        validate_eosmat_document(document)
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(document, eosmat_schema())
    if value is not None:
        with pytest.raises(MaterialError, match="family_id"):
            replace(get_material("gold"), family_id=value)


@pytest.mark.parametrize(
    "document",
    [
        {},
        [1],
        [{}],
        [{"identifier": "bad id", "name": "X", "description": "Scope"}],
        [{"identifier": "x", "name": "", "description": "Scope"}],
        [{"identifier": "x", "name": "X", "description": 3}],
        [{"identifier": "x", "name": "X", "description": "Scope", "formula": ""}],
        [{"identifier": "x", "name": "X", "description": "Scope"}] * 2,
    ],
)
def test_registry_rejects_malformed_or_duplicate_definitions(document):
    with pytest.raises(MaterialError):
        _parse_family_registry(document)


def test_bundled_catalog_rejects_unknown_membership_even_on_structure_only_card(
    monkeypatch,
):
    import peritheos.catalog as catalog

    document = get_material_document("gold")
    document.update(family_id="missing_family", eos_records=[])
    monkeypatch.setattr(catalog, "list_material_documents", lambda: ("gold",))
    monkeypatch.setattr(catalog, "get_material_document", lambda _: document)
    with pytest.raises(MaterialError, match="references unknown family"):
        catalog._catalog_index.__wrapped__()
