from collections import Counter

import numpy as np
import pytest

from peritheos import (
    EOSRecord,
    Material,
    get_eos_record,
    get_material,
    get_material_document,
    search_eos_records,
    search_materials,
)
from peritheos.catalog_compat import LEGACY_RECORD_AUDIT
from peritheos.materials import _EOS_RECORD_CATALOG
from peritheos.materials import (
    search_eos_records as materials_search_eos_records,
)
from peritheos.materials import (
    search_materials as materials_search_materials,
)


def identifiers(items):
    return tuple(item.identifier for item in items)


def test_canonical_catalog_objects_execute_directly_from_normal_api():
    material = get_material("aragonite")
    record = get_eos_record("aragonite_martinez_1996_bm2_2")

    assert isinstance(material, Material)
    assert isinstance(record, EOSRecord)
    assert record in material.eos_records
    assert np.isfinite(record.pressure(record.reference_volume, 873.0))
    assert materials_search_materials(formula="MgO") == search_materials(formula="MgO")
    assert materials_search_eos_records(formula="MgO") == search_eos_records(
        formula="MgO"
    )


def test_every_historical_record_has_an_explicit_audited_disposition():
    canonical = set(identifiers(search_eos_records()))

    assert set(LEGACY_RECORD_AUDIT) == set(_EOS_RECORD_CATALOG)
    assert Counter(entry.relationship for entry in LEGACY_RECORD_AUDIT.values()) == {
        "equivalent": 9,
        "different": 17,
        "absent": 11,
    }
    assert all(
        entry.canonical_identifier in canonical
        for entry in LEGACY_RECORD_AUDIT.values()
        if entry.canonical_identifier is not None
    )


def test_material_text_name_formula_phase_and_alias_discovery():
    assert identifiers(search_materials(name="lithium fluoride")) == ("lif_b1",)
    assert identifiers(search_materials(formula="MgO")) == ("mgo",)
    assert "zinc_oxide_wurtzite" in identifiers(search_materials(phase="wurtzite"))
    assert "mgo" in identifiers(search_materials(alias="periclase"))
    assert identifiers(search_materials("post aragonite Pmmn")) == (
        "calcium_carbonate_post_aragonite",
    )
    assert get_material("Periclase") is get_material("mgo")
    assert search_materials("") == tuple(
        sorted(search_materials(), key=lambda x: x.identifier)
    )


def test_record_model_reference_and_capability_discovery():
    vinet = search_eos_records(model_family="Vinet")
    dewaele = search_eos_records(doi="https://doi.org/10.1103/PhysRevB.78.104102")

    assert vinet
    assert all(isinstance(record, EOSRecord) for record in vinet)
    assert identifiers(dewaele) == (
        "nickel_dewaele_2008_vinet_1",
        "silver_dewaele_2008_vinet_1",
    )
    assert (
        search_eos_records(doi=" HTTP://DX.DOI.ORG/10.1103/PhysRevB.78.104102 ")
        == dewaele
    )
    assert "rhenium_anzellini_2014_vinet_1" in identifiers(
        search_eos_records(author="Anzellini")
    )
    assert "silver_dewaele_2008_vinet_1" in identifiers(
        search_eos_records(reference="Table II Ni and Ag rows")
    )
    assert search_eos_records(thermal=True)
    assert search_eos_records(caloric=True)
    assert all(record.is_thermal for record in search_eos_records(thermal=True))
    assert all(
        record.parameter_errors or record.parameter_covariance
        for record in search_eos_records(uncertainty=True)
    )


def test_search_filters_are_combined_on_one_record_per_material():
    assert identifiers(
        search_materials(
            formula="C",
            model_family="double debye",
            caloric=True,
            pressure_gpa=500.0,
            temperature_k=5000.0,
        )
    ) == ("diamond",)
    assert not search_materials(formula="C", model_family="double debye", thermal=False)
    # These authors occur in separate gold records; text cannot leak across
    # records merely because both belong to the same material.
    assert search_materials("Anderson Fei") == ()


def test_calibration_range_queries_are_closed_and_have_explicit_semantics():
    zinc = "zinc_oxide_wurtzite_hanna_2011_bm2_1"

    assert zinc in identifiers(
        search_eos_records(pressure_gpa=(8.0, 8.3), temperature_k=300.0)
    )
    assert zinc not in identifiers(
        search_eos_records(pressure_gpa=(8.0, 9.0), temperature_k=300.0)
    )
    assert zinc in identifiers(
        search_eos_records(
            pressure_gpa=(8.0, 9.0),
            temperature_k=300.0,
            range_semantics="overlaps",
        )
    )
    assert zinc not in identifiers(search_eos_records(temperature_k=301.0))


def test_missing_calibration_ranges_do_not_claim_unbounded_coverage():
    assert "gold_sokolova_2013_holzapfel_4" not in identifiers(
        search_eos_records(pressure_gpa=0.0)
    )
    assert "gold_fei_2007_vinet_2" not in identifiers(
        search_eos_records(temperature_k=2000.0)
    )
    assert search_eos_records(pressure_gpa=1.0e9) == ()


def test_isothermal_reference_temperature_is_searchable_without_rewriting_metadata():
    identifier = "alumina_dewaele_2013_vinet_1"
    record = get_eos_record(identifier)

    assert identifier in identifiers(search_eos_records(temperature_k=300.0))
    assert identifier not in identifiers(search_eos_records(temperature_k=301.0))
    assert np.isinf(record.validity.temperature_k[1])
    exported = get_material("alumina").to_eosmat()["eos_records"][0]
    assert "experimental_temperature_range_k" not in exported
    assert "temperature_k" not in exported.get("validity", {})


@pytest.mark.parametrize(
    ("options", "message"),
    [
        ({"pressure_gpa": (2.0, 1.0)}, "ordered"),
        ({"temperature_k": 0.0}, "greater than zero"),
        ({"range_semantics": "inside"}, "contains.*overlaps"),
        ({"validation_status": "reviewed"}, "validation status"),
        ({"thermal": 1}, "thermal must be true"),
        ({"pressure_gpa": object()}, "finite number or two-value range"),
    ],
)
def test_search_rejects_ambiguous_or_invalid_range_options(options, message):
    with pytest.raises(ValueError, match=message):
        search_eos_records(**options)


def test_validation_status_filter_and_ordering_are_deterministic():
    records = search_eos_records(validation_status=("primary_source_validated",))
    assert len(records) == 150
    assert identifiers(records) == tuple(sorted(identifiers(records)))
    assert search_eos_records(validation_status="deferred") == ()


def test_lookup_errors_offer_close_identifier_matches():
    with pytest.raises(KeyError, match="zirconium_alpha") as material_error:
        get_material("zirconum_alpha")
    with pytest.raises(KeyError, match="gold_fei_2007_vinet_2"):
        get_eos_record("gold_fei_2007_vinet")
    with pytest.raises(KeyError, match="zirconium_alpha") as document_error:
        get_material_document("zirconum_alpha")
    assert material_error.value.operation == "lookup_material"
    assert document_error.value.operation == "lookup_material_document"
    assert document_error.value.context == {"identifier": "zirconum_alpha"}


def test_document_aliases_survive_executable_round_trip():
    material = get_material("mgo")
    payload = material.to_eosmat()

    assert "Periclase" in material.aliases
    assert payload["aliases"] == list(material.aliases)
    assert Material.from_eosmat(payload).aliases == material.aliases

    payload["eos_records"][0]["aliases"] = ["historical-mgo-record"]
    loaded = Material.from_eosmat(payload)
    aliased_record = loaded.get_eos_record("historical-mgo-record")
    assert aliased_record.aliases == ("historical-mgo-record",)
    assert loaded.to_eosmat()["eos_records"][0]["aliases"] == ["historical-mgo-record"]
