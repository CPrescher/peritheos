"""Parameter origin stays distinct from paper, model and validation categories."""

import copy

import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    eosmat_schema,
    get_eos_record,
    get_material_document,
    list_material_documents,
    search_eos_records,
    search_materials,
    validate_eosmat_document,
)
from peritheos.errors import EosmatError, MaterialError

METHODS = ("experimental", "theoretical", "hybrid", "unknown")


def test_every_bundled_record_has_an_explicit_method_and_rationale():
    for identifier in list_material_documents():
        document = get_material_document(identifier)
        for record in document["eos_records"]:
            assert record["determination_method"] in METHODS, record["identifier"]
            assert record["parameter_provenance"]["determination_method"].strip()
        if not document["eos_records"]:
            continue
        # Include the two deferred source records in metadata interchange checks.
        material = Material.from_eosmat(document, require_primary_validation=False)
        exported = material.to_eosmat()
        by_id = {r["identifier"]: r for r in document["eos_records"]}
        for record in material.eos_records:
            assert (
                record.determination_method
                == by_id[record.identifier]["determination_method"]
            )
        for record in exported["eos_records"]:
            assert (
                record["determination_method"]
                == by_id[record["identifier"]]["determination_method"]
            )
            assert (
                record["parameter_provenance"]["determination_method"]
                == by_id[record["identifier"]]["parameter_provenance"][
                    "determination_method"
                ]
            )


@pytest.mark.parametrize("method", METHODS)
def test_methods_validate_and_round_trip(method):
    document = get_material_document("srsio3_cubic_perovskite")
    document["eos_records"][0]["determination_method"] = method
    Draft202012Validator(eosmat_schema()).validate(document)
    validate_eosmat_document(document)
    material = Material.from_eosmat(document)
    identifier = document["eos_records"][0]["identifier"]
    assert material.get_eos_record(identifier).determination_method == method
    result = next(
        r for r in material.to_eosmat()["eos_records"] if r["identifier"] == identifier
    )
    assert result["determination_method"] == method


@pytest.mark.parametrize(
    "invalid", ["Experimental", "computed", "", None, True, 1, [], {}]
)
def test_invalid_method_is_rejected_by_schema_and_loader(invalid):
    document = get_material_document("srsio3_cubic_perovskite")
    document["eos_records"][0]["determination_method"] = invalid
    errors = list(Draft202012Validator(eosmat_schema()).iter_errors(document))
    assert any(list(error.path)[-1] == "determination_method" for error in errors)
    with pytest.raises(EosmatError, match="determination_method"):
        validate_eosmat_document(document)
    with pytest.raises(EosmatError, match="determination_method"):
        Material.from_eosmat(document)


def test_missing_legacy_method_is_unknown_without_inference_or_mutation():
    document = get_material_document("srsio3_cubic_perovskite")
    for record in document["eos_records"]:
        record.pop("determination_method")
        record["parameter_provenance"].pop("determination_method")
    original = copy.deepcopy(document)
    Draft202012Validator(eosmat_schema()).validate(document)
    validate_eosmat_document(document)
    material = Material.from_eosmat(document)
    assert all(r.determination_method == "unknown" for r in material.eos_records)
    assert all(
        "determination_method" not in r for r in material.to_eosmat()["eos_records"]
    )
    assert document == original


@pytest.mark.parametrize(
    "identifier, expected",
    [
        ("srsio3_cubic_xiao_2013_experimental_bm2_1", "experimental"),
        ("srsio3_cubic_xiao_2013_gga_bm2_2", "theoretical"),
        ("palladium_baty_2024_bm3_1", "experimental"),
        ("palladium_baty_2024_bm3_dft_2", "theoretical"),
        ("diamond_correa_2008_double_debye_log_moment_5", "theoretical"),
        ("diamond_correa_2008_dewaele_anchored", "hybrid"),
        ("mgo_oganov_2003_paw_small_core_static_bm3", "theoretical"),
        ("mgo_oganov_2003_pressure_corrected_298k_bm3", "hybrid"),
        ("ca_perovskite_ono_2013_bm3_log_thermal", "hybrid"),
        ("kcl_b2_dewaele_2012_vinet_3", "hybrid"),
        ("iron_dewaele_2006_bm3", "experimental"),
        ("iron_dewaele_2006_vinet_thermal", "hybrid"),
        ("iron_zhang_2025_fit1_birch_murnaghan_3_mgd", "hybrid"),
        ("platinum_holmes_1989_vinet_1", "hybrid"),
        ("bridgmanite_hamahata_2000_md_300k_bm3", "hybrid"),
        ("mgsio3_post_perovskite_komabayashi_2008_thermal_bm3", "hybrid"),
        ("mgo_sokolova_2013_holzapfel_4", "experimental"),
        ("mgo_li_2006_bm3_absolute_acoustic", "experimental"),
        ("coesite_iv_bykova_2018_am05_static_bm3_refit", "theoretical"),
        ("neon_fcc_hemley_1989_bm3_refit", "experimental"),
        ("bridgmanite_holland_2013_mpv_modified_tait", "unknown"),
        # A theoretical pressure calibrant does not reclassify measured sample data.
        ("nacl_b2_sakai_2011_holmes_pt_bm3", "experimental"),
    ],
)
def test_source_specific_classification_boundaries(identifier, expected):
    assert get_eos_record(identifier).determination_method == expected


@pytest.mark.parametrize("method", METHODS)
def test_catalog_method_filters_partition_records_and_match_one_record(method):
    records = search_eos_records(determination_method=method)
    expected = tuple(
        r for r in search_eos_records() if r.determination_method == method
    )
    assert records == expected
    assert records
    assert all(
        any(r.determination_method == method for r in material.eos_records)
        for material in search_materials(determination_method=method)
    )
    # The experimental sibling must not satisfy the filter for a named GGA record.
    matches = search_materials(
        text="srsio3_cubic_xiao_2013_gga_bm2_2", determination_method=method
    )
    assert bool(matches) is (method == "theoretical")


@pytest.mark.parametrize("search", [search_eos_records, search_materials])
@pytest.mark.parametrize("invalid", ["Experimental", "computed", "", True, [], {}])
def test_invalid_method_filter_fails_even_without_matching_materials(search, invalid):
    with pytest.raises(MaterialError, match="determination method"):
        search(formula="no_such_material", determination_method=invalid)
