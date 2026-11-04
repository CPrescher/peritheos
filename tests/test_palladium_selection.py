"""An unresolved Pd fit cannot become a default, including by fallback."""

import pytest

from peritheos import Material, get_eos_record, get_material, get_material_document
from peritheos.errors import MaterialError


def test_palladium_default_is_on_hold_after_round_trip():
    material = get_material("palladium")
    assert material.default_record() is None
    assert material.default_equilibrium_record() is None
    assert not any(r.is_default for r in material.eos_records)

    record = get_eos_record("palladium_frost_2023_vinet_1")
    # Baty supplementary Table S3 independently tabulates Frost's Vinet
    # pressure at 14.0 A^3/atom; its two-decimal pressure sets the tolerance.
    assert record.pressure(4.0 * 14.0, 300.0) == pytest.approx(10.04, abs=0.005)
    assert record.volume(record.pressure(52.0, 300.0), 300.0) == pytest.approx(52.0)

    for document in (get_material_document("palladium"), material.to_eosmat()):
        reloaded = Material.from_eosmat(document)
        assert reloaded.default_record() is None
        assert reloaded.eosmat_metadata["recommendation_hold"]["scope"] == (
            "equilibrium_default"
        )
        assert not any(r.is_default for r in reloaded.eos_records)


@pytest.mark.parametrize("suffix", ["vinet_1", "bm3_2"])
def test_frost_failed_reproduction_is_explicitly_not_recommended(suffix):
    identifier = f"palladium_frost_2023_{suffix}"
    doc = get_material_document("palladium")
    raw = next(r for r in doc["eos_records"] if r["identifier"] == identifier)
    record = get_eos_record(identifier)
    assert record.scientific_validation_status == "not_reproduced"
    assert not record.is_default
    assert raw["catalog_access"] == "explicit_selection"
    assert raw["scientific_validation"]["source_equation_status"] == "reproduced"
    assert raw["scientific_validation"]["usage_recommendation"] == (
        "not_recommended_for_quantitative_use"
    )
    assert identifier not in {
        r.identifier for r in Material.from_eosmat(doc).eos_records
    }
    with pytest.raises(MaterialError, match="not_reproduced"):
        Material.from_eosmat(doc, record_identifiers=[identifier])
    accepted = Material.from_eosmat(
        doc, record_identifiers=[identifier], require_primary_validation=False
    )
    assert accepted.default_record() is None
    assert accepted.eos_records[0].eosmat_metadata["eos"] == raw["eos"]


def test_baty_experimental_eos_remains_explicitly_not_recommended():
    identifier = "palladium_baty_2024_bm3_1"
    document = get_material_document("palladium")
    raw = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    record = get_eos_record(identifier)

    assert record.scientific_validation_status == "not_reproduced"
    assert not record.is_default
    assert raw["catalog_access"] == "explicit_selection"
    assert raw["scientific_validation"]["usage_recommendation"] == (
        "not_recommended_for_quantitative_use"
    )
    assert identifier not in {
        r.identifier for r in Material.from_eosmat(document).eos_records
    }


def test_baty_dft_is_retained_as_a_separate_nondefault_zero_kelvin_record():
    record = get_eos_record("palladium_baty_2024_bm3_dft_2")

    assert record in get_material("palladium").eos_records
    assert not record.is_default
    assert record.reference_temperature == 0.0
    assert record.reference_volume == 58.4668
    assert record.pressure(record.reference_volume, 0.0) == pytest.approx(
        0.0, abs=1e-10
    )
