import copy

import numpy as np
import pytest

import peritheos.eosmat as eosmat_module
import peritheos.pressure_calibrations as calibration_module
from peritheos import (
    get_material_document,
    get_pressure_calibration,
    get_pressure_calibration_document,
    list_material_documents,
    list_pressure_calibrations,
    list_xrd_pressure_standards,
    recalculate_ruby_pressure,
    recalculate_ruby_to_xrd_pressure,
    validate_pressure_calibration_references,
    xrd_standard_pressure,
)
from peritheos.errors import EosmatError, MaterialLookupError, ValidationError
from peritheos.materials import Material

EXPECTED_RUBY_SCALES = {
    "ruby_mao_1978",
    "ruby_mao_1986",
    "ruby_dewaele_2004",
    "ruby_holzapfel_2005",
    "ruby_dorogokupets_oganov_2007",
}


def test_bundled_ruby_calibration_documents_are_executable():
    assert set(list_pressure_calibrations()) == EXPECTED_RUBY_SCALES
    for identifier in list_pressure_calibrations():
        document = get_pressure_calibration_document(identifier)
        calibration = get_pressure_calibration(identifier)
        assert document["kind"] == "ruby_fluorescence"
        assert calibration.identifier == identifier
        assert calibration.pressure_from_ratio(1.0) == pytest.approx(0.0)


@pytest.mark.parametrize("identifier", sorted(EXPECTED_RUBY_SCALES))
def test_ruby_scales_round_trip_pressure_ratio_and_wavelength(identifier):
    calibration = get_pressure_calibration(identifier)
    pressure = np.array([0.0, 1.0, 50.0, 150.0, 300.0])
    ratio = calibration.wavelength_ratio(pressure)
    assert calibration.pressure_from_ratio(ratio) == pytest.approx(pressure)

    wavelength = calibration.wavelength_from_pressure(pressure)
    assert calibration.pressure_from_wavelength(wavelength) == pytest.approx(pressure)


def test_pressure_recalculation_between_ruby_scales_is_reversible():
    source_pressure = np.array([0.0, 10.0, 50.0, 100.0])
    recalculated = recalculate_ruby_pressure(
        source_pressure,
        "ruby_mao_1986",
        "ruby_dorogokupets_oganov_2007",
    )
    assert recalculated[-1] > source_pressure[-1]
    assert recalculate_ruby_pressure(
        recalculated,
        "ruby_dorogokupets_oganov_2007",
        "ruby_mao_1986",
    ) == pytest.approx(source_pressure)


def test_ruby_shift_api_and_input_validation():
    calibration = get_pressure_calibration("ruby_mao_1986")
    shift = calibration.wavelength_from_pressure(50.0) - calibration.reference_wavelength_nm
    assert calibration.pressure_from_shift(shift) == pytest.approx(50.0)
    with pytest.raises(ValidationError):
        calibration.pressure_from_ratio(0.99)
    with pytest.raises(ValidationError):
        calibration.wavelength_ratio(-1.0)
    with pytest.raises(ValidationError):
        calibration.wavelength_ratio(np.inf)
    with pytest.raises(ValidationError):
        calibration.pressure_from_wavelength(700.0, reference_wavelength_nm=0.0)
    with pytest.raises(ValidationError):
        calibration.pressure_from_shift(1.0, reference_wavelength_nm=np.nan)
    with pytest.raises(ValidationError):
        calibration.wavelength_from_pressure(1.0, reference_wavelength_nm=-1.0)
    with pytest.raises(ValidationError):
        get_pressure_calibration("ruby_holzapfel_2005").wavelength_ratio(2_000.0)


def test_every_identified_ruby_use_links_to_an_executable_calibration():
    ruby_methods = []
    for material_identifier in list_material_documents():
        for record in get_material_document(material_identifier)["eos_records"]:
            for method in record["pressure_calibration"]["methods"]:
                if method["kind"] == "ruby_fluorescence":
                    ruby_methods.append(method)

    assert len(ruby_methods) == 33
    assert all(
        method["reference_calibration_record"] in EXPECTED_RUBY_SCALES
        for method in ruby_methods
    )
    validate_pressure_calibration_references()


def test_holzapfel_2005_reference_is_not_the_unrelated_aip_doi():
    document = get_pressure_calibration_document("ruby_holzapfel_2005")
    assert document["parameters"] == {
        "A_gpa": 1845.0,
        "B": 14.7,
        "C": 7.5,
        "reference_wavelength_nm": 694.24,
    }
    assert document["reference"]["doi"] == "10.1080/09511920500147501"


def test_xrd_standard_pressure_uses_bundled_reference_eos():
    assert "gold_dewaele_2004_vinet_5" in list_xrd_pressure_standards()
    gold = Material.from_eosmat(
        get_material_document("gold"),
        record_identifiers=["gold_dewaele_2004_vinet_5"],
    ).get_eos_record("gold_dewaele_2004_vinet_5")
    volumes = gold.volume(np.array([0.0, 50.0, 100.0]))
    assert xrd_standard_pressure(
        "gold_dewaele_2004_vinet_5", volumes
    ) == pytest.approx([0.0, 50.0, 100.0])


def test_ruby_pressures_can_be_replaced_by_paired_gold_observations():
    gold = Material.from_eosmat(
        get_material_document("gold"),
        record_identifiers=["gold_dewaele_2004_vinet_5"],
    ).get_eos_record("gold_dewaele_2004_vinet_5")
    source_pressure = np.array([10.0, 50.0, 90.0])
    gold_pressure = np.array([10.5, 52.0, 95.0])
    result = recalculate_ruby_to_xrd_pressure(
        source_pressure,
        "ruby_mao_1986",
        "gold_dewaele_2004_vinet_5",
        gold.volume(gold_pressure),
    )
    assert result.target_pressure_gpa == pytest.approx(gold_pressure)
    assert result.pressure_difference_gpa == pytest.approx(
        gold_pressure - source_pressure
    )
    assert result.source_wavelength_ratio == pytest.approx(
        get_pressure_calibration("ruby_mao_1986").wavelength_ratio(source_pressure)
    )


def test_calibration_lookup_rejects_unknown_duplicate_and_non_ruby(monkeypatch):
    with pytest.raises(MaterialLookupError):
        get_pressure_calibration_document("unknown")

    mao = get_pressure_calibration_document("ruby_mao_1986")
    monkeypatch.setattr(
        calibration_module,
        "pressure_calibration_library",
        lambda: {"calibrations": [mao, mao]},
    )
    with pytest.raises(ValidationError, match="Duplicate"):
        get_pressure_calibration_document("ruby_mao_1986")

    non_ruby = dict(mao, kind="other_optical_gauge")
    monkeypatch.setattr(
        calibration_module,
        "get_pressure_calibration_document",
        lambda _identifier: non_ruby,
    )
    with pytest.raises(ValidationError, match="not a ruby scale"):
        get_pressure_calibration("not_ruby")


def test_xrd_recalculation_rejects_missing_standard_and_shape_mismatch():
    with pytest.raises(MaterialLookupError):
        xrd_standard_pressure("missing_eos", 10.0)
    with pytest.raises(ValidationError, match="broadcast-compatible"):
        recalculate_ruby_to_xrd_pressure(
            [10.0, 20.0],
            "ruby_mao_1986",
            "gold_dewaele_2004_vinet_5",
            [60.0, 58.0, 56.0],
        )


def test_eosmat_rejects_invalid_calibration_links_and_missing_targets(monkeypatch):
    document = get_material_document("gold")
    method = document["eos_records"][0]["pressure_calibration"]["methods"][0]
    method["reference_calibration_record"] = ""
    with pytest.raises(EosmatError, match="must be a non-empty string"):
        eosmat_module.validate_eosmat_document(document)

    method["reference_calibration_record"] = "ruby_mao_1986"
    method["kind"] = "other"
    with pytest.raises(EosmatError, match="requires a ruby_fluorescence method"):
        eosmat_module.validate_eosmat_document(document)

    monkeypatch.setattr(eosmat_module, "list_pressure_calibrations", lambda: ())
    with pytest.raises(EosmatError, match="references missing pressure calibration"):
        eosmat_module.validate_pressure_calibration_references()


def test_eosmat_rejects_other_malformed_pressure_calibration_fields():
    base = get_material_document("gold")

    documents = []
    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["methods"][0]["kind"] = (
        "invalid"
    )
    documents.append(document)

    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["methods"][0][
        "source_location"
    ] = ""
    documents.append(document)

    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["methods"][0][
        "reference"
    ] = 3
    documents.append(document)

    document = copy.deepcopy(base)
    method = document["eos_records"][0]["pressure_calibration"]["methods"][0]
    method["reference_eos_record"] = ""
    documents.append(document)

    document = copy.deepcopy(base)
    method = document["eos_records"][0]["pressure_calibration"]["methods"][0]
    method["kind"] = "equation_of_state"
    method.pop("reference_calibration_record")
    method.pop("reference")
    documents.append(document)

    document = copy.deepcopy(base)
    document["eos_records"][0]["pressure_calibration"]["recalculation"][
        "status"
    ] = "invalid"
    documents.append(document)

    for document in documents:
        with pytest.raises(EosmatError):
            eosmat_module.validate_eosmat_document(document)
