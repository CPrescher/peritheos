import pytest

from peritheos import get_material_document
from peritheos.materials import Material
from scripts.reproduce_noguchi_2013_casio3 import isothermal_bulk_modulus

RECORD_ID = "ca_perovskite_noguchi_2013_bm2_mgd_1"
TARGET_DOI = "10.1007/s00269-012-0549-1"


def _source_and_record():
    document = get_material_document("ca_perovskite")
    source = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    record = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]
    return document, source, record


def test_noguchi_2013_is_one_distinct_nondefault_700_k_record():
    document, source, record = _source_and_record()
    matching = [
        item
        for item in document["eos_records"]
        if item["reference"].get("doi", "").lower() == TARGET_DOI
    ]

    assert matching == [source]
    assert source["default"] is False
    assert source["record_kind"] == "published"
    assert source["equation_kind"] == "thermal"
    assert source["temperature_ref"] == 700.0
    assert record.reference_temperature == 700.0
    assert source["phase_selection"]["assigned_phase"] == ("cubic CaSiO3 perovskite")
    assert source["phase_selection"]["cell_basis"] == (
        "one-formula-unit pseudo-cubic cell (Z=1)"
    )
    assert source["fit_datasets"] == []

    default_record = next(
        item for item in document["eos_records"] if item.get("default", False)
    )
    assert default_record["identifier"] == "ca_perovskite_shim_2000_bm3_1"


def test_noguchi_2013_published_bm2_mgd_parameters_and_errors():
    _, source, record = _source_and_record()

    assert source["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 46.5, "K0": 207.0},
    }
    assert source["parameter_errors"] == {"V0": 0.1, "K0": 4.0}
    assert source["implicit_parameters"] == {"K0_prime": 4.0}
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert source["thermal"] == {
        "type": "MieGruneisenDebye",
        "model": "mie_gruneisen_debye",
        "debye_temperature_law": "integrated_gruneisen",
        "parameters": {
            "Tr": 700.0,
            "theta0": 1300.0,
            "gamma0": 2.7,
            "q": 1.2,
            "n": 5.0,
        },
        "parameter_errors": {
            "Tr": None,
            "theta0": 500.0,
            "gamma0": 0.3,
            "q": 0.8,
            "n": None,
        },
        "fixed_parameters": ["Tr", "n"],
    }
    assert record.reference_volume == 46.5
    assert record.pressure(46.5, 700.0) == pytest.approx(0.0, abs=1.0e-12)
    assert record.thermal_pressure_increment(46.5, 700.0) == pytest.approx(
        0.0, abs=1.0e-12
    )


def test_noguchi_2013_primary_scope_calibration_and_refit_parity():
    _, source, record = _source_and_record()

    assert source["experimental_pressure_range_gpa"] == [51.3, 127.2]
    assert source["experimental_temperature_range_k"] == [700.0, 2300.0]
    assert record.validity.pressure_gpa == (51.3, 127.2)
    assert record.validity.temperature_k == (700.0, 2300.0)
    calibration = source["pressure_calibration"]
    assert calibration["methods"][0]["material"] == "Pt"
    assert calibration["methods"][0]["reference"]["doi"] == (
        "10.1016/j.pepi.2003.09.018"
    )
    assert calibration["recalculation"]["status"] == ("reference_eos_not_bundled")
    validation = source["scientific_validation"]
    assert validation["primary_data_check"]["status"] == "parameterization_only"
    assert validation["independent_refit"]["result"] == "parity"
    assert len(validation["excluded_alternatives"]) == 3


def test_noguchi_2013_reproduces_published_300_k_extrapolation():
    _, _, record = _source_and_record()
    volume = record.volume(0.0, 300.0)
    modulus = isothermal_bulk_modulus(record, volume, 300.0)

    assert volume == pytest.approx(45.8, abs=0.02)
    # The paper's 225 GPa is calculated from coefficients printed at limited
    # precision; the rounded production coefficients give 226.47 GPa.
    assert modulus == pytest.approx(225.0, abs=1.6)
    assert record.pressure(34.0, 700.0) == pytest.approx(121.4445932889)
    assert record.pressure(34.0, 1600.0) == pytest.approx(129.1431629601)
    assert record.pressure(34.0, 2100.0) == pytest.approx(134.2503898744)

    with pytest.raises(ValueError, match="outside the published calibration/data"):
        record.pressure(46.5, 700.0, check_validity=True)
