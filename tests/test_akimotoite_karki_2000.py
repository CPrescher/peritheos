from math import sqrt

import pytest

from peritheos import get_material_document
from peritheos.materials import Material

IDENTIFIER = "akimotoite_karki_2000_bm3_4"


def _source_and_record():
    document = get_material_document("akimotoite")
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == IDENTIFIER
    )
    record = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).eos_records[0]
    return document, source, record


def test_karki_2000_parameters_static_reference_and_scope():
    document, source, record = _source_and_record()

    assert source["record_kind"] == "published"
    assert source["reference"]["doi"] == "10.2138/am-2000-2-309"
    assert source["eos"] == {
        "type": "BM3",
        "parameters": {
            "V0": 252.75,
            "K0": 224.0,
            "K0_prime": 4.18,
        },
        "model": "birch_murnaghan_3",
    }
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": None,
        "K0_prime": None,
    }
    assert source["fixed_parameters"] == []
    assert source["fit_datasets"] == []
    assert source["experimental_pressure_range_gpa"] == [0.0, 40.0]
    assert source["pressure_range_status"] == "theoretical"
    assert source["temperature_ref"] == 0.0
    assert document["formula_units_per_cell"] == 6
    assert record.reference_temperature == 0.0
    assert record.reference_volume == pytest.approx(252.75)
    assert (
        source["scientific_validation"]["primary_data_check"]["status"]
        == "theoretical_parameterization_only"
    )


def test_karki_2000_hexagonal_cell_volume_basis():
    _, _, record = _source_and_record()

    table_1_volume = sqrt(3.0) / 2.0 * 4.686**2 * 13.291
    assert table_1_volume == pytest.approx(252.75, abs=0.002)
    assert record.reference_volume == pytest.approx(table_1_volume, abs=0.002)


def test_karki_2000_published_bm3_curve_and_round_trip():
    _, _, record = _source_and_record()

    assert record.pressure(252.75, 0.0, check_validity=True) == pytest.approx(
        0.0, abs=1.0e-12
    )
    assert record.pressure(0.9 * 252.75, 0.0, check_validity=True) == pytest.approx(
        29.429006104947266, rel=1.0e-12
    )

    volume_at_40_gpa = record.volume(40.0, 0.0, check_validity=True)
    assert volume_at_40_gpa / 252.75 == pytest.approx(0.8740149072261841, rel=1.0e-10)
    assert record.pressure(volume_at_40_gpa, 0.0, check_validity=True) == pytest.approx(
        40.0, rel=1.0e-10
    )
