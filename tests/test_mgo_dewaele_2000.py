import csv
from importlib import resources

import pytest

from peritheos import get_material_document
from peritheos.materials import Material

RECORD_ID = "mgo_dewaele_2000_bm3_mgd_5"
DATASET_ID = "mgo_dewaele_2000_table2_pvt"


def _source_record():
    document = get_material_document("mgo")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    return document, record


def test_dewaele_2000_has_one_preferred_bm3_mgd_and_four_cold_sensitivities():
    document, record = _source_record()
    source_records = [
        item
        for item in document["eos_records"]
        if item["reference"].get("doi", "").lower() == "10.1029/1999jb900364"
    ]

    assert len(source_records) == 5
    assert source_records[0] == record
    assert [item["identifier"] for item in source_records[1:]] == [
        "mgo_dewaele_2000_bm2_sensitivity_1",
        "mgo_dewaele_2000_murnaghan_sensitivity_2",
        "mgo_dewaele_2000_vinet_sensitivity_3",
        "mgo_dewaele_2000_natural_strain3_sensitivity_4",
    ]
    assert record["reference"]["authors"][0] == "Dewaele"
    assert record["eos"]["type"] == "BM3"
    assert record["eos"]["parameters"] == {
        "V0": 74.71,
        "K0": 161.0,
        "K0_prime": 3.94,
    }
    assert record["thermal"]["type"] == "MieGruneisenDebye"
    assert record["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert record["thermal"]["fixed_parameters"] == [
        "Tr",
        "theta0",
        "gamma0",
        "n",
    ]
    assert record["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 800.0,
        "gamma0": 1.45,
        "q": 0.8,
        "n": 2.0,
    }


def test_dewaele_2000_reproduces_published_bm3_extrapolation():
    document, record = _source_record()
    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    eos_record = loaded.eos_records[0]
    volume = record["eos"]["parameters"]["V0"] * 0.667

    # The discussion following Table 3 reports 145 GPa for BM3 at V/V0=0.667.
    assert eos_record.pressure(volume, 300.0, check_validity=False) == pytest.approx(
        145.0, abs=0.1
    )


def test_dewaele_2000_table3_model_sensitivities_are_executable():
    document, _ = _source_record()
    identifiers = [
        "mgo_dewaele_2000_bm2_sensitivity_1",
        "mgo_dewaele_2000_murnaghan_sensitivity_2",
        "mgo_dewaele_2000_vinet_sensitivity_3",
        "mgo_dewaele_2000_natural_strain3_sensitivity_4",
    ]
    loaded = Material.from_eosmat(document, record_identifiers=identifiers)
    by_id = {record.identifier: record for record in loaded.eos_records}
    volume = 74.71 * 0.667

    assert by_id[identifiers[0]].pressure(volume) == pytest.approx(145.17149164)
    assert by_id[identifiers[1]].pressure(volume) == pytest.approx(153.09675297)
    assert by_id[identifiers[2]].pressure(volume) == pytest.approx(141.29063259)
    assert by_id[identifiers[3]].pressure(volume) == pytest.approx(138.91909057)


def test_dewaele_2000_table2_resource_contains_all_printed_rows():
    document, record = _source_record()
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert dataset["used_by_eos_records"] == [
        record["identifier"],
        "mgo_dewaele_2000_bm2_sensitivity_1",
        "mgo_dewaele_2000_murnaghan_sensitivity_2",
        "mgo_dewaele_2000_vinet_sensitivity_3",
        "mgo_dewaele_2000_natural_strain3_sensitivity_4",
    ]
    assert len(rows) == 61
    assert sum(row["series"] == "heated" for row in rows) == 41
    assert min(float(row["pressure_gpa"]) for row in rows) == 0.0
    assert max(float(row["pressure_gpa"]) for row in rows) == 53.0
    assert max(float(row["temperature_k"]) for row in rows) == 2474.0
