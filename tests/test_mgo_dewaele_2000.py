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


def test_dewaele_2000_is_one_preferred_bm3_mgd_record():
    document, record = _source_record()
    source_records = [
        item
        for item in document["eos_records"]
        if item["reference"].get("doi", "").lower() == "10.1029/1999jb900364"
    ]

    assert source_records == [record]
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


def test_dewaele_2000_table2_resource_contains_all_printed_rows():
    document, record = _source_record()
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))

    assert dataset["used_by_eos_records"] == [record["identifier"]]
    assert len(rows) == 61
    assert sum(row["series"] == "heated" for row in rows) == 41
    assert min(float(row["pressure_gpa"]) for row in rows) == 0.0
    assert max(float(row["pressure_gpa"]) for row in rows) == 53.0
    assert max(float(row["temperature_k"]) for row in rows) == 2474.0
