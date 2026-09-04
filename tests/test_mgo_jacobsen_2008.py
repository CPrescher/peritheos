import csv
import hashlib
from importlib import resources

import numpy as np
import pytest

from peritheos import get_eos_record_document, get_material_document
from peritheos.eos.rt import BM3
from peritheos.fitting import fit_rt_eos
from peritheos.materials import Material

DOI = "10.2138/am.2008.2988"
KCL_RECORD_ID = "mgo_jacobsen_2008_bm3_kcl_mao1978"
HELIUM_RECORD_ID = "mgo_jacobsen_2008_bm3_helium_mao1986"
KCL_DATASET_ID = "mgo_jacobsen_2008_table2_kcl_compression"
HELIUM_DATASET_ID = "mgo_jacobsen_2008_table1_helium_compression"


def _document_and_records():
    document = get_material_document("mgo")
    records = {
        record["identifier"]: record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == DOI
    }
    loaded = Material.from_eosmat(
        document,
        record_identifiers=[KCL_RECORD_ID, HELIUM_RECORD_ID],
    )
    return document, records, loaded


def _dataset_rows(document, identifier):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == identifier
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    return dataset, payload, rows


def test_jacobsen_2008_preserves_two_pressure_scale_specific_bm3_records():
    _, records, loaded = _document_and_records()

    assert set(records) == {KCL_RECORD_ID, HELIUM_RECORD_ID}
    assert records[KCL_RECORD_ID]["eos"] == {
        "type": "BM3",
        "parameters": {"V0": 74.698, "K0": 164.1, "K0_prime": 4.05},
        "model": "birch_murnaghan_3",
    }
    assert records[HELIUM_RECORD_ID]["eos"] == {
        "type": "BM3",
        "parameters": {"V0": 74.687, "K0": 159.6, "K0_prime": 3.74},
        "model": "birch_murnaghan_3",
    }
    assert records[KCL_RECORD_ID]["parameter_errors"] == {
        "V0": 0.007,
        "K0": 0.9,
        "K0_prime": 0.04,
    }
    assert records[HELIUM_RECORD_ID]["parameter_errors"] == {
        "V0": 0.006,
        "K0": 0.6,
        "K0_prime": 0.03,
    }
    assert all(record["fixed_parameters"] == [] for record in records.values())
    assert all(record["parameter_covariance"] is None for record in records.values())
    assert all(record["temperature_ref"] == 300.0 for record in records.values())
    assert len(loaded.eos_records) == 2
    assert get_eos_record_document(KCL_RECORD_ID)["identifier"] == KCL_RECORD_ID
    assert get_eos_record_document(HELIUM_RECORD_ID)["identifier"] == HELIUM_RECORD_ID

    kcl_method = records[KCL_RECORD_ID]["pressure_calibration"]["methods"][0]
    assert kcl_method["reference_calibration_record"] == "ruby_mao_1978"
    helium_methods = records[HELIUM_RECORD_ID]["pressure_calibration"]["methods"]
    assert helium_methods[0]["reference_calibration_record"] == "ruby_mao_1986"
    assert helium_methods[1]["kind"] == "other_optical_gauge"


def test_jacobsen_2008_tables_are_complete_and_checksummed():
    document, records, _ = _document_and_records()
    helium_dataset, helium_payload, helium_rows = _dataset_rows(
        document, HELIUM_DATASET_ID
    )
    kcl_dataset, kcl_payload, kcl_rows = _dataset_rows(document, KCL_DATASET_ID)

    assert (
        hashlib.sha256(helium_payload).hexdigest()
        == helium_dataset["resource"]["sha256"]
    )
    assert hashlib.sha256(kcl_payload).hexdigest() == kcl_dataset["resource"]["sha256"]
    assert len(helium_rows) == 52
    assert sum(row["run"] == "1" for row in helium_rows) == 19
    assert sum(row["run"] == "2" for row in helium_rows) == 7
    assert sum(row["run"] == "3" for row in helium_rows) == 25
    assert (
        sum(row["pressure_marker"] == "diamond_raman_sun_2005" for row in helium_rows)
        == 6
    )
    assert helium_rows[0]["unit_cell_volume_a3"] == "74.698"
    assert helium_rows[-1]["reported_pressure_gpa"] == "111.0"
    assert helium_rows[-1]["zha_2000_mgo_scale_pressure_gpa"] == "118.1"
    assert len(kcl_rows) == 26
    assert kcl_rows[0]["reported_pressure_gpa"] == "3.6"
    assert kcl_rows[-1] == {
        "reported_pressure_gpa": "86.6",
        "reported_pressure_standard_deviation_gpa": "0.1",
        "unit_cell_volume_a3": "55.792",
        "unit_cell_volume_standard_deviation_a3": "0.006",
    }
    assert set(helium_dataset["used_by_eos_records"]) == set(records)
    assert kcl_dataset["used_by_eos_records"] == [KCL_RECORD_ID]


def test_jacobsen_2008_curves_reproduce_high_pressure_table_rows_and_invert():
    _, _, loaded = _document_and_records()
    kcl = loaded.get_eos_record(KCL_RECORD_ID)
    helium = loaded.get_eos_record(HELIUM_RECORD_ID)

    assert kcl.pressure(55.792, 300.0) == pytest.approx(86.6, abs=0.08)
    assert helium.pressure(52.239, 300.0) == pytest.approx(111.0, abs=0.25)
    for record, pressure in ((kcl, 86.6), (helium, 111.0)):
        volume = record.volume(pressure, 300.0, check_validity=True)
        assert record.pressure(volume, 300.0, check_validity=True) == pytest.approx(
            pressure, rel=1.0e-11
        )


@pytest.mark.parametrize(
    ("dataset_identifier", "expected_k0", "expected_k0_prime"),
    [
        (KCL_DATASET_ID, 164.1035648585, 4.0523892945),
        (HELIUM_DATASET_ID, 159.3535994115, 3.7451730152),
    ],
)
def test_jacobsen_2008_independent_errors_in_variables_refit(
    dataset_identifier, expected_k0, expected_k0_prime
):
    document = get_material_document("mgo")
    _, _, rows = _dataset_rows(document, dataset_identifier)
    if dataset_identifier == HELIUM_DATASET_ID:
        rows = rows[1:]

    volumes = np.array([float(row["unit_cell_volume_a3"]) for row in rows])
    pressures = np.array([float(row["reported_pressure_gpa"]) for row in rows])
    volume_sigma = np.array(
        [float(row["unit_cell_volume_standard_deviation_a3"]) for row in rows]
    )
    pressure_sigma = np.array(
        [float(row["reported_pressure_standard_deviation_gpa"]) for row in rows]
    )
    fit = fit_rt_eos(
        BM3,
        volumes,
        pressures,
        {"K0": 160.0, "K0_prime": 4.0},
        fixed={"V0": 74.698},
        pressure_sigma=pressure_sigma,
        volume_sigma=volume_sigma,
        absolute_sigma=True,
    )

    assert fit.success
    assert fit.parameters["K0"] == pytest.approx(expected_k0, abs=2.0e-7)
    assert fit.parameters["K0_prime"] == pytest.approx(expected_k0_prime, abs=2.0e-9)

    source_id = (
        KCL_RECORD_ID if dataset_identifier == KCL_DATASET_ID else HELIUM_RECORD_ID
    )
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == source_id
    )
    assert (
        abs(fit.parameters["K0"] - source["eos"]["parameters"]["K0"])
        < source["parameter_errors"]["K0"]
    )
    assert (
        abs(fit.parameters["K0_prime"] - source["eos"]["parameters"]["K0_prime"])
        < source["parameter_errors"]["K0_prime"]
    )


def test_jacobsen_2008_resolves_abstract_results_v0_difference_explicitly():
    _, records, _ = _document_and_records()
    helium = records[HELIUM_RECORD_ID]
    inconsistency = helium["scientific_validation"]["reported_inconsistencies"]

    assert inconsistency[0]["field"] == "preferred all-data helium V0"
    assert inconsistency[0]["abstract"] == "74.697(6) A^3"
    assert inconsistency[0]["results"] == "74.687(6) A^3"
    assert inconsistency[1]["field"] == "fixed-K0 sensitivity-fit V0 unit"
    assert inconsistency[1]["results"] == "V0=74.695(6) GPa"
    assert helium["eos"]["parameters"]["V0"] == 74.687
    assert (
        helium["scientific_validation"]["reported_parameterizations"][1]["V0_a3"]
        == 74.697
    )
    assert (
        helium["scientific_validation"]["reported_parameterizations"][2]["V0_a3"]
        == 74.695
    )
