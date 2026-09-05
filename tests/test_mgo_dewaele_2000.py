import csv
import hashlib
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.fitting import fit_joint_eos
from peritheos.materials import Material

RECORD_ID = "mgo_dewaele_2000_bm3_mgd_5"
DATASET_ID = "mgo_dewaele_2000_table2_pvt"
FEI_DATASET_ID = "mgo_fei_1999_table1_1100k"


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
    assert record["fit_datasets"] == [DATASET_ID, FEI_DATASET_ID]


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


def test_fei_1999_resource_contains_complete_1100_k_isotherm():
    document, record = _source_record()
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == FEI_DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))

    assert hashlib.sha256(payload).hexdigest() == dataset["resource"]["sha256"]
    assert dataset["used_by_eos_records"] == [record["identifier"]]
    assert dataset["provenance"]["source_pdf_sha256"] == (
        "98f44f1a44c5b0bcfd18baa84a460e36d53d46f5512e67b636b085e7bf9c0ac3"
    )
    assert len(rows) == 16
    assert {float(row["temperature_k"]) for row in rows} == {1100.0}
    gold_pressure = [float(row["gold_pressure_gpa"]) for row in rows]
    molar_volume = [float(row["mgo_molar_volume_cm3_mol"]) for row in rows]
    assert [gold_pressure[0], gold_pressure[-1]] == [0.0, 25.51]
    assert [molar_volume[0], molar_volume[-1]] == [11.63, 10.19]
    assert rows[0]["gold_pressure_uncertainty_gpa"] == ""
    assert rows[0]["mgo_molar_volume_uncertainty_cm3_mol"] == ""


def test_dewaele_2000_complete_thermal_input_refit_reproduces_q():
    document, record = _source_record()
    datasets = {item["identifier"]: item for item in document["datasets"]}

    dewaele_path = resources.files("peritheos.data").joinpath(
        datasets[DATASET_ID]["resource"]["path"]
    )
    with dewaele_path.open(newline="", encoding="utf-8") as handle:
        dewaele_rows = [
            row for row in csv.DictReader(handle) if row["series"] == "heated"
        ]

    fei_path = resources.files("peritheos.data").joinpath(
        datasets[FEI_DATASET_ID]["resource"]["path"]
    )
    with fei_path.open(newline="", encoding="utf-8") as handle:
        fei_rows = list(csv.DictReader(handle))

    pressure = np.array(
        [float(row["pressure_gpa"]) for row in dewaele_rows]
        + [float(row["gold_pressure_gpa"]) for row in fei_rows]
    )
    temperature = np.array(
        [float(row["temperature_k"]) for row in dewaele_rows]
        + [float(row["temperature_k"]) for row in fei_rows]
    )
    conventional_cell_volume = np.array(
        [float(row["mgo_unit_cell_volume_a3"]) for row in dewaele_rows]
        + [
            float(row["mgo_molar_volume_cm3_mol"]) * 4.0 / 0.602214076
            for row in fei_rows
        ]
    )

    loaded = Material.from_eosmat(document, record_identifiers=[RECORD_ID])
    eos_record = loaded.eos_records[0]
    executable = eos_record.eos
    volume = conventional_cell_volume * eos_record.volume_scale
    result = fit_joint_eos(
        type(executable),
        type(executable.rt_eos),
        volume=volume,
        temperature=temperature,
        pressure=pressure,
        initial={"q": record["thermal"]["parameters"]["q"]},
        fixed={
            "rt_eos.V0": record["eos"]["parameters"]["V0"] * eos_record.volume_scale,
            "rt_eos.K0": record["eos"]["parameters"]["K0"],
            "rt_eos.K0_prime": record["eos"]["parameters"]["K0_prime"],
            "Tr": record["thermal"]["parameters"]["Tr"],
            "theta0": record["thermal"]["parameters"]["theta0"],
            "gamma0": record["thermal"]["parameters"]["gamma0"],
            "n": record["thermal"]["parameters"]["n"],
        },
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={"q": (1.0e-9, 10.0)},
        max_nfev=5000,
    )

    published_pressure = np.asarray(executable.pressure(volume, temperature))
    published_rmse = np.sqrt(np.mean((published_pressure - pressure) ** 2))
    refit_rmse = np.sqrt(np.mean(np.asarray(result.residuals) ** 2))

    assert pressure.size == 57
    assert result.parameters["q"] == pytest.approx(0.785888544, abs=1.0e-6)
    assert result.standard_errors["q"] == pytest.approx(0.124307071, abs=1.0e-6)
    assert (
        abs(result.parameters["q"] - 0.8) < record["thermal"]["parameter_errors"]["q"]
    )
    assert published_rmse == pytest.approx(1.030760335, abs=1.0e-6)
    assert refit_rmse == pytest.approx(1.030641245, abs=1.0e-6)
