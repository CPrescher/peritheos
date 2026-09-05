import csv
import hashlib
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.materials import Material
from peritheos.units import molar_volume_to_cell_volume

DOI = "10.2138/am-1999-0308"
HYDRO_RECORD_ID = "mgo_fei_1999_bm3_300k_ne"
HOT_RECORD_ID = "mgo_fei_1999_bm3_1100k_nacl"
NONHYDRO_RECORD_ID = "mgo_fei_1999_bm3_300k_nonhydrostatic"
HYDRO_DATASET_ID = "mgo_fei_1999_table1_300k_ne"
HOT_DATASET_ID = "mgo_fei_1999_table1_1100k"
NONHYDRO_DATASET_ID = "mgo_fei_1999_table1_300k_nonhydrostatic"


def _document():
    return get_material_document("mgo")


def _records(document):
    return {
        record["identifier"]: record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == DOI.lower()
    }


def _dataset(document, identifier):
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == identifier
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    return dataset, payload, list(csv.DictReader(payload.decode("utf-8").splitlines()))


def test_fei_1999_registers_three_distinct_published_bm3_isotherms():
    document = _document()
    records = _records(document)

    assert set(records) == {
        HYDRO_RECORD_ID,
        HOT_RECORD_ID,
        NONHYDRO_RECORD_ID,
    }
    assert records[HYDRO_RECORD_ID]["eos"]["parameters"] == {
        "V0": pytest.approx(74.79068),
        "K0": 160.0,
        "K0_prime": 4.15,
    }
    assert records[HOT_RECORD_ID]["eos"]["parameters"] == {
        "V0": pytest.approx(77.248277),
        "K0": 135.0,
        "K0_prime": 4.2,
    }
    assert records[NONHYDRO_RECORD_ID]["eos"]["parameters"] == {
        "V0": pytest.approx(74.79068),
        "K0": 185.0,
        "K0_prime": 4.0,
    }
    assert records[HYDRO_RECORD_ID]["eos"]["parameters"]["V0"] == pytest.approx(
        molar_volume_to_cell_volume(11.26, 4, from_unit="cm^3/mol"), abs=5.0e-7
    )
    assert records[HOT_RECORD_ID]["eos"]["parameters"]["V0"] == pytest.approx(
        molar_volume_to_cell_volume(11.63, 4, from_unit="cm^3/mol"), abs=5.0e-7
    )
    assert all(
        record["eos"]["model"] == "birch_murnaghan_3" for record in records.values()
    )
    assert all(
        record["fixed_parameters"] == ["V0", "K0_prime"] for record in records.values()
    )
    assert records[NONHYDRO_RECORD_ID]["scientific_validation"][
        "usage_recommendation"
    ] == ("not_recommended_for_quantitative_use")


def test_fei_1999_table1_datasets_are_complete_and_checksummed():
    document = _document()
    expected = {
        HYDRO_DATASET_ID: (
            19,
            "d596636a19359483fcc4778ce0440d0d617a2e500923546eb59d0ee5a50404c5",
        ),
        HOT_DATASET_ID: (
            16,
            "78821a17b70db4612ba8f6c3780b1a479bf99af7ef18a2ca9b4c2726d0bc0ca6",
        ),
        NONHYDRO_DATASET_ID: (
            27,
            "292aa2343bd31678abd9bf95c8fbb3668246633b525b2284afda2ee8bf51f009",
        ),
    }

    total_rows = 0
    for identifier, (row_count, digest) in expected.items():
        dataset, payload, rows = _dataset(document, identifier)
        assert hashlib.sha256(payload).hexdigest() == digest
        assert dataset["resource"]["sha256"] == digest
        assert dataset["provenance"]["source_pdf_sha256"] == (
            "98f44f1a44c5b0bcfd18baa84a460e36d53d46f5512e67b636b085e7bf9c0ac3"
        )
        assert len(rows) == row_count
        total_rows += len(rows)

    assert total_rows == 62

    _, _, hydro = _dataset(document, HYDRO_DATASET_ID)
    assert hydro[0]["nacl_pressure_gpa"] == "0"
    assert hydro[0]["mgo_molar_volume_cm3_mol"] == "11.26"
    assert hydro[-1]["nacl_pressure_gpa"] == "23.22"

    _, _, nonhydro = _dataset(document, NONHYDRO_DATASET_ID)
    assert nonhydro[15]["nacl_pressure_gpa"] == "27.62"
    assert nonhydro[16]["gold_pressure_gpa"] == "30.56"
    assert nonhydro[-1]["nacl_pressure_gpa"] == ""
    assert nonhydro[-1]["gold_pressure_gpa"] == "65.32"


def _bm3_shape(molar_volume, v0, k0_prime):
    eta = (v0 / molar_volume) ** (1.0 / 3.0)
    return 1.5 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))


@pytest.mark.parametrize(
    (
        "record_id",
        "dataset_id",
        "v0",
        "k0_prime",
        "published_k0",
        "expected_k0",
        "expected_error",
        "published_rmse",
        "refit_rmse",
    ),
    [
        (
            HYDRO_RECORD_ID,
            HYDRO_DATASET_ID,
            11.26,
            4.15,
            160.0,
            159.013715338,
            0.748283375,
            0.296511166,
            0.283161068,
        ),
        (
            HOT_RECORD_ID,
            HOT_DATASET_ID,
            11.63,
            4.2,
            135.0,
            136.581065899,
            0.780252062,
            0.364102633,
            0.322614115,
        ),
        (
            NONHYDRO_RECORD_ID,
            NONHYDRO_DATASET_ID,
            11.26,
            4.0,
            185.0,
            184.832371646,
            0.768696882,
            0.679922305,
            0.679301371,
        ),
    ],
)
def test_fei_1999_complete_observation_refits_have_parameter_parity(
    record_id,
    dataset_id,
    v0,
    k0_prime,
    published_k0,
    expected_k0,
    expected_error,
    published_rmse,
    refit_rmse,
):
    document = _document()
    record = _records(document)[record_id]
    _, _, rows = _dataset(document, dataset_id)
    volume = np.asarray([float(row["mgo_molar_volume_cm3_mol"]) for row in rows])

    if dataset_id == NONHYDRO_DATASET_ID:
        pressure = np.asarray(
            [float(row["nacl_pressure_gpa"]) for row in rows[:16]]
            + [float(row["gold_pressure_gpa"]) for row in rows[16:]]
        )
    else:
        pressure = np.asarray([float(row["nacl_pressure_gpa"]) for row in rows])

    shape = _bm3_shape(volume, v0, k0_prime)
    fitted_k0 = float(shape @ pressure / (shape @ shape))
    residuals = fitted_k0 * shape - pressure
    standard_error = float(
        np.sqrt(np.sum(residuals**2) / (len(rows) - 1) / (shape @ shape))
    )
    calculated_published_rmse = float(
        np.sqrt(np.mean((published_k0 * shape - pressure) ** 2))
    )
    calculated_refit_rmse = float(np.sqrt(np.mean(residuals**2)))

    stored = record["scientific_validation"]
    assert fitted_k0 == pytest.approx(expected_k0, abs=1.0e-9)
    assert standard_error == pytest.approx(expected_error, abs=1.0e-9)
    assert calculated_published_rmse == pytest.approx(published_rmse, abs=1.0e-9)
    assert calculated_refit_rmse == pytest.approx(refit_rmse, abs=1.0e-9)
    assert abs(fitted_k0 - published_k0) < record["parameter_errors"]["K0"]
    assert stored["independent_refit"]["parameters"]["K0"] == pytest.approx(fitted_k0)
    assert stored["independent_refit"]["standard_errors"]["K0"] == pytest.approx(
        standard_error
    )
    assert stored["numerical_reproduction"][
        "published_curve_rmse_gpa"
    ] == pytest.approx(calculated_published_rmse)
    assert stored["numerical_reproduction"]["refit_curve_rmse_gpa"] == pytest.approx(
        calculated_refit_rmse
    )


def test_fei_1999_records_execute_on_their_published_isotherms():
    document = _document()
    records = _records(document)

    for identifier, temperature in (
        (HYDRO_RECORD_ID, 300.0),
        (HOT_RECORD_ID, 1100.0),
        (NONHYDRO_RECORD_ID, 300.0),
    ):
        loaded = Material.from_eosmat(document, record_identifiers=[identifier])
        executable = loaded.eos_records[0]
        v0 = records[identifier]["eos"]["parameters"]["V0"] * executable.volume_scale
        assert executable.pressure(v0, temperature) == pytest.approx(0.0, abs=1.0e-12)
        volume = executable.volume(10.0, temperature, check_validity=True)
        assert executable.pressure(
            volume, temperature, check_validity=True
        ) == pytest.approx(10.0, abs=1.0e-10)


def test_fei_1999_preserves_pressure_scale_and_temperature_scope():
    records = _records(_document())

    assert records[HYDRO_RECORD_ID]["pressure_calibration"]["methods"][0]["reference"][
        "authors"
    ] == ["Birch"]
    assert len(records[HOT_RECORD_ID]["pressure_calibration"]["methods"]) == 2
    assert len(records[NONHYDRO_RECORD_ID]["pressure_calibration"]["methods"]) == 2
    assert records[HOT_RECORD_ID]["temperature_ref"] == 1100.0
    assert records[HOT_RECORD_ID]["validity"]["temperature_k"] == [1100.0, 1100.0]
    assert (
        "single high-temperature reference isotherm"
        in records[HOT_RECORD_ID]["validity"]["notes"][1]
    )
