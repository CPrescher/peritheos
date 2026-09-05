import csv
import hashlib
import math
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document, validate_eosmat_document
from peritheos.eosmat import validate_pressure_calibration_references

ROOT = Path(__file__).parents[1]
MATERIAL_ID = "mg05fe05al05si05o3_bridgmanite"
DATASET_ID = "mg05fe05al05si05o3_bridgmanite_zhu_2020_zenodo_latticeparameters"
PREFERRED_ID = "mg05fe05al05si05o3_bridgmanite_zhu_2020_preferred_bm2_1"
POST_ANNEAL_ID = "mg05fe05al05si05o3_bridgmanite_zhu_2020_post_anneal_bm2_2"
COLD_ID = "mg05fe05al05si05o3_bridgmanite_zhu_2020_cold_compressed_bm2_3"
KOEMETS_ID = "mg05fe05al05si05o3_bridgmanite_koemets_2023_bm2_4"
ZHU_RECORD_IDS = {PREFERRED_ID, POST_ANNEAL_ID, COLD_ID}
DOI = "10.1029/2020JB019964"
DATASET_SHA256 = "ed88e17fb3028d1abe2c5fa57bd4bc38b7c61d6fe11b74c0ee952436bcaba3ce"


def _document_records_dataset_rows():
    document = get_material_document(MATERIAL_ID)
    validate_eosmat_document(document)
    records = {record["identifier"]: record for record in document["eos_records"]}
    dataset = next(
        dataset
        for dataset in document["datasets"]
        if dataset["identifier"] == DATASET_ID
    )
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    payload = path.read_bytes()
    assert hashlib.sha256(payload).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    return document, records, dataset, rows


def _bm2_pressure(parameters, volume):
    v0, k0 = parameters
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5)


def test_fa50_material_identity_volume_basis_and_source_scope():
    document, records, _, _ = _document_records_dataset_rows()

    assert document["formula"] == "Mg0.5Fe0.5Al0.5Si0.5O3"
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    assert document["units"]["volume"] == "angstrom^3/conventional_unit_cell"
    assert "20-atom conventional-cell basis" in document["cell_contents"]
    assert "lattice" not in document
    assert "atom_sites" not in document
    assert "LiNbO3-type at ambient" in document["notes"]
    assert "no independently fitted, fully low-spin EOS branch" in document["notes"]

    assert set(records) == ZHU_RECORD_IDS | {KOEMETS_ID}
    assert records[PREFERRED_ID]["default_for"] == "equilibrium"
    assert "default_for" not in records[POST_ANNEAL_ID]
    assert "default_for" not in records[COLD_ID]
    assert all(
        records[identifier]["reference"]["doi"] == DOI for identifier in ZHU_RECORD_IDS
    )
    assert all(
        record["scientific_validation"]["status"] == "primary_source_validated"
        for record in records.values()
    )


@pytest.mark.parametrize(
    "record_id,parameters,errors,pressure_range,volume_ratio",
    [
        (
            PREFERRED_ID,
            {"V0": 172.1, "K0": 229.0},
            {"V0": 0.4, "K0": 4.0},
            [26.8, 102.2],
            [0.7660564497, 0.9056263335],
        ),
        (
            POST_ANNEAL_ID,
            {"V0": 172.0, "K0": 230.0},
            {"V0": 0.7, "K0": 6.0},
            [29.0, 102.2],
            [0.7665018314, 0.9019469477],
        ),
        (
            COLD_ID,
            {"V0": 170.7, "K0": 240.0},
            {"V0": 0.3, "K0": 3.0},
            [27.5, 103.6],
            [0.7724757411, 0.9077366315],
        ),
    ],
)
def test_fa50_published_bm2_parameters_and_executable_envelopes(
    record_id, parameters, errors, pressure_range, volume_ratio
):
    document, records, _, _ = _document_records_dataset_rows()
    source = records[record_id]
    executable = Material.from_eosmat(
        document, record_identifiers=[record_id]
    ).get_eos_record(record_id)

    assert source["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": parameters,
    }
    assert source["parameter_errors"] == errors
    assert source["parameter_error_confidence"] is None
    assert source["fixed_parameters"] == []
    assert source["parameter_covariance"] is None
    assert source["temperature_ref"] == 300.0
    assert source["experimental_pressure_range_gpa"] == pressure_range
    assert source["validity"]["volume_ratio"] == volume_ratio
    assert source["fit_datasets"] == [DATASET_ID]
    assert (
        source["parameter_provenance"]["implicit_K0_prime"].endswith("K0'=4.")
        or "K0'=4" in source["parameter_provenance"]["implicit_K0_prime"]
    )

    assert executable.pressure(parameters["V0"]) == pytest.approx(0.0, abs=1e-12)
    midpoint = sum(pressure_range) / 2.0
    volume = executable.volume(midpoint, check_validity=True)
    assert executable.pressure(volume, check_validity=True) == pytest.approx(
        midpoint, rel=1e-11
    )
    with pytest.raises(ValueError, match="outside the published calibration"):
        executable.volume(pressure_range[0] - 0.1, check_validity=True)


def test_fa50_official_dataset_is_complete_and_selections_are_exact():
    _, _, dataset, rows = _document_records_dataset_rows()

    assert len(rows) == 44
    assert Counter(row["sample_number"] for row in rows) == {"1": 23, "2": 21}
    assert sum(int(row["used_in_preferred_sample1_fit"]) for row in rows) == 23
    assert sum(int(row["used_in_post_anneal_subset_fit"]) for row in rows) == 11
    assert sum(int(row["used_in_cold_compressed_sample2_fit"]) for row in rows) == 21
    assert all(float(row["pressure_uncertainty_gpa"]) == 0.1 for row in rows)
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert dataset["used_by_eos_records"] == [PREFERRED_ID, POST_ANNEAL_ID, COLD_ID]

    post_anneal_pressures = [
        float(row["pressure_gpa"])
        for row in rows
        if row["used_in_post_anneal_subset_fit"] == "1"
    ]
    assert post_anneal_pressures == [
        29.0,
        35.2,
        43.8,
        53.2,
        62.2,
        68.0,
        74.8,
        82.2,
        91.8,
        95.8,
        102.2,
    ]
    assert all(
        row["measured_after_current_step_laser_anneal"]
        == row["used_in_post_anneal_subset_fit"]
        for row in rows
        if row["sample_number"] == "1"
    )
    assert all(
        row["thermal_history"] == "cold_compressed"
        for row in rows
        if row["sample_number"] == "2"
    )
    assert all(
        abs(
            float(row["lattice_a_angstrom"])
            * float(row["lattice_b_angstrom"])
            * float(row["lattice_c_angstrom"])
            - float(row["volume_a3_conventional_cell"])
        )
        < 0.05
        for row in rows
    )


@pytest.mark.parametrize(
    "record_id,flag,expected_count,expected_refit,expected_refit_rmse,expected_published_rmse",
    [
        (
            PREFERRED_ID,
            "used_in_preferred_sample1_fit",
            23,
            [172.0737222855, 229.4797393590],
            0.9545433388,
            0.9577399457,
        ),
        (
            POST_ANNEAL_ID,
            "used_in_post_anneal_subset_fit",
            11,
            [172.0370217144, 230.0026422447],
            0.9935784693,
            0.9990515601,
        ),
        (
            COLD_ID,
            "used_in_cold_compressed_sample2_fit",
            21,
            [170.7468795339, 240.0511433144],
            0.6707208695,
            0.6882434168,
        ),
    ],
)
def test_fa50_unweighted_pressure_refits_reproduce_published_coefficients(
    record_id,
    flag,
    expected_count,
    expected_refit,
    expected_refit_rmse,
    expected_published_rmse,
):
    _, records, _, rows = _document_records_dataset_rows()
    source = records[record_id]
    fit_rows = [row for row in rows if row[flag] == "1"]
    assert len(fit_rows) == expected_count
    pressures = np.array([float(row["pressure_gpa"]) for row in fit_rows])
    volumes = np.array([float(row["volume_a3_conventional_cell"]) for row in fit_rows])

    result = least_squares(
        lambda candidate: _bm2_pressure(candidate, volumes) - pressures,
        x0=[source["eos"]["parameters"]["V0"], source["eos"]["parameters"]["K0"]],
        xtol=1e-14,
        ftol=1e-14,
        gtol=1e-14,
        max_nfev=10000,
    )
    refit_rmse = math.sqrt(
        float(np.mean((_bm2_pressure(result.x, volumes) - pressures) ** 2))
    )
    published = [source["eos"]["parameters"]["V0"], source["eos"]["parameters"]["K0"]]
    published_rmse = math.sqrt(
        float(np.mean((_bm2_pressure(published, volumes) - pressures) ** 2))
    )
    stored = source["scientific_validation"]["independent_refit"]

    assert result.x == pytest.approx(expected_refit, abs=2e-6)
    assert result.x == pytest.approx(
        [stored["refit_parameters"]["V0"], stored["refit_parameters"]["K0"]],
        abs=2e-6,
    )
    assert refit_rmse == pytest.approx(expected_refit_rmse, abs=5e-10)
    assert published_rmse == pytest.approx(expected_published_rmse, abs=5e-10)
    assert stored["result"] == "parity"


def test_fa50_spin_and_thermal_history_are_not_overinterpreted_as_branches():
    _, records, _, _ = _document_records_dataset_rows()
    preferred = records[PREFERRED_ID]
    post_anneal = records[POST_ANNEAL_ID]
    cold = records[COLD_ID]

    assert "mixed-spin population" in preferred["iron_speciation"]
    assert "not a fully low-spin branch" in preferred["iron_speciation"]
    assert "separate low-spin branch" in post_anneal["iron_speciation"]
    assert "metastable thermal-history comparison" in cold["iron_speciation"]
    assert "high spin on the A site" in cold["iron_speciation"]
    assert all("low_spin" not in record["identifier"] for record in records.values())

    unresolved = preferred["scientific_validation"]["unresolved_issues"]
    assert any(
        "DFT workbook" in issue and "no executable DFT EOS" in issue
        for issue in unresolved
    )
    assert any("two analyzed synthesis compositions" in issue for issue in unresolved)
    assert len(records) == 4


def test_fa50_pressure_scale_is_linked_but_rowwise_reduction_is_unavailable():
    _, records, _, _ = _document_records_dataset_rows()

    for record_id in ZHU_RECORD_IDS:
        source = records[record_id]
        calibration = source["pressure_calibration"]
        assert calibration["status"] == "partially_resolved"
        assert calibration["methods"][0]["material"] == "Au"
        assert calibration["methods"][0]["reference"]["doi"] == (
            "10.1073/pnas.0609013104"
        )
        assert calibration["methods"][0]["reference_eos_record"] == (
            "gold_fei_2007_vinet_2"
        )
        assert calibration["recalculation"]["status"] == (
            "missing_calibrant_observations"
        )
    validate_pressure_calibration_references()


def test_fa50_official_workbook_checksums_are_pinned():
    document, _, dataset, _ = _document_records_dataset_rows()
    files = {
        entry["name"]: entry["sha256"]
        for entry in document["source"]["official_data_deposit"]["files"]
    }
    assert files == {
        "latticeparameters.xlsx": "d1a94effe4d202284d166957679616c65e1e4a9d3d89906ccd424acc9ac73cee",
        "EoS-Ffplot-bulkv.xlsx": "aee87799ba994cce1dc455f77b8e65dca3b142164e0bb1fe6b12e8854efabe55",
        "EPMA.xlsx": "266091fdd085fdd144da51aae85e13bbef00f7c8ace7cc8ad98c632f34f03f4d",
        "DFTEoS.xlsx": "2deb83e0b1e9ea8c94f0c02f7f602e2588a50cb5e0104de34c5885864b716268",
    }
    assert (
        dataset["provenance"]["official_lattice_workbook_sha256"]
        == files["latticeparameters.xlsx"]
    )
    assert (
        dataset["provenance"]["official_eos_workbook_sha256"]
        == files["EoS-Ffplot-bulkv.xlsx"]
    )
