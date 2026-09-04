import csv
import hashlib
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document, validate_eosmat_document
from peritheos.eosmat import validate_pressure_calibration_references

ROOT = Path(__file__).parents[1]
MATERIAL_ID = "mg09fe01sio3_bridgmanite"
RECORD_ID = "mg09fe01sio3_bridgmanite_mao_2015_low_spin_bm3_1"
DATASET_ID = "mg09fe01sio3_bridgmanite_mao_2015_table_s1_pv"
DOI = "10.1002/2015GL064400"
DATASET_SHA256 = "6778bbc5ce58ea671fd1af7332d237343545fb5cd238a85d6e36c41d5c15598b"


def _document_record_and_rows():
    document = get_material_document(MATERIAL_ID)
    validate_eosmat_document(document)
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_ID
    )
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
    executable = Material.from_eosmat(
        document, record_identifiers=[RECORD_ID]
    ).get_eos_record(RECORD_ID)
    return document, source, dataset, rows, executable


def _bm3_pressure(parameters, volume):
    v0, k0 = parameters
    eta = (v0 / volume) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5)


def test_mao_material_identity_structure_basis_and_source_scope():
    document, source, _, _, _ = _document_record_and_rows()

    assert document["formula"] == "Mg0.9Fe0.1SiO3"
    assert document["space_group"] == "Pbnm"
    assert document["space_group_number"] == 62
    assert document["formula_units_per_cell"] == 4
    lattice = document["lattice"]
    assert lattice["a"] * lattice["b"] * lattice["c"] == pytest.approx(164.2516902961)

    proxy_contents = Counter()
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        proxy_contents[site["element"]] += multiplicity * site["occupancy"]
    assert proxy_contents == {"Mg": 4.0, "Si": 4.0, "O": 12.0}
    assert "topology proxy" in document["cell_contents"]
    assert "do not encode" in document["cell_contents"]
    assert document["source"]["structure_reference"]["doi"] == ("10.2138/am.2012.4000")
    assert document["source"]["topology_proxy_reference"]["doi"] == (
        "10.1007/BF00308114"
    )

    assert len(document["eos_records"]) == 1
    assert source["reference"]["doi"] == DOI
    assert "low-spin" in source["label"]
    assert "B-site Fe3+ is low spin" in source["iron_speciation"]
    assert source["scientific_validation"]["status"] == ("primary_source_validated")
    unresolved = " ".join(source["scientific_validation"]["unresolved_issues"])
    assert "high-spin EOS" in unresolved
    assert "does not publish crystallographic site occupancies" in unresolved


def test_mao_published_bm3_parameters_reference_state_and_validity():
    _, source, _, _, record = _document_record_and_rows()

    assert source["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 162.6, "K0": 284.0, "K0_prime": 4.0},
    }
    assert source["parameter_errors"] == {
        "V0": 0.3,
        "K0": 4.0,
        "K0_prime": None,
    }
    assert source["parameter_error_confidence"] is None
    assert source["fixed_parameters"] == ["K0_prime"]
    assert source["parameter_covariance"] is None
    assert source["temperature_ref"] == 300.0
    assert source["experimental_pressure_range_gpa"] == [25.4, 125.6]
    assert source["validity"]["volume_ratio"] == [0.7691266913, 0.9257687577]
    assert source["fit_datasets"] == [DATASET_ID]

    assert record.pressure(162.6) == pytest.approx(0.0, abs=1.0e-12)
    assert record.eos.bulk_modulus(162.6) == pytest.approx(284.0)
    for pressure in (26.0, 80.7, 125.6):
        volume = record.volume(pressure, check_validity=True)
        assert record.pressure(volume, check_validity=True) == pytest.approx(
            pressure, rel=1.0e-11
        )
    with pytest.raises(ValueError, match="outside the published calibration"):
        record.volume(25.0, check_validity=True)


def test_mao_table_s1_transcription_selection_and_checksum():
    _, _, dataset, rows, _ = _document_record_and_rows()

    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 97
    assert Counter(row["run_number"] for row in rows) == {
        "1": 46,
        "2": 21,
        "3": 30,
    }
    assert Counter(row["spin_regime"] for row in rows) == {
        "high_spin": 24,
        "mixed_spin_transition": 10,
        "low_spin": 63,
    }
    assert sum(int(row["used_in_published_low_spin_fit"]) for row in rows) == 63
    assert sum(int(row["laser_annealed"]) for row in rows) == 28
    assert all(
        (float(row["pressure_gpa"]) > 25.0)
        == bool(int(row["used_in_published_low_spin_fit"]))
        for row in rows
    )

    assert rows[0] == {
        "source_order": "1",
        "run_number": "1",
        "used_in_published_low_spin_fit": "0",
        "spin_regime": "high_spin",
        "pressure_gpa": "0.0001",
        "pressure_standard_deviation_gpa": "0",
        "volume_a3_conventional_cell": "163.27",
        "volume_standard_deviation_a3": "0.04",
        "laser_annealed": "0",
        "ambient_high_spin_anchor": "1",
    }
    duplicate_pressure = [
        row
        for row in rows
        if row["run_number"] == "3" and row["pressure_gpa"] == "86.6"
    ]
    assert len(duplicate_pressure) == 2
    assert {row["volume_a3_conventional_cell"] for row in duplicate_pressure} == {
        "132.57"
    }
    assert {row["laser_annealed"] for row in duplicate_pressure} == {"0", "1"}


def test_mao_published_curve_reproduces_low_spin_table_s1():
    _, source, _, rows, record = _document_record_and_rows()
    fit_rows = [row for row in rows if row["used_in_published_low_spin_fit"] == "1"]
    volumes = np.array([float(row["volume_a3_conventional_cell"]) for row in fit_rows])
    pressures = np.array([float(row["pressure_gpa"]) for row in fit_rows])
    calculated = np.asarray(record.pressure(volumes), dtype=float)
    assert math.sqrt(float(np.mean((calculated - pressures) ** 2))) == pytest.approx(
        0.4607032343, abs=5.0e-10
    )

    reproduction = source["scientific_validation"]["numerical_reproduction"]
    for benchmark in reproduction["states"]:
        assert record.pressure(benchmark["source_volume_a3"]) == pytest.approx(
            benchmark["calculated_pressure_gpa"], abs=5.0e-10
        )
        assert (
            benchmark["absolute_difference_gpa"]
            <= benchmark["source_pressure_standard_deviation_gpa"]
        )


def test_mao_diagnostic_effective_variance_refit_has_parameter_parity():
    _, source, _, rows, _ = _document_record_and_rows()
    fit_rows = [row for row in rows if row["used_in_published_low_spin_fit"] == "1"]
    pressures = np.array([float(row["pressure_gpa"]) for row in fit_rows])
    pressure_sigmas = np.array(
        [float(row["pressure_standard_deviation_gpa"]) for row in fit_rows]
    )
    volumes = np.array([float(row["volume_a3_conventional_cell"]) for row in fit_rows])
    volume_sigmas = np.array(
        [float(row["volume_standard_deviation_a3"]) for row in fit_rows]
    )

    def effective_sigmas(parameters):
        step = 1.0e-5
        derivative = (
            _bm3_pressure(parameters, volumes + step)
            - _bm3_pressure(parameters, volumes - step)
        ) / (2.0 * step)
        return np.sqrt(pressure_sigmas**2 + (derivative * volume_sigmas) ** 2)

    parameters = np.array([162.6, 284.0])
    for _ in range(10):
        sigmas = effective_sigmas(parameters)
        result = least_squares(
            lambda candidate: (_bm3_pressure(candidate, volumes) - pressures) / sigmas,
            x0=parameters,
            xtol=1.0e-14,
            ftol=1.0e-14,
            gtol=1.0e-14,
        )
        if np.max(np.abs(result.x - parameters)) < 1.0e-12:
            break
        parameters = result.x

    stored = source["scientific_validation"]["independent_refit"]
    assert result.x == pytest.approx(
        [stored["refit_parameters"]["V0"], stored["refit_parameters"]["K0"]],
        abs=5.0e-7,
    )
    assert result.x == pytest.approx([162.5849494078, 283.9449500278], abs=5.0e-7)
    assert abs(result.x[0] - 162.6) < 0.3
    assert abs(result.x[1] - 284.0) < 4.0


def test_mao_pressure_scale_recalculation_limit_is_explicit():
    _, source, _, _, _ = _document_record_and_rows()
    calibration = source["pressure_calibration"]
    assert calibration["status"] == "partially_resolved"
    assert calibration["methods"][0]["material"] == "Pt"
    assert calibration["methods"][0]["reference"]["doi"] == ("10.1073/pnas.0609013104")
    assert calibration["recalculation"]["status"] == ("missing_calibrant_observations")
    assert (
        "row-wise Pt lattice parameters or volumes"
        in calibration["recalculation"]["notes"]
    )
    validate_pressure_calibration_references()
