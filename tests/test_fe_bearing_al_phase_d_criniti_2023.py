"""Primary-source checks for Criniti et al. (2023) Fe-bearing Al-phase D."""

import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document
from peritheos.eos.rt import BM3, Vinet

ROOT = Path(__file__).resolve().parents[1]
DATASET_SHA256 = "b87e61f6d441d6685e677443b2527be112b400b14637c176301a1b4edd0901f3"
BM3_ID = "fe_bearing_al_phase_d_criniti_2023_bm3_1"
VINET_ID = "fe_bearing_al_phase_d_criniti_2023_vinet_2"


def _document_and_rows():
    document = get_material_document("fe_bearing_al_phase_d")
    dataset = document["datasets"][0]
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, dataset, path, rows


def _effective_variance_refit(model_class, rows):
    selected = [row for row in rows if row["used_in_room_temperature_fit"] == "1"]
    pressure = np.array([float(row["pressure_gpa"]) for row in selected])
    pressure_sigma = np.array(
        [float(row["pressure_uncertainty_gpa"]) for row in selected]
    )
    volume = np.array([float(row["volume_a3"]) for row in selected])
    volume_sigma = np.array([float(row["volume_uncertainty_a3"]) for row in selected])
    parameters = np.array([83.68, 166.0, 4.5])

    for _ in range(100):
        model = model_class(V0=parameters[0], K0=parameters[1], K0_prime=parameters[2])
        dp_dv = -np.asarray(model.bulk_modulus(volume), dtype=float) / volume
        sigma_effective = np.sqrt(pressure_sigma**2 + (dp_dv * volume_sigma) ** 2)

        def residual(candidate):
            candidate_model = model_class(
                V0=candidate[0], K0=candidate[1], K0_prime=candidate[2]
            )
            return (
                np.asarray(candidate_model.pressure(volume), dtype=float) - pressure
            ) / sigma_effective

        optimization = least_squares(
            residual,
            parameters,
            bounds=([82.0, 100.0, 0.0], [85.0, 250.0, 10.0]),
            x_scale="jac",
        )
        if np.allclose(optimization.x, parameters, rtol=1.0e-13, atol=1.0e-14):
            return optimization.x
        parameters = optimization.x

    raise AssertionError("effective-variance refit did not converge")


def test_criniti_2023_material_identity_structure_and_volume_basis():
    document, _, _, _ = _document_and_rows()

    assert document["formula"] == "Al1.53Fe0.22Si0.86O6H3.33"
    assert document["space_group"] == "P63 22"
    assert document["space_group_number"] == 182
    assert document["formula_units_per_cell"] == 1
    assert document["source"]["structure_file_sha256"] == (
        "8448888b776c5c786f33b68a1735613921330a90bdf25aafe33a2f3fb0fa677e"
    )

    lattice = document["lattice"]
    cell_volume = (
        lattice["a"]
        * lattice["b"]
        * lattice["c"]
        * math.sin(math.radians(lattice["gamma"]))
    )
    assert cell_volume == pytest.approx(83.703, abs=0.001)

    cell_contents = {}
    for site in document["atom_sites"]:
        cell_contents[site["element"]] = cell_contents.get(site["element"], 0.0) + (
            site["multiplicity"] * site["occupancy"]
        )
    assert cell_contents == pytest.approx(
        {"Al": 1.506, "Si": 0.8468, "Fe": 0.2166, "O": 6.0}
    )
    assert "Hydrogen" in document["notes"]
    assert "P-31m refinement" in document["notes"]


def test_criniti_2023_table1_transcription_and_fit_selection():
    _, dataset, path, rows = _document_and_rows()

    assert hashlib.sha256(path.read_bytes()).hexdigest() == DATASET_SHA256
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 26
    assert sum(int(row["used_in_room_temperature_fit"]) for row in rows) == 21
    assert {int(row["run_number"]) for row in rows} == {1, 2}
    assert rows[0] == {
        "source_order": "1",
        "run_number": "1",
        "used_in_room_temperature_fit": "1",
        "pressure_gpa": "0.35",
        "pressure_uncertainty_gpa": "0.02",
        "lattice_a_angstrom": "4.7430",
        "lattice_a_uncertainty_angstrom": "0.0005",
        "lattice_c_angstrom": "4.2873",
        "lattice_c_uncertainty_angstrom": "0.0003",
        "volume_a3": "83.525",
        "volume_uncertainty_a3": "0.017",
    }
    assert rows[20]["pressure_gpa"] == "33.30"
    assert rows[20]["used_in_room_temperature_fit"] == "1"
    assert rows[21]["pressure_gpa"] == "40.75"
    assert rows[21]["used_in_room_temperature_fit"] == "0"
    assert rows[-1]["volume_a3"] == "67.91"
    assert dataset["used_by_eos_records"] == [BM3_ID, VINET_ID]


def test_criniti_2023_published_records_execute_and_reproduce_table1():
    document, _, _, _ = _document_and_rows()
    material = Material.from_eosmat(document)
    records = {record.identifier: record for record in material.eos_records}

    assert set(records) == {BM3_ID, VINET_ID}
    assert material.default_record().identifier == BM3_ID
    assert records[BM3_ID].reference_temperature == 293.0
    assert records[VINET_ID].reference_temperature == 293.0

    # Independent, non-reference-state checks against printed Table 1 rows.
    assert records[BM3_ID].pressure(80.98) == pytest.approx(5.84, abs=0.04)
    assert records[VINET_ID].pressure(80.98) == pytest.approx(5.84, abs=0.02)
    assert records[BM3_ID].pressure(71.302) == pytest.approx(37.92, abs=0.20)
    assert records[VINET_ID].pressure(71.302) == pytest.approx(37.92, abs=0.15)

    for record in records.values():
        for volume in (83.525, 75.01, 71.302):
            assert record.volume(record.pressure(volume)) == pytest.approx(volume)


def test_criniti_2023_effective_variance_refits_recover_published_coefficients():
    document, _, _, rows = _document_and_rows()
    by_identifier = {record["identifier"]: record for record in document["eos_records"]}

    bm3_refit = _effective_variance_refit(BM3, rows)
    vinet_refit = _effective_variance_refit(Vinet, rows)
    assert bm3_refit == pytest.approx(
        [83.68041718, 166.25494346, 4.45597358], abs=2.0e-6
    )
    assert vinet_refit == pytest.approx(
        [83.68347936, 165.47809134, 4.62083697], abs=2.0e-6
    )

    bm3_published = by_identifier[BM3_ID]["eos"]["parameters"]
    vinet_published = by_identifier[VINET_ID]["eos"]["parameters"]
    assert bm3_refit == pytest.approx(
        [bm3_published["V0"], bm3_published["K0"], bm3_published["K0_prime"]],
        abs=0.05,
    )
    assert vinet_refit == pytest.approx(
        [
            vinet_published["V0"],
            vinet_published["K0"],
            vinet_published["K0_prime"],
        ],
        abs=0.05,
    )


def test_criniti_2023_preserves_bm3_vinet_choice_and_source_limitations():
    document, _, _, _ = _document_and_rows()
    records = {record["identifier"]: record for record in document["eos_records"]}

    assert records[BM3_ID]["default_for"] == "equilibrium"
    assert "default_for" not in records[VINET_ID]
    assert records[BM3_ID]["fixed_parameters"] == []
    assert records[VINET_ID]["fixed_parameters"] == []
    assert records[BM3_ID]["parameter_covariance"] is None
    assert records[VINET_ID]["parameter_covariance"] is None
    assert records[BM3_ID]["validity"]["pressure_gpa"] == [0.0, 38.0]
    assert records[VINET_ID]["validity"]["pressure_gpa"] == [0.0, 38.0]
    assert records[BM3_ID]["validity"]["volume_ratio"][0] == pytest.approx(
        71.302 / 83.68
    )
    assert records[VINET_ID]["validity"]["volume_ratio"][0] == pytest.approx(
        71.302 / 83.68
    )
    assert records[BM3_ID]["pressure_calibration"]["recalculation"]["status"] == (
        "missing_calibrant_observations"
    )
    assert any(
        "Low-spin state EOS" in issue
        for issue in records[BM3_ID]["scientific_validation"]["unresolved_issues"]
    )
