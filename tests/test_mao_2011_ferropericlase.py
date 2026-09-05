import csv
import hashlib
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.eosmat import validate_pressure_calibration_references
from scripts.reproduce_mao_2011_ferropericlase import (
    PUBLISHED,
    branch_mask,
    fit_static_branch,
    fit_thermal_branch,
    load_rows,
    published_curve_residuals,
    reconstructed_low_spin_fraction,
)

ROOT = Path(__file__).parents[1]
MATERIAL_ID = "mg075fe025o_ferropericlase"
HS_ID = "mg075fe025o_ferropericlase_mao_2011_high_spin_bm3_thermal"
LS_ID = "mg075fe025o_ferropericlase_mao_2011_low_spin_bm3_thermal"
DATASET_ID = "mg075fe025o_ferropericlase_mao_2011_table_s1_pvt"
DATASET_SHA256 = "45f6a2e1c8f275cab6837d6525c18df7b9c7bde7198235a37da23802efa4355c"
DOI = "10.1029/2011GL049915"


def _document_records_dataset_rows():
    document = get_material_document(MATERIAL_ID)
    validate_eosmat_document(document)
    records = {record["identifier"]: record for record in document["eos_records"]}
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    payload = path.read_bytes()
    rows = list(csv.DictReader(payload.decode("utf-8").splitlines()))
    executable = Material.from_eosmat(document)
    return document, records, dataset, path, payload, rows, executable


def test_mao_2011_material_identity_structure_and_pure_branch_scope():
    document, records, _, _, _, _, executable = _document_records_dataset_rows()

    assert document["formula"] == "Mg0.75Fe0.25O"
    assert document["phase"].startswith("B1 rocksalt ferropericlase")
    assert document["space_group"] == "Fm-3m"
    assert document["space_group_number"] == 225
    assert document["formula_units_per_cell"] == 4
    assert document["lattice"]["a"] ** 3 == pytest.approx(76.3383389845)

    contents = Counter()
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        contents[site["element"]] += multiplicity * site["occupancy"]
    assert contents == {"Mg": 3.0, "Fe": 1.0, "O": 4.0}

    assert set(records) == {HS_ID, LS_ID}
    for identifier in (HS_ID, LS_ID):
        source = records[identifier]
        assert source["reference"]["doi"] == DOI
        assert source["phase_scope"].endswith(
            "mixed-spin observations are outside this executable model."
        )
        assert source["scientific_validation"]["status"] == ("primary_source_validated")
        excluded = source["scientific_validation"]["excluded_parameterizations"]
        assert "mixed-spin crossover" in excluded[0]["scope"]
        assert "No surrogate crossover record" in excluded[0]["reason"]

    assert len(executable.eos_records) == 2
    for record in executable.eos_records:
        assert isinstance(record.eos, ThermalReferenceStateEOS)
        assert isinstance(record.eos.rt_eos, BM3)
        assert record.eos.reference_volume_law == "integrated_expansivity"
        assert record.eos.thermal_expansion_law == "constant"


def test_mao_2011_corrected_parameters_fixed_free_status_and_pressure_scale():
    _, records, _, _, _, _, executable = _document_records_dataset_rows()
    high = records[HS_ID]
    low = records[LS_ID]

    assert high["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 76.34, "K0": 162.0, "K0_prime": 4.0},
    }
    assert high["parameter_errors"] == {
        "V0": 0.01,
        "K0": 1.0,
        "K0_prime": None,
    }
    assert high["fixed_parameters"] == ["V0", "K0_prime"]
    assert low["eos"]["parameters"] == {
        "V0": 74.4,
        "K0": 166.0,
        "K0_prime": 4.0,
    }
    assert low["parameter_errors"] == {
        "V0": 0.6,
        "K0": 7.0,
        "K0_prime": None,
    }
    assert low["fixed_parameters"] == ["K0_prime"]
    assert low["scientific_validation"]["correction"] == {
        "doi": "10.1029/2011GL050814",
        "applied_changes": [
            "Use corrected low-spin V0=74.4 A^3 rather than 74.44 A^3.",
            "Use the corrected logistic expression for n_LS.",
            "Treat the static P-V fit as the 300 K fit.",
        ],
    }

    for source in (high, low):
        assert source["parameter_error_confidence"] is None
        assert source["parameter_covariance"] is None
        assert source["thermal"]["fixed_parameters"] == ["Tr"]
        assert source["thermal"]["parameter_errors"] == {
            "Tr": None,
            "alpha0": 5e-7,
            "dK_dT": 0.002,
        }
        calibration = source["pressure_calibration"]
        assert calibration["methods"][0]["reference_eos_record"] == (
            "gold_fei_2007_vinet_2"
        )
        assert calibration["methods"][0]["reference"]["doi"] == (
            "10.1073/pnas.0609013104"
        )
        assert calibration["recalculation"]["status"] == (
            "missing_calibrant_observations"
        )
    validate_pressure_calibration_references()

    by_id = {record.identifier: record for record in executable.eos_records}
    assert by_id[HS_ID].reference_volume == pytest.approx(76.34)
    assert by_id[LS_ID].reference_volume == pytest.approx(74.4)
    assert by_id[HS_ID].pressure(76.34, 300.0) == pytest.approx(0.0, abs=1e-12)
    assert by_id[LS_ID].pressure(74.4, 300.0) == pytest.approx(0.0, abs=1e-12)


def test_mao_2011_complete_table_checksum_and_reconstructed_ranges():
    _, records, dataset, _, payload, rows, _ = _document_records_dataset_rows()

    assert hashlib.sha256(payload).hexdigest() == DATASET_SHA256
    assert dataset["resource"]["sha256"] == DATASET_SHA256
    assert [column["name"] for column in dataset["columns"]] == list(rows[0])
    assert len(rows) == 195
    assert Counter(row["temperature_k"] for row in rows) == {
        "300": 58,
        "1200": 38,
        "1500": 35,
        "1800": 33,
        "2000": 31,
    }
    assert Counter(row["inferred_spin_regime"] for row in rows) == {
        "high_spin": 97,
        "mixed_spin_crossover": 54,
        "low_spin": 44,
    }
    assert sum(int(row["used_in_diagnostic_high_spin_refit"]) for row in rows) == 97
    assert sum(int(row["used_in_diagnostic_low_spin_refit"]) for row in rows) == 27
    assert not any(
        row["temperature_k"] == "2000" and row["inferred_spin_regime"] == "low_spin"
        for row in rows
    )
    assert not any(
        row["temperature_k"] == "2000"
        and row["used_in_diagnostic_low_spin_refit"] == "1"
        for row in rows
    )

    ranges = dataset["spin_classification"]["temperature_ranges_gpa"]
    assert ranges == {
        "300": {
            "high_spin_through": 54.6,
            "crossover_from": 55.2,
            "low_spin_from": 81.5,
        },
        "1200": {
            "high_spin_through": 68.5,
            "crossover_from": 73.7,
            "low_spin_from": 100.0,
        },
        "1500": {
            "high_spin_through": 78.3,
            "crossover_from": 81.9,
            "low_spin_from": 123.3,
        },
        "1800": {
            "high_spin_through": 82.7,
            "crossover_from": 85.2,
            "low_spin_from": 124.6,
        },
        "2000": {
            "high_spin_through": 85.4,
            "crossover_from": 86.8,
            "low_spin_from": None,
        },
    }
    assert records[HS_ID]["spin_state_scope"][
        "inferred_last_high_spin_pressure_gpa"
    ] == {key: value["high_spin_through"] for key, value in ranges.items()}
    assert records[LS_ID]["spin_state_scope"][
        "inferred_first_low_spin_pressure_gpa"
    ] == {key: value["low_spin_from"] for key, value in ranges.items()}


def test_mao_2011_spin_fractions_and_row_flags_are_reproducible():
    rows = load_rows()
    calculated = reconstructed_low_spin_fraction(rows)

    assert calculated == pytest.approx(rows["stored_low_spin_fraction"], abs=5.1e-9)
    assert (
        branch_mask(rows, "high_spin").tolist()
        == rows["high_spin_flag"].astype(bool).tolist()
    )
    assert (
        branch_mask(rows, "low_spin").tolist()
        == rows["low_spin_flag"].astype(bool).tolist()
    )
    assert calculated.min() == pytest.approx(-0.1019949307, abs=5e-10)
    assert calculated.max() == pytest.approx(1.0581007251, abs=5e-10)


def test_mao_2011_static_and_thermal_diagnostic_refits_have_parity():
    _, records, _, _, _, _, _ = _document_records_dataset_rows()
    rows = load_rows()

    high_static = fit_static_branch(rows, "high_spin")
    low_static = fit_static_branch(rows, "low_spin")
    assert high_static == pytest.approx([162.39174473], abs=5e-7)
    assert low_static == pytest.approx([74.50083421, 165.00280482], abs=5e-7)
    assert abs(high_static[0] - 162.0) < 1.0
    assert abs(low_static[0] - 74.4) < 0.6
    assert abs(low_static[1] - 166.0) < 7.0

    for branch, identifier in (("high_spin", HS_ID), ("low_spin", LS_ID)):
        result = fit_thermal_branch(rows, branch)
        source = records[identifier]
        stored = source["scientific_validation"]["independent_refit"]
        expected = np.array(
            [
                stored["refined_parameters"]["alpha0"],
                stored["refined_parameters"]["dK_dT"],
            ]
        )
        assert result.parameters == pytest.approx(expected, abs=5e-11)
        assert (
            result.observations == source["spin_state_scope"]["selected_observations"]
        )
        assert abs(result.parameters[0] - PUBLISHED[branch]["alpha0"]) < 5e-7
        assert abs(result.parameters[1] - PUBLISHED[branch]["dK_dT"]) < 0.002

        residuals = published_curve_residuals(rows, branch)
        assert np.sqrt(np.mean(residuals**2)) == pytest.approx(
            stored["published_curve_pressure_rmse_gpa"], abs=5e-10
        )
        assert np.max(np.abs(residuals)) == pytest.approx(
            stored["published_curve_max_absolute_pressure_residual_gpa"], abs=5e-10
        )
