import csv
import hashlib
import json
import subprocess
import sys
from collections import Counter
from importlib import resources
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document

DATASET_ID = "iron_zhang_2025_tables_s1_s3_s4_pvt"
RECORD_IDS = {
    "iron_zhang_2025_fit1_birch_murnaghan_3_mgd",
    "iron_zhang_2025_fit2_birch_murnaghan_3_mgd",
    "iron_zhang_2025_fit5_vinet_mgd",
}
EXPECTED_SHA256 = "763a371196e3b4660d43cfaa924d1a96cf236bccef937cb1e9970b9eadb3c614"


def _dataset_path() -> Path:
    return Path(
        resources.files("peritheos").joinpath(
            "data", "datasets", "iron-zhang-2025-tables-s1-s3-s4-pvt.csv"
        )
    )


def _rows():
    with _dataset_path().open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_zhang_source_selection_and_provenance_are_explicit():
    rows = _rows()
    assert len(rows) == 1313
    assert Counter(row["source_table"] for row in rows) == {
        "Table S1": 1078,
        "Table S3": 132,
        "Table S4": 103,
    }
    assert Counter(row["data_class"] for row in rows) == {
        "static_experiment": 1078,
        "dynamic_experiment": 132,
        "ab_initio": 103,
    }
    assert Counter(row["temperature_provenance"] for row in rows) == {
        "experimental_reported": 1078,
        "theory_assigned_zhuang_2021_hugoniot": 62,
        "theory_assigned_zhuang_2021_isentrope": 60,
        "model_calculated_huang_2022": 10,
        "ab_initio_calculated": 103,
    }
    static = [row for row in rows if row["source_table"] == "Table S1"]
    assert Counter(row["pressure_correction_applied"] for row in static) == {
        "1": 641,
        "0": 437,
    }
    assert Counter(row["pressure_scale"] for row in static) == {
        "Fei et al. (2007)": 499,
        "Ye et al. (2018)": 142,
        "": 437,
    }
    anomalous = [row for row in rows if row["source_cell_note"]]
    assert len(anomalous) == 1
    assert anomalous[0]["source_row"] == "461"
    assert anomalous[0]["molar_volume_cm3_mol"] == "5.4983362575646"
    assert all(row["used_in_fit"] == "1" for row in rows)

    dynamic = [row for row in rows if row["source_table"] == "Table S3"]
    zhuang_hugoniot = [
        row
        for row in dynamic
        if row["dynamic_path"] == "Hugoniot"
        and not row["source_reference_as_published"].startswith("Huang")
    ]
    zhuang_ramp = [row for row in dynamic if row["dynamic_path"] == "Ramp"]
    hugoniot_pressure = np.asarray(
        [float(row["pressure_gpa_fit"]) for row in zhuang_hugoniot]
    )
    ramp_pressure = np.asarray([float(row["pressure_gpa_fit"]) for row in zhuang_ramp])
    assert np.asarray([float(row["temperature_k"]) for row in zhuang_hugoniot]) == (
        pytest.approx(
            0.01491 * hugoniot_pressure**2 + 21.77176 * hugoniot_pressure - 395.42764,
            abs=3e-9,
        )
    )
    assert np.asarray([float(row["temperature_k"]) for row in zhuang_ramp]) == (
        pytest.approx(
            4.22e-7 * ramp_pressure**3
            - 0.00155 * ramp_pressure**2
            + 2.58948 * ramp_pressure
            + 920.817,
            abs=5e-9,
        )
    )


def test_zhang_dataset_resource_and_record_metadata():
    path = _dataset_path()
    assert hashlib.sha256(path.read_bytes()).hexdigest() == EXPECTED_SHA256
    document = get_material_document("iron")
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    assert dataset["resource"]["sha256"] == EXPECTED_SHA256
    assert dataset["license"] == "CC BY 4.0"
    assert dataset["provenance"]["source_workbook_sha256"] == (
        "a3654047a5cabf03633c9ce280a264fa2c7bfae1c86d4e4c04028e715a0cd4a2"
    )
    quality_control = dataset["provenance"]["quality_control"]
    assert {
        key: quality_control[key]
        for key in (
            "rows",
            "table_s1_rows",
            "table_s3_rows",
            "table_s4_rows",
            "table_s2_excluded_rows",
            "s1_corrected_pressure_rows",
            "s1_retained_pressure_rows",
        )
    } == {
        "rows": 1313,
        "table_s1_rows": 1078,
        "table_s3_rows": 132,
        "table_s4_rows": 103,
        "table_s2_excluded_rows": 325,
        "s1_corrected_pressure_rows": 641,
        "s1_retained_pressure_rows": 437,
    }
    relations = quality_control["inferred_dynamic_temperature_relations"]
    assert relations["status"] == (
        "exact_relationship_in_deposited_values_not_published_methodology"
    )
    assert relations["zhuang_hugoniot_rows"] == 62
    assert relations["zhuang_ramp_rows"] == 60

    records = [
        item for item in document["eos_records"] if item["identifier"] in RECORD_IDS
    ]
    assert {item["identifier"] for item in records} == RECORD_IDS
    for record in records:
        assert record["fit_datasets"] == [DATASET_ID]
        assert record["validity"]["pressure_gpa"] == [10.7, 1374.0]
        assert record["validity"]["temperature_k"] == [298.0, 12000.0]
        assert record["pressure_calibration"]["status"] == "partially_resolved"
        assert record["pressure_calibration"]["recalculation"]["status"] == (
            "missing_calibrant_observations"
        )
        assert record["fit_provenance"]["selection"]["included_rows"] == 1313
        assert record["fit_provenance"]["reproduction"]["scope_status"] == (
            "final_input_parity_upstream_reduction_partial"
        )
        assert record["scientific_validation"]["primary_data_check"]["status"] == (
            "bundled"
        )
        assert record["parameter_provenance"]["covariance"].startswith(
            "Supplementary Table S5"
        )
        assert record["thermal"]["fixed_parameters"] == ["Tr", "n"]
        assert record["thermal"]["gas_constant_j_mol_k"] == 8.314


def test_zhang_records_are_executable_and_independent_refit_passes():
    document = get_material_document("iron")
    material = Material.from_eosmat(document, record_identifiers=sorted(RECORD_IDS))
    for wrapped in material.eos_records:
        raw = next(
            item
            for item in document["eos_records"]
            if item["identifier"] == wrapped.identifier
        )
        v0 = raw["eos"]["parameters"]["V0"]
        assert wrapped.pressure(v0, 300.0) == pytest.approx(0.0, abs=1.0e-9)
        assert wrapped.pressure(v0, 2000.0) > 0.0

    root = Path(__file__).resolve().parents[1]
    completed = subprocess.run(
        [
            sys.executable,
            str(root / "scripts" / "reproduce_zhang_2025_hcp_iron.py"),
            "--check",
        ],
        cwd=root,
        check=True,
        capture_output=True,
        text=True,
    )
    outcome = json.loads(completed.stdout)
    assert outcome["fits"]["1"]["observations"] == 1313
    assert outcome["fits"]["1"]["deposited_residual_reconstruction_max_abs_gpa"] < 2e-6
