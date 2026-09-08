import hashlib
import json
from pathlib import Path

import pytest

from peritheos import get_material_document
from scripts.reproduce_anderson_1989_gold_thermal_eos import reproduce

ROOT = Path(__file__).resolve().parents[1]
REPORT = ROOT / "docs" / "data" / "anderson-1989-gold-thermal-eos-reproduction.json"
RECORD_IDENTIFIER = "gold_anderson_1989_bm3_1"
EXPECTED_DATASETS = {
    "gold_anderson_1989_table1_thermodynamic_inputs": (
        20,
        "6a460edad1925133feb9edbf2d0bad9d27dd1677a2e14821fb0f3a65f174253e",
    ),
    "gold_anderson_1989_table2_temperature_derivatives": (
        10,
        "be8233887659814ea16d52910b5f5772f2761a1c65b69191bead2096b4784e6e",
    ),
    "gold_anderson_1989_table3_bulk_moduli": (
        20,
        "aaa2e6f0e18eb927d4b4ee3c6d0b90635c041d9aebcc3e1a53f44f36f6320d06",
    ),
    "gold_anderson_1989_table4_anharmonic_diagnostics": (
        19,
        "348a86c38ad512ad1716ac62f1cda3dced28f9e17efc280632b373d78db0c18f",
    ),
    "gold_anderson_1989_table5_pressure_grid": (
        126,
        "08ca602e42c8e56ce47c78dba3057b9b1e0cd219de122b25351e5b6db5908be7",
    ),
}


def _record_and_datasets():
    document = get_material_document("gold")
    record = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == RECORD_IDENTIFIER
    )
    datasets = {
        item["identifier"]: item
        for item in document["datasets"]
        if item["identifier"] in EXPECTED_DATASETS
    }
    return record, datasets


def test_anderson_reproduction_report_is_current():
    assert json.loads(REPORT.read_text(encoding="utf-8")) == reproduce()


def test_anderson_source_tables_are_complete_and_checksummed():
    record, datasets = _record_and_datasets()
    assert set(datasets) == set(EXPECTED_DATASETS)
    for identifier, (expected_rows, expected_hash) in EXPECTED_DATASETS.items():
        dataset = datasets[identifier]
        path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
        assert dataset["resource"]["sha256"] == expected_hash
        assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_hash
        with path.open(encoding="utf-8") as stream:
            assert sum(1 for _ in stream) - 1 == expected_rows
        assert dataset["used_by_eos_records"] == [RECORD_IDENTIFIER]

    check = record["scientific_validation"]["primary_data_check"]
    assert check["status"] == "bundled_indirect"
    assert check["dataset_identifiers"] == list(EXPECTED_DATASETS)
    assert "fit_datasets" not in record


def test_anderson_300k_cp_print_is_preserved_but_not_silently_used():
    audit = reproduce()["table1_300k_cp_audit"]
    assert audit["printed_0p1_j_g_k"] == pytest.approx(1.1288)
    assert audit["equation4_0p1_j_g_k"] == pytest.approx(1.28820, abs=5e-5)
    assert audit["equation12_inverted_0p1_j_g_k"] == pytest.approx(1.28843, abs=5e-5)
    assert audit["gamma_from_printed_cp"] == pytest.approx(3.39458, abs=5e-5)
    assert audit["gamma_from_reconstructed_cp"] == pytest.approx(
        audit["table4_gamma"], abs=0.001
    )


def test_anderson_staged_regressions_recover_published_coefficients():
    report = reproduce()
    ambient = report["table2_ambient_kt_temperature_derivatives"]
    assert ambient["low_temperature"]["slope_gpa_k"] == pytest.approx(
        ambient["low_temperature"]["published_slope_gpa_k"], abs=6e-4
    )
    assert ambient["high_temperature"]["slope_gpa_k"] == pytest.approx(
        ambient["high_temperature"]["published_slope_gpa_k"], abs=6e-4
    )

    thermal = report["table4_thermal_pressure_linear_fit"]
    assert thermal["slope_gpa_k"] == pytest.approx(0.00714, abs=5e-6)
    assert thermal["intercept_gpa"] == pytest.approx(-0.40, abs=0.005)

    trials = report["table3_constant_volume_linear_fits"]["trials"]
    for trial in trials.values():
        assert trial["slope_gpa_k"] == pytest.approx(
            trial["published_slope_gpa_k"], abs=7e-5
        )
    assert trials["5.5"]["slope_gpa_k"] == pytest.approx(-0.0115, abs=2e-5)


def test_anderson_equation29_reproduces_every_table5_top_row():
    report = reproduce()
    table5 = report["table5_equation29_reproduction"]
    assert table5["states"] == 126
    assert table5["max_abs_residual_gpa"] <= table5["rounding_tolerance_gpa"]
    assert table5["rmse_gpa"] < 0.003


def test_anderson_global_fit_is_explicitly_not_defined():
    record, _ = _record_and_datasets()
    report = reproduce()
    conclusion = report["conclusion"]
    assert conclusion["global_fit_status"] == "not_defined_by_source"
    assert conclusion["classification"] == "not_refittable"
    assert record["pressure_calibration"]["status"] == "not_applicable"
    assert record["scientific_validation"]["reproduction"]["table5_states"] == 126
    assert "global P-V-T regression" in record["notes"]
