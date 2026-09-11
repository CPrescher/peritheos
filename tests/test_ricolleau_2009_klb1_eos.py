from __future__ import annotations

import csv
import hashlib
import importlib.util
import json
from collections import Counter
from pathlib import Path

import pytest

from peritheos import get_material_document

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts" / "audit_ricolleau_2009_klb1_eos.py"
RAW = ROOT / "peritheos/data/datasets/klb1-ricolleau-2009-table-s1.tsv"
CSV = ROOT / "peritheos/data/datasets/klb1-ricolleau-2009-table-s1-pvt.csv"
AUDIT = ROOT / "docs/data/ricolleau-2009-klb1-eos-refit.json"


@pytest.fixture(scope="module")
def reproduction():
    spec = importlib.util.spec_from_file_location("ricolleau_2009_klb1", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module.reproduce()


def test_official_table_s1_transcription_is_complete_and_lossless():
    assert hashlib.sha256(RAW.read_bytes()).hexdigest() == (
        "e3f504dcb8953dbd6ef0cc0d341090a05116d32c865b748b9c95173fb94e7aa1"
    )
    assert hashlib.sha256(CSV.read_bytes()).hexdigest() == (
        "9d33c2a4c8ca7b4e4bc9d83e8881d89d7a0f41f42a252150e1d2fb360d77380f"
    )
    with CSV.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 153
    assert Counter(row["temperature_k"] for row in rows) == {
        "300": 17,
        **Counter(
            row["temperature_k"] for row in rows if row["temperature_k"] != "300"
        ),
    }
    assert sum(row["temperature_k"] != "300" for row in rows) == 136
    assert sum(bool(row["pressure_medium_volume_2_a3"]) for row in rows) == 8
    assert all(
        row[column]
        for row in rows
        for column in row
        if column
        not in {
            "pressure_medium_volume_2_a3",
            "pressure_medium_volume_2_sigma_a3",
        }
    )
    assert rows[0]["temperature_k"] == "1600"
    assert rows[0]["pressure_gpa"] == "32.75"
    assert rows[0]["mg_perovskite_volume_sigma_a3"] == "0.053"
    assert rows[-1]["temperature_k"] == "300"
    assert rows[-1]["pressure_gpa"] == "98.36"
    assert rows[-1]["gold_volume_a3"] == "51.747"


def test_staged_spin_and_thermal_refits_are_reproducible(
    reproduction, assert_audit_close
):
    assert_audit_close(reproduction, json.loads(AUDIT.read_text(encoding="utf-8")))
    assert reproduction["observations"] == {
        "total": 153,
        "room_temperature": 17,
        "heated": 136,
        "secondary_pressure_medium_volumes": 8,
    }
    recalculation = reproduction["pressure_recalculation"]
    assert recalculation["reference_eos_record"] == "gold_fei_2007_vinet_2"
    assert recalculation["observations"] == 153
    assert recalculation["recalculated_minus_published_mean_gpa"] == pytest.approx(
        0.6688880182
    )
    assert recalculation["recalculated_minus_published_rmse_gpa"] == pytest.approx(
        0.7013335207
    )
    assert recalculation["recalculated_minus_published_max_abs_gpa"] == pytest.approx(
        0.9410996839
    )
    selection = reproduction["spin_selection"]
    assert selection["high_spin_limiting_branch"]["observations"] == 7
    assert selection["unassigned_crossover"]["pressures_gpa"] == [
        52.22,
        53.0,
        54.89,
    ]
    assert selection["low_spin_limiting_branch"]["observations"] == 7

    refits = reproduction["record_refits"]
    high = refits["klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt"][
        "static_fit"
    ]
    low = refits["klb1_ferropericlase_ricolleau_2009_low_spin_300k_bm2"]["fit"]
    assert high["parameters"]["V0"] == pytest.approx(76.4382790)
    assert high["standard_errors"]["V0"] == pytest.approx(0.0177118)
    assert low["parameters"]["V0"] == pytest.approx(74.0428803)
    assert low["standard_errors"]["V0"] == pytest.approx(0.0253752)

    mg = refits["klb1_mg_perovskite_ricolleau_2009_bm2_alphakt"]["fit"]
    ca = refits["klb1_ca_perovskite_ricolleau_2009_bm2_alphakt"]["fit"]
    fp = refits["klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt"][
        "thermal_fit"
    ]
    assert mg["parameters"] == pytest.approx(
        {
            "K0": 245.1917171,
            "dK_dT": -0.03637101,
            "alpha0": 3.4591171e-5,
            "alpha1": 6.5642355e-9,
        },
        rel=5e-6,
    )
    assert ca["parameters"] == pytest.approx(
        {
            "K0": 243.6059204,
            "dK_dT": -0.03512956,
            "alpha0": 3.4593130e-5,
            "alpha1": 5.6142331e-9,
        },
        rel=5e-6,
    )
    assert fp["parameters"] == pytest.approx(
        {
            "dK_dT": -0.03401873,
            "alpha0": 2.7597290e-5,
            "alpha1": 3.0874780e-8,
        },
        rel=5e-6,
    )
    assert "row exclusions" in reproduction["conclusion"]["irreducible_blocker"]


def test_four_records_register_source_rows_and_recalculation_readiness():
    expected = {
        "klb1_mg_perovskite": 1,
        "klb1_ca_perovskite": 1,
        "klb1_ferropericlase": 2,
    }
    for material_identifier, count in expected.items():
        document = get_material_document(material_identifier)
        assert len(document["datasets"]) == 1
        dataset = document["datasets"][0]
        assert dataset["resource"]["path"].endswith(
            "klb1-ricolleau-2009-table-s1-pvt.csv"
        )
        assert len(dataset["columns"]) == 16
        assert len(document["eos_records"]) == count
        for record in document["eos_records"]:
            assert record["fit_datasets"] == [dataset["identifier"]]
            assert record["fit_provenance"]["dataset"] == dataset["identifier"]
            assert record["pressure_calibration"]["recalculation"]["status"] == (
                "ready"
            )
            check = record["scientific_validation"]["primary_data_check"]
            assert check["status"] == "bundled"
            assert check["dataset_identifiers"] == [dataset["identifier"]]


def test_common_ledger_reports_two_parity_and_two_similar_results():
    ledger = json.loads(
        (ROOT / "docs/data/primary-eos-refits.json").read_text(encoding="utf-8")
    )
    outcomes = {
        row["record_identifier"]: row
        for row in ledger["records"]
        if row["record_identifier"].startswith("klb1_")
    }
    assert {identifier: row["status"] for identifier, row in outcomes.items()} == {
        "klb1_ca_perovskite_ricolleau_2009_bm2_alphakt": "similar",
        "klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt": "similar",
        "klb1_ferropericlase_ricolleau_2009_low_spin_300k_bm2": "parity",
        "klb1_mg_perovskite_ricolleau_2009_bm2_alphakt": "parity",
    }
    assert (
        outcomes["klb1_ferropericlase_ricolleau_2009_high_spin_bm2_alphakt"][
            "parameters"
        ][1]["within_combined_2sigma"]
        is False
    )
