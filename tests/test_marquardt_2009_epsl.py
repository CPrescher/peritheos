import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material
from scripts.reproduce_marquardt_2009_epsl import (
    DATASET,
    LS_RECORDS,
    RECORD,
    bm3_pressure,
    reproduce,
    reproduce_supplement,
)

ROOT = Path(__file__).resolve().parents[1]
DOI = "10.1016/j.epsl.2009.08.017"


def document():
    return json.loads(
        (ROOT / "peritheos/data/materials/mg090fe010o.eosmat").read_text()
    )


def test_epsl_source_has_one_hs_and_twelve_ls_owned_published_records():
    records = [
        record
        for path in (ROOT / "peritheos/data/materials").glob("*.eosmat")
        for record in json.loads(path.read_text())["eos_records"]
        if record["reference"].get("doi") == DOI
    ]
    assert [r["identifier"] for r in records] == [RECORD, *LS_RECORDS]
    record = records[0]
    assert record["record_kind"] == "published"
    assert record["eos"]["type"] == "BM3"
    assert record["equation_kind"] == "isothermal"
    assert record["eos"]["parameters"] == dict(V0=75.62, K0=158.2, K0_prime=3.98)
    assert record["parameter_errors"] == dict(V0=0.06, K0=2, K0_prime=0.14)
    assert record["fixed_parameters"] == ["V0"]
    assert record["volume_basis"]["formula_units"] == 4
    assert record["experimental_pressure_range_gpa"] == [0.0001, 43.8]
    assert record["spin_crossover_context"]["mixed_spin_pressure_gpa"] == [45, 63]
    assert record["fit_datasets"] == [DATASET]
    assert "adiabatic" in record["parameter_provenance"]["Brillouin_constraint"]
    assert record["pressure_calibration"]["status"] == "partially_resolved"
    audit = json.loads((ROOT / "peritheos/data/primary-source-audit.json").read_text())
    assert [r["record"] for r in audit["records"] if r.get("doi") == DOI] == [
        RECORD,
        *LS_RECORDS,
    ]


def test_complete_table2_and_source_specific_selection():
    dataset = document()["datasets"][0]
    resource = dataset["resource"]
    path = ROOT / "peritheos/data" / resource["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == resource["sha256"]
    assert (
        resource["sha256"]
        == "ac504ad2b8175dc1c0f199309cc1c143ec2744e4299305f699dd9a4410cb3b24"
    )
    assert dataset["reference"]["doi"] == DOI
    assert dataset["used_by_eos_records"] == [RECORD, *LS_RECORDS]
    with path.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 29
    assert Counter(r["spin_region"] for r in rows) == dict(
        high_spin=14, mixed_spin=9, low_spin=6
    )
    assert [int(r["source_order"]) for r in rows] == list(range(1, 30))
    for row in rows:
        p = float(row["pressure_gpa"])
        assert row["spin_region"] == (
            "high_spin" if p < 45 else "low_spin" if p > 63 else "mixed_spin"
        )
    assert [float(rows[i]["pressure_gpa"]) for i in [0, 13, 14, 22, 23, 28]] == [
        0.0001,
        43.8,
        45.7,
        62.7,
        63.9,
        77.4,
    ]


def test_published_curve_matches_independent_bm3_and_validation_stays_separate():
    eos = Material.from_eosmat(document()).get_eos_record(RECORD)
    volumes = np.array([75.62, 72.78, 67.22, 62.72])
    assert np.asarray(eos.pressure(volumes)) == pytest.approx(
        bm3_pressure(volumes, 75.62, 158.2, 3.98)
    )
    result = reproduce()
    assert result["solver_success"]
    assert result["source_orders"] == list(range(1, 15))
    assert result["purpose"] == "validation_only_pv_refit_without_Brillouin_constraint"
    assert result["published_pressure_rmse_gpa"] == pytest.approx(0.7158442283)
    assert result["published_volume_rmse_a3"] == pytest.approx(0.2286494131)
    assert result["validation_parameters"]["K0"] == pytest.approx(138.24256, abs=0.001)
    assert result["validation_parameters"]["K0_prime"] == pytest.approx(
        5.65047, abs=0.0001
    )
    # Later coupled fits retain their own source and values.
    later = [
        r
        for r in document()["eos_records"]
        if r["reference"]["doi"] == "10.2138/am-2016-5510"
    ]
    assert [r["eos"]["parameters"] for r in later] == [
        dict(V0=75.55, K0=159, K0_prime=3.96),
        dict(V0=74.59, K0=159, K0_prime=4),
    ]


def test_ls_family_preserves_all_source_triplets_and_reproduces_kt():
    expected = [
        (74.85, 144, 4.5),
        (74.10, 153, 4.5),
        (73.34, 163, 4.5),
        (72.59, 175, 4.5),
        (74.85, 156, 4.0),
        (74.10, 164, 4.0),
        (73.34, 175, 4.0),
        (72.59, 186, 4.0),
        (74.85, 167, 3.5),
        (74.10, 177, 3.5),
        (73.34, 187, 3.5),
        (72.59, 199, 3.5),
    ]
    doc = document()
    records = {r["identifier"]: r for r in doc["eos_records"]}
    material = Material.from_eosmat(doc)
    for index, (rid, parameters) in enumerate(zip(LS_RECORDS, expected), 1):
        r = records[rid]
        assert (
            tuple(r["eos"]["parameters"][key] for key in ["V0", "K0", "K0_prime"])
            == parameters
        )
        assert r["record_kind"] == "published"
        assert r["fixed_parameters"] == ["V0", "K0_prime"]
        assert r["parameter_errors"] == dict(V0=None, K0=None, K0_prime=None)
        assert r["volume_basis"]["formula_units"] == 4
        assert r["experimental_pressure_range_gpa"] == [63.9, 77.4]
        assert r["fit_datasets"] == [DATASET]
        assert r["constrained_fit_family"]["source_row"] == index
        assert not r["constrained_fit_family"]["preferred_member"]
        assert r["reference"]["doi"] == DOI
        volumes = np.array([58.51, 57.30, 56.26])
        assert np.asarray(
            material.get_eos_record(rid).pressure(volumes)
        ) == pytest.approx(bm3_pressure(volumes, *parameters))
    result = reproduce_supplement()
    assert result["source_orders"] == list(range(24, 30))
    assert len(result["models"]) == 12
    for model in result["models"]:
        assert model["solver_success"]
        assert abs(model["validation_k0_gpa"] / model["published_k0_gpa"] - 1) < 0.01
        for point in model["ls_kt_checkpoints"]:
            assert abs(point["calculated_kt_gpa"] - point["printed_kt_gpa"]) < 0.5
    assert result["ambient_table_s2_kt_conversion_gpa"] == pytest.approx(162.6794364)


def test_supplement_tables_complete_checksummed_and_not_pv_inputs():
    datasets = {d["identifier"]: d for d in document()["datasets"]}
    for suffix, count in [("s1_ls_models", 15), ("s2_elasticity", 26)]:
        d = datasets["mg090fe010o_marquardt_2009_table_" + suffix]
        resource = ROOT / "peritheos/data" / d["resource"]["path"]
        assert (
            hashlib.sha256(resource.read_bytes()).hexdigest() == d["resource"]["sha256"]
        )
        assert (
            d["source_document"]["sha256"]
            == "78c3a49a41795450b2ae39e97882268fab7da08cc89630ae0eaaded026160f04"
        )
        assert len(d["source_document"]["source_filenames"]) == 2
        with resource.open(newline="") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == count
        if suffix.startswith("s1"):
            assert [r["source_row"] for r in rows] == [str(i) for i in range(1, 13)] + [
                "min",
                "max",
                "avg",
            ]
            assert rows[-1]["v0_ls_a3_conventional_cell"] == "73.72"
            assert rows[-1]["k0_gpa"] == "171"
        else:
            assert rows[0]["ks_gpa_printed"] == "165(3)"
            assert rows[0]["ks_gpa_uncertainty"] == "3"
            assert rows[-1]["pressure_gpa"] == "81.2"
            assert rows[-1]["vbulk_km_s_printed"] == "9.48(27)"
            assert float(rows[-1]["vbulk_km_s_uncertainty"]) == 0.27
    assert all(
        r["fit_datasets"] == [DATASET]
        for r in document()["eos_records"]
        if r["reference"]["doi"] == DOI
    )


def test_global_refits_select_six_ls_rows_and_do_not_claim_uncertainty_parity():
    ledger = json.loads((ROOT / "docs/data/primary-eos-refits.json").read_text())
    records = {r["record_identifier"]: r for r in ledger["records"]}
    assert records[RECORD]["status"] == "not_refittable"
    for rid in LS_RECORDS:
        r = records[rid]
        assert r["status"] == "similar"
        assert r["observations"] == 6
        assert r["selection"] == "spin_region=low_spin"
        assert r["fixed_parameters"] == ["V0", "K0_prime"]
        assert r["free_parameters"] == ["K0"]
        assert r["parameters"][0]["published_error"] is None
