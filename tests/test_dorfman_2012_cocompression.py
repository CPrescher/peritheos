import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import pytest

from peritheos import get_material_document
from scripts.reproduce_dorfman_2012_cocompression import (
    DATASET_LICENSE,
    SOURCE_SHA256,
    extract_observations,
    fit,
    load_dataset,
    vinet_pressure,
)

ROOT = Path(__file__).resolve().parents[1]
SUMMARY_PATH = ROOT / "docs" / "data" / "dorfman-2012-cocompression-refit.json"
DATASET_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "dorfman-2012-tables-s1-s6-cocompression.csv"
)
LICENSE_PATH = DATASET_PATH.with_suffix(".LICENSE.md")
DATASET_ID = "dorfman_2012_tables_s1_s6_cocompression"


def load_summary():
    return json.loads(SUMMARY_PATH.read_text(encoding="utf-8"))


def test_dorfman_summary_records_source_custody_and_complete_available_scope():
    summary = load_summary()

    assert summary["source"]["doi"] == "10.1029/2012JB009292"
    assert summary["source"]["correction_doi"] == "10.1029/2012JB009800"
    assert summary["source"]["supporting_sha256"] == SOURCE_SHA256
    assert "publisher PDF" in summary["source"]["license"]
    assert summary["dataset_release"]["license"] == DATASET_LICENSE
    assert (
        summary["dataset_release"]["sha256"]
        == hashlib.sha256(DATASET_PATH.read_bytes()).hexdigest()
    )
    assert "third-party rights are excluded" in summary["dataset_release"]["scope"]
    assert summary["reproduction"]["status"] == (
        "source_rows_refitted_but_published_coefficients_not_reproduced"
    )
    assert summary["selection"]["observation_rows"] == 165
    assert summary["selection"]["volume_values"] == 368
    assert summary["selection"]["pair_count"] == 241
    assert summary["selection"]["single_peak_mgo_rows_included"] == 3
    assert "AN012" in summary["selection"]["missing_run"]
    assert "observations" not in summary


def test_dorfman_bundled_csv_is_complete_lossless_and_refittable():
    with DATASET_PATH.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))

    assert len(rows) == 368
    assert len({row["observation_index"] for row in rows}) == 165
    assert Counter(row["material"] for row in rows) == {
        "Au": 22,
        "MgO": 108,
        "Mo": 17,
        "NaCl B2": 80,
        "Ne": 44,
        "Pt": 97,
    }
    assert sum(row["single_peak_mgo"] == "1" for row in rows) == 3
    assert rows[0]["volume_source_token"] == "26.35(2)"
    assert rows[0]["reported_pressure_source_token"] == "61.9"
    assert rows[-1]["volume_source_token"] == "17.46(7)"

    observations = load_dataset(DATASET_PATH)
    summary = load_summary()
    assert len(observations) == 165
    for mode in ("fixed", "free"):
        replay = fit(observations, mode)
        assert replay["objective"] == pytest.approx(summary["fits"][mode]["objective"])
        for material, parameters in replay["parameters"].items():
            assert parameters == pytest.approx(
                summary["fits"][mode]["parameters"][material]
            )

    license_notice = LICENSE_PATH.read_text(encoding="utf-8")
    assert "CC0 1.0 Universal" in license_notice
    assert "does not apply to" in license_notice
    assert "publisher's auxiliary PDF" in license_notice


def test_dorfman_coupled_refits_converge_without_using_printed_pressures():
    summary = load_summary()

    assert "never as fit observations" in summary["reproduction"]["interpretation"]
    assert summary["published_pressure_transcription_check"][
        "max_abs_difference_gpa"
    ] == pytest.approx(0.18238546658528776)
    for mode in ("fixed", "free"):
        fit = summary["fits"][mode]
        assert fit["success"]
        assert fit["pair_count"] == 241
        assert fit["objective"] < fit["published_objective"]

    assert summary["fits"]["fixed"]["parameters"]["Mo"]["K0_prime"] == pytest.approx(
        4.13809805589023
    )
    assert summary["fits"]["free"]["parameters"]["Au"]["K0"] == pytest.approx(
        178.7131955926287
    )
    assert summary["fits"]["free"]["parameters"]["Pt"]["K0_prime"] == pytest.approx(
        4.898844953030691
    )


def test_dorfman_reproducer_rejects_an_unverified_source_file(tmp_path):
    source = tmp_path / "support.pdf"
    source.write_bytes(b"not the official source")

    with pytest.raises(ValueError, match="unexpected supporting-PDF SHA-256"):
        extract_observations(source)


def test_dorfman_records_link_coupled_refit_and_tange_anchor():
    expected = {
        "gold": "Au",
        "molybdenum": "Mo",
        "platinum": "Pt",
    }
    for material, source_name in expected.items():
        document = get_material_document(material)
        records = [
            record
            for record in document["eos_records"]
            if "_dorfman_2012_tange_mgo_k0_" in record["identifier"]
        ]
        assert len(records) == 2
        for record in records:
            check = record["scientific_validation"]["primary_data_check"]
            method = record["pressure_calibration"]["methods"][0]
            assert check["status"] == "bundled"
            assert check["dataset_identifiers"] == [DATASET_ID]
            assert record["fit_datasets"] == [DATASET_ID]
            assert check["pair_count"] == 241
            assert check["objective"] < check["published_objective"]
            assert source_name in record["label"]
            assert set(check["record_parameter_comparison"]) <= {"K0", "K0_prime"}
            assert record["pressure_calibration"]["recalculation"]["status"] == "ready"
            assert method["reference"]["doi"] == "10.1029/2008JB005813"
            assert method["reference_eos_record"] == "mgo_b1_tange_2009_vinet"
            assert record["pressure_range_status"] == "reference_parameterization"

        dataset = next(
            item for item in document["datasets"] if item["identifier"] == DATASET_ID
        )
        assert dataset["license"] == "CC0-1.0"
        assert "excludes the article" in dataset["license_scope"]
        assert (
            dataset["resource"]["sha256"]
            == hashlib.sha256(DATASET_PATH.read_bytes()).hexdigest()
        )
        assert dataset["used_by_eos_records"] == [
            f"{material}_dorfman_2012_tange_mgo_k0_fixed_vinet",
            f"{material}_dorfman_2012_tange_mgo_k0_free_vinet",
        ]

    ledger = json.loads(
        (ROOT / "docs" / "data" / "primary-eos-refits.json").read_text(encoding="utf-8")
    )
    outcomes = [
        item
        for item in ledger["records"]
        if "_dorfman_2012_tange_mgo_k0_" in item["record_identifier"]
    ]
    assert len(outcomes) == 6
    assert {item["status"] for item in outcomes} == {"parity_not_achieved"}
    assert {item["comparison_scope"] for item in outcomes} == {"joint_system"}


def test_dorfman_vinet_reference_state_is_zero_pressure():
    parameters = {"V0": 67.85, "K0": 167.0, "K0_prime": 5.88}
    assert vinet_pressure(parameters["V0"], parameters) == pytest.approx(0.0)
