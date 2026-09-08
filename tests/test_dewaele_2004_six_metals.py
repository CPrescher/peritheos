import csv
import hashlib
import math
from pathlib import Path

import pytest

from peritheos import get_material_document
from scripts.validate_primary_eos_refits import validate_all

ROOT = Path(__file__).resolve().parents[1]
DOI = "10.1103/physrevb.70.094112"
SOURCE_PDF_SHA256 = "1f9e02e4a77f49a868b1a6261290465bdb00597b1660691b1790c7502ec9aeeb"
CASES = {
    "aluminum": (
        40,
        4,
        "61fa04908b3a4827fb6566e4ceea29706ac52cbaae1cee5bb0590d1af44423b7",
    ),
    "copper": (
        42,
        4,
        "0f4535c7e29960690c098408711f63ae559be1d73bba18da73047d75f266120e",
    ),
    "gold": (37, 4, "c3a52a211b94cef2cd46198577e18efe4fc66befd7d38aa4f20063d290e97438"),
    "platinum": (
        36,
        4,
        "b0230e9fdd1f8a36cecc35cb79a6d56bd4ee170007a465f1842ef76544c47273",
    ),
    "tantalum": (
        36,
        2,
        "3c69e1a9943dc1f882d60f63d01cf6708c114a4569a20e04a8b2278e5c09bd78",
    ),
    "tungsten": (
        42,
        2,
        "d73f1d9114c974be80175c6ade256bf17cb4a529e911480e6cd1e337394ef523",
    ),
}
NEW_RECORDS = {
    "aluminum_dewaele_2004_mao_ruby_vinet": ([0.0, 144.3], ["V0"]),
    "copper_dewaele_2004_mao_ruby_vinet": ([0.0, 144.3], ["V0"]),
    "gold_dewaele_2004_mao_ruby_vinet": ([0.0, 90.0], ["V0"]),
    "platinum_dewaele_2004_mao_ruby_vinet": ([0.0, 90.0], ["V0"]),
    "platinum_dewaele_2004_revised_ruby_vinet": (
        [0.0, 93.6],
        ["V0", "K0"],
    ),
    "tantalum_dewaele_2004_mao_ruby_vinet": ([0.0, 90.0], ["V0"]),
    "tantalum_dewaele_2004_revised_ruby_vinet": (
        [0.0, 93.6],
        ["V0", "K0"],
    ),
    "tungsten_dewaele_2004_mao_ruby_vinet": ([0.0, 144.3], ["V0"]),
}
EXISTING_REVISED_RECORDS = {
    "aluminum_dewaele_2004_vinet_1": ([0.0, 153.0], ["V0"]),
    "copper_dewaele_2004_vinet_1": ([0.0, 153.0], ["V0"]),
    "gold_dewaele_2004_vinet_5": ([0.0, 93.6], ["V0", "K0"]),
    "tungsten_dewaele_2004_vinet_2": ([0.0, 153.0], ["V0"]),
}


def _dataset(material_identifier):
    document = get_material_document(material_identifier)
    dataset_id = f"{material_identifier}_dewaele_2004_table1_compression"
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == dataset_id
    )
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, dataset, path, rows


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_dewaele_table1_transcriptions_are_complete_and_checksummed(
    material_identifier, case
):
    expected_rows, _, expected_sha256 = case
    _, dataset, path, rows = _dataset(material_identifier)

    assert len(rows) == expected_rows
    assert hashlib.sha256(path.read_bytes()).hexdigest() == expected_sha256
    assert dataset["resource"]["sha256"] == expected_sha256
    assert dataset["transcription"]["source_artifact"]["sha256"] == (SOURCE_PDF_SHA256)
    assert "not asserted to be openly licensed" in dataset["license"]
    assert set(rows[0]) == {
        "ruby_pressure_classical_gpa",
        "ruby_pressure_revised_gpa",
        "atomic_volume_a3",
        "atomic_volume_uncertainty_a3",
    }
    assert {float(row["atomic_volume_uncertainty_a3"]) for row in rows} == {0.01}


def test_dewaele_table1_preserves_printed_blanks_and_pressure_pair_anomaly():
    _, _, _, platinum = _dataset("platinum")
    _, _, _, tantalum = _dataset("tantalum")
    _, _, _, copper = _dataset("copper")

    assert not any(
        float(row["ruby_pressure_classical_gpa"]) == 15.04 for row in platinum
    )
    assert not any(
        float(row["ruby_pressure_classical_gpa"]) == 78.0 for row in tantalum
    )
    anomalous = next(
        row for row in copper if float(row["ruby_pressure_classical_gpa"]) == 37.0
    )
    assert anomalous["ruby_pressure_revised_gpa"] == "37.1"

    classical = float(anomalous["ruby_pressure_classical_gpa"])
    wavelength_ratio = (1.0 + 7.665 * classical / 1904.0) ** (1.0 / 7.665)
    recalculated = 1904.0 / 9.5 * (wavelength_ratio**9.5 - 1.0)
    assert recalculated == pytest.approx(37.6362915365)
    assert abs(float(anomalous["ruby_pressure_revised_gpa"]) - recalculated) > 0.5


@pytest.mark.parametrize(("material_identifier", "case"), CASES.items())
def test_dewaele_records_preserve_cell_normalization_and_source_metadata(
    material_identifier, case
):
    _, formula_units_per_cell, _ = case
    document, dataset, _, _ = _dataset(material_identifier)
    assert document["formula_units_per_cell"] == formula_units_per_cell

    expected_ids = {
        identifier
        for identifier in NEW_RECORDS | EXISTING_REVISED_RECORDS
        if identifier.startswith(f"{material_identifier}_")
    }
    if material_identifier != "gold":
        expected_ids.update(
            {
                f"{material_identifier}_dewaele_2019_mao_vinet",
                f"{material_identifier}_dewaele_2019_dor_vinet",
            }
        )
    assert set(dataset["used_by_eos_records"]) == expected_ids
    records = {
        record["identifier"]: record
        for record in document["eos_records"]
        if record["identifier"] in NEW_RECORDS | EXISTING_REVISED_RECORDS
    }
    for identifier, record in records.items():
        expected_range, fixed = (NEW_RECORDS | EXISTING_REVISED_RECORDS)[identifier]
        assert record["reference"]["doi"].lower() == DOI
        assert record["temperature_ref"] == 298.0
        assert record["experimental_pressure_range_gpa"] == expected_range
        assert record["validity"]["pressure_gpa"] == expected_range
        assert record["fixed_parameters"] == fixed
        assert record["parameter_error_confidence"] == 0.95
        assert record["pressure_calibration"]["recalculation"]["status"] == "ready"
        atomic_v0 = record["eos"]["parameters"]["V0"] / formula_units_per_cell
        expected_atomic_v0 = {
            "aluminum": 16.573,
            "copper": 11.810,
            "gold": 16.962,
            "platinum": 15.095,
            "tantalum": 18.035,
            "tungsten": 15.862,
        }[material_identifier]
        assert atomic_v0 == pytest.approx(expected_atomic_v0)


def test_all_twelve_dewaele_source_fits_reach_uncertainty_parity():
    ledger = validate_all()
    expected_ids = set(NEW_RECORDS | EXISTING_REVISED_RECORDS)
    results = {
        item["record_identifier"]: item
        for item in ledger["records"]
        if item["record_identifier"] in expected_ids
    }

    assert set(results) == expected_ids
    assert {item["status"] for item in results.values()} == {"parity"}
    assert all(
        comparison["within_combined_2sigma"]
        for item in results.values()
        for comparison in item["parameters"]
    )
    assert results["platinum_dewaele_2004_mao_ruby_vinet"]["observations"] == 36
    assert results["tantalum_dewaele_2004_mao_ruby_vinet"]["observations"] == 36
    assert results["gold_dewaele_2004_mao_ruby_vinet"]["columns"]["pressure"] == (
        "ruby_pressure_classical_gpa"
    )
    assert results["gold_dewaele_2004_vinet_5"]["columns"]["pressure"] == (
        "ruby_pressure_revised_gpa"
    )


def test_printed_pressure_columns_follow_the_two_ruby_scales_with_rounding_caveat():
    for material_identifier in CASES:
        _, _, _, rows = _dataset(material_identifier)
        differences = []
        for row in rows:
            classical = float(row["ruby_pressure_classical_gpa"])
            revised = float(row["ruby_pressure_revised_gpa"])
            ratio = (1.0 + 7.665 * classical / 1904.0) ** (1.0 / 7.665)
            from_classical = 1904.0 / 9.5 * (ratio**9.5 - 1.0)
            differences.append(revised - from_classical)
        assert all(math.isfinite(value) for value in differences)
        assert max(map(abs, differences)) < 0.54
