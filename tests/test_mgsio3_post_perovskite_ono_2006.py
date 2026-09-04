import csv
import hashlib
import math
import re
from collections import Counter
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.eosmat import validate_pressure_calibration_references

ROOT = Path(__file__).parents[1]
DOI = "10.2138/am.2006.2118"
DATASET_ID = "mgsio3_post_perovskite_ono_2006_table2_compression"
RESOURCE = "mgsio3-post-perovskite-ono-2006-table2-compression.csv"
SHA256 = "6ebd1e402c43d31054990344aa387730b488ef702b152d57f3376b940e49042d"
CASES = {
    "anderson": {
        "identifier": "mgsio3_post_perovskite_ono_2006_anderson_bm2_3",
        "pressure_column": "pressure_anderson_1989_gpa",
        "pressure_uncertainty_column": "pressure_anderson_1989_uncertainty_gpa",
        "k0": 237.0,
        "range": [116.1, 144.4],
        "rmse": 1.4316084628770336,
        "maximum_residual": 2.468353125587896,
        "eiv_k0": 237.18599456965256,
        "reference_eos": "gold_anderson_1989_bm3_1",
    },
    "jamieson": {
        "identifier": "mgsio3_post_perovskite_ono_2006_jamieson_bm2_4",
        "pressure_column": "pressure_jamieson_1982_gpa",
        "pressure_uncertainty_column": "pressure_jamieson_1982_uncertainty_gpa",
        "k0": 226.0,
        "range": [111.5, 137.0],
        "rmse": 1.74311488136343,
        "maximum_residual": 2.900201714695612,
        "eiv_k0": 226.81325417192315,
        "reference_eos": None,
    },
    "dewaele": {
        "identifier": "mgsio3_post_perovskite_ono_2006_dewaele_bm2_5",
        "pressure_column": "pressure_dewaele_2004_gpa",
        "pressure_uncertainty_column": "pressure_dewaele_2004_uncertainty_gpa",
        "k0": 248.0,
        "range": [121.4, 151.1],
        "rmse": 1.4859096031658499,
        "maximum_residual": 2.6365045364801745,
        "eiv_k0": 248.0464476576452,
        "reference_eos": "gold_dewaele_2004_vinet_5",
    },
}


def _document_records_and_rows():
    document = get_material_document("mgsio3_post_perovskite")
    records = {
        record["identifier"]: record
        for record in document["eos_records"]
        if "_ono_2006_" in record["identifier"]
    }
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = ROOT / "peritheos" / "data" / dataset["resource"]["path"]
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return document, records, dataset, path, rows


def _executable(document, identifier):
    return Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).get_eos_record(identifier)


def test_ono_2006_phase_cell_basis_and_distinct_published_records():
    document, records, _, _, _ = _document_records_and_rows()

    assert document["formula"] == "MgSiO3"
    assert "post-perovskite" in document["phase"]
    assert document["space_group"] == "Cmcm"
    assert document["space_group_number"] == 63
    assert document["formula_units_per_cell"] == 4

    contents = Counter()
    for site in document["atom_sites"]:
        multiplicity = int(re.match(r"\d+", site["wyckoff"]).group())
        contents[site["element"]] += multiplicity * site["occupancy"]
    assert contents == {"Mg": 4.0, "Si": 4.0, "O": 12.0}

    assert set(records) == {case["identifier"] for case in CASES.values()}
    assert (
        len(
            {tuple(record["eos"]["parameters"].values()) for record in records.values()}
        )
        == 3
    )
    assert all(record["reference"]["doi"].lower() == DOI for record in records.values())
    assert not any(
        record["reference"]["doi"].lower() == DOI
        for record in document["eos_records"]
        if record["identifier"] not in records
    )


@pytest.mark.parametrize("case", CASES.values(), ids=CASES)
def test_ono_2006_bm2_parameters_reference_state_and_scale_range(case):
    document, records, _, _, _ = _document_records_and_rows()
    source = records[case["identifier"]]
    record = _executable(document, case["identifier"])

    assert source["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 162.86, "K0": case["k0"]},
    }
    assert source["parameter_errors"] == {"V0": None, "K0": 1.0}
    assert source["parameter_error_confidence"] is None
    assert source["parameter_covariance"] is None
    assert source["fixed_parameters"] == ["V0"]
    assert source["temperature_ref"] == 300.0
    assert source["experimental_pressure_range_gpa"] == case["range"]
    assert source["fit_datasets"] == [DATASET_ID]

    method = source["pressure_calibration"]["methods"][0]
    assert method.get("reference_eos_record") == case["reference_eos"]
    assert source["scientific_validation"]["status"] == "primary_source_validated"

    assert record.pressure(162.86) == pytest.approx(0.0, abs=1.0e-12)
    assert record.eos.bulk_modulus(162.86) == pytest.approx(case["k0"])
    pressure_min, pressure_max = case["range"]
    for pressure in np.linspace(pressure_min + 1.0e-6, pressure_max - 1.0e-6, 5):
        volume = record.volume(pressure, check_validity=True)
        assert record.pressure(volume, check_validity=True) == pytest.approx(
            pressure, rel=1.0e-11
        )


def test_ono_2006_table2_transcription_checksum_and_printed_typo():
    _, records, dataset, path, rows = _document_records_and_rows()

    assert dataset["resource"]["path"].endswith(RESOURCE)
    assert hashlib.sha256(path.read_bytes()).hexdigest() == SHA256
    assert dataset["resource"]["sha256"] == SHA256
    assert len(rows) == 6
    assert len(rows[0]) == 19
    assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
    assert rows[0]["gold_lattice_a_angstrom"] == "3.6862"
    assert rows[0]["pressure_anderson_1989_gpa"] == "116.1"
    assert rows[-1]["pressure_dewaele_2004_gpa"] == "151.1"
    assert all(row["temperature_k"] == "300" for row in rows)
    assert all(row["fit_included"] == "1" for row in rows)
    assert dataset["used_by_eos_records"] == list(records)

    typographical_row = rows[2]
    assert typographical_row["mgsio3_lattice_a_angstrom"] == "2.250"
    lattice_product = math.prod(
        float(typographical_row[column])
        for column in (
            "mgsio3_lattice_a_angstrom",
            "mgsio3_lattice_b_angstrom",
            "mgsio3_lattice_c_angstrom",
        )
    )
    assert lattice_product == pytest.approx(110.2832595)
    assert float(typographical_row["mgsio3_unit_cell_volume_a3"]) == 120.11
    assert (
        "2.250(3)" in dataset["quality_control"]["reported_typographical_inconsistency"]
    )


@pytest.mark.parametrize("case", CASES.values(), ids=CASES)
def test_ono_2006_published_curve_and_diagnostic_refit(case):
    document, records, _, _, rows = _document_records_and_rows()
    source = records[case["identifier"]]
    record = _executable(document, case["identifier"])

    volume = np.array([float(row["mgsio3_unit_cell_volume_a3"]) for row in rows])
    volume_error = np.array(
        [float(row["mgsio3_unit_cell_volume_uncertainty_a3"]) for row in rows]
    )
    pressure = np.array([float(row[case["pressure_column"]]) for row in rows])
    pressure_error = np.array(
        [float(row[case["pressure_uncertainty_column"]]) for row in rows]
    )

    calculated = np.asarray(record.pressure(volume), dtype=float)
    residual = calculated - pressure
    assert math.sqrt(float(np.mean(residual**2))) == pytest.approx(
        case["rmse"], abs=5.0e-12
    )
    assert float(np.max(np.abs(residual))) == pytest.approx(
        case["maximum_residual"], abs=5.0e-12
    )

    v0 = 162.86
    basis = 1.5 * ((v0 / volume) ** (7.0 / 3.0) - (v0 / volume) ** (5.0 / 3.0))
    k0 = case["k0"]
    for _ in range(100):
        step = 1.0e-5
        basis_plus = 1.5 * (
            (v0 / (volume + step)) ** (7.0 / 3.0)
            - (v0 / (volume + step)) ** (5.0 / 3.0)
        )
        basis_minus = 1.5 * (
            (v0 / (volume - step)) ** (7.0 / 3.0)
            - (v0 / (volume - step)) ** (5.0 / 3.0)
        )
        derivative = k0 * (basis_plus - basis_minus) / (2.0 * step)
        effective_error = np.sqrt(pressure_error**2 + (derivative * volume_error) ** 2)
        next_k0 = np.sum(basis * pressure / effective_error**2) / np.sum(
            basis**2 / effective_error**2
        )
        if abs(next_k0 - k0) < 1.0e-13:
            break
        k0 = next_k0

    assert k0 == pytest.approx(case["eiv_k0"], abs=5.0e-10)
    stored = source["scientific_validation"]["independent_refit"]["parameters"]
    assert stored["K0"] == pytest.approx(k0, abs=5.0e-8)
    assert abs(k0 - case["k0"]) <= 1.0


def test_ono_2006_gold_calibrant_recalculation_and_global_links():
    document, records, _, _, rows = _document_records_and_rows()
    gold_document = get_material_document("gold")
    gold_a = np.array([float(row["gold_lattice_a_angstrom"]) for row in rows])

    checks = {
        "anderson": ("pressure_anderson_1989_gpa", 0.75),
        "dewaele": ("pressure_dewaele_2004_gpa", 0.07),
    }
    for scale, (column, tolerance) in checks.items():
        identifier = CASES[scale]["reference_eos"]
        gold = _executable(gold_document, identifier)
        recalculated = np.asarray(gold.pressure(gold_a**3, 300.0), dtype=float)
        published = np.array([float(row[column]) for row in rows])
        assert float(np.max(np.abs(recalculated - published))) < tolerance
        assert (
            records[CASES[scale]["identifier"]]["pressure_calibration"]["methods"][0][
                "reference_eos_record"
            ]
            == identifier
        )

    jamieson = records[CASES["jamieson"]["identifier"]]["pressure_calibration"]
    assert jamieson["status"] == "partially_resolved"
    assert jamieson["recalculation"]["status"] == "reference_eos_not_bundled"
    validate_pressure_calibration_references()


def test_ono_2006_source_inconsistencies_are_explicit():
    _, records, _, _, _ = _document_records_and_rows()
    anderson = records[CASES["anderson"]["identifier"]]["scientific_validation"]
    jamieson = records[CASES["jamieson"]["identifier"]]["scientific_validation"]
    dewaele = records[CASES["dewaele"]["identifier"]]["scientific_validation"]

    assert "236(1)" in anderson["reported_inconsistencies"][0]
    assert "225(1)" in jamieson["reported_inconsistencies"][0]
    assert "151.1 GPa" in dewaele["reported_inconsistencies"][0]
    assert all(
        any(
            "Z=2" in issue and "Z=4" in issue
            for issue in validation["reported_inconsistencies"]
        )
        for validation in (anderson, jamieson, dewaele)
    )
