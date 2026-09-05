import csv
import hashlib
import math
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_criniti_2021_bridgmanite import (
    cell_volume_from_density,
)

ROOT = Path(__file__).parents[1]
DOI = "10.1029/2020JB020967"
RECORD = "bridgmanite_criniti_2021_absolute_bm3"
DATASET = "bridgmanite_criniti_2021_table1_density"
RESOURCE = "bridgmanite-criniti-2021-table1-density.csv"
CHECKSUM = "f989ce7a1a715cedd17b7a9c013c8a4144b4b7ffd64f53c63bd108cf725ab9a8"


def _rows():
    path = ROOT / "peritheos" / "data" / "datasets" / RESOURCE
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def test_criniti_absolute_pressure_record_is_volumetric_and_executable():
    document = get_material_document("bridgmanite")
    stored = next(row for row in document["eos_records"] if row["identifier"] == RECORD)
    assert stored["reference"]["doi"].lower() == DOI.lower()
    assert stored["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 162.453274, "K0": 254.5, "K0_prime": 3.73},
    }
    assert stored["fixed_parameters"] == ["V0"]
    assert stored["pressure_calibration"]["methods"][0]["kind"] == "self_consistent"
    assert not any(
        row["identifier"] == "bridgmanite_criniti_2021_jacobsen_ruby_bm3"
        for row in document["eos_records"]
    )

    material = Material.from_eosmat(document, record_identifiers=[RECORD])
    record = material.get_eos_record(RECORD)
    assert float(record.pressure(np.array([162.453274]))[0]) == pytest.approx(0.0)


def test_criniti_density_data_and_conversion_are_bundled_exactly():
    document = get_material_document("bridgmanite")
    dataset = next(row for row in document["datasets"] if row["identifier"] == DATASET)
    assert dataset["resource"]["sha256"] == CHECKSUM
    path = ROOT / "peritheos" / "data" / "datasets" / RESOURCE
    assert hashlib.sha256(path.read_bytes()).hexdigest() == CHECKSUM
    rows = _rows()
    assert len(rows) == 14
    density = np.array([float(row["density_g_cm3"]) for row in rows])
    volume = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    assert cell_volume_from_density(density) == pytest.approx(volume, abs=5.0e-7)


def test_criniti_published_absolute_curve_matches_rounded_table_series():
    document = get_material_document("bridgmanite")
    material = Material.from_eosmat(document, record_identifiers=[RECORD])
    record = material.get_eos_record(RECORD)
    rows = _rows()
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    volume = np.array([float(row["volume_a3_conventional_cell"]) for row in rows])
    residual = np.asarray(record.pressure(volume)) - pressure
    assert math.sqrt(float(np.mean(residual**2))) == pytest.approx(
        0.0658748178, abs=5.0e-10
    )
    assert float(np.max(np.abs(residual))) == pytest.approx(0.2434557113, abs=5.0e-10)


def test_criniti_audit_disposes_all_five_litcurate_candidates():
    audit = (
        ROOT / "docs" / "literature-reproductions" / "criniti-2021-bridgmanite.md"
    ).read_text(encoding="utf-8")
    identifiers = {
        "litcurate_7c709a358383b534",
        "litcurate_1473f5b46c728a08",
        "litcurate_9bb38c0a21fcf50b",
        "litcurate_19a5bbcd8ab9c0e6",
        "litcurate_9c32311783c63f9c",
    }
    assert all(audit.count(identifier) == 1 for identifier in identifiers)
