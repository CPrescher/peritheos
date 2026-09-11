import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_ono_2006_mgal2o4 import DATA, bm2_pressure, reproduce

ROOT = Path(__file__).resolve().parents[1]
RECORD = "mgal2o4_cati2o4_ono_2006_bm2_2"
DOI = "10.1007/s00269-006-0068-z"


def test_ono_published_record_identity_uncertainties_and_calibration():
    document = get_material_document("mgal2o4_cati2o4")
    record = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert document["formula"] == "MgAl2O4"
    assert document["space_group"] == "Cmcm"
    assert document["formula_units_per_cell"] == 4
    assert document["source"]["structure_reference"]["year"] == 1998
    assert len(document["eos_records"]) == 2
    assert record["default"] is False
    assert record["reference"]["doi"] == DOI
    assert record["eos"] == {
        "type": "BM2",
        "model": "birch_murnaghan_2",
        "parameters": {"V0": 238.9, "K0": 219.0},
    }
    assert record["parameter_errors"] == {"V0": 0.9, "K0": 6.0}
    assert record["parameter_error_confidence"] is None
    assert record["parameter_covariance"] is None
    assert record["temperature_ref"] == 300.0
    assert record["experimental_pressure_range_gpa"] == [0.0, 91.2]
    assert record["experimental_temperature_range_k"] == [300.0, 300.0]
    assert "thermal" not in record
    calibration = record["pressure_calibration"]
    assert calibration["methods"][0]["reference_eos_record"] == (
        "platinum_holmes_1989_vinet_1"
    )
    assert calibration["recalculation"]["status"] == "missing_calibrant_observations"
    assert record["scientific_validation"]["status"] == "primary_source_validated"
    restored = Material.from_eosmat(document).to_eosmat()
    restored_record = next(
        r for r in restored["eos_records"] if r["identifier"] == RECORD
    )
    for field in (
        "parameter_provenance",
        "scientific_validation",
        "pressure_calibration",
    ):
        assert restored_record[field] == record[field]


def test_ono_table3_transcription_and_zero_rounded_errors():
    document = get_material_document("mgal2o4_cati2o4")
    dataset = document["datasets"][-1]
    assert dataset["used_by_eos_records"] == [RECORD]
    assert (
        hashlib.sha256(DATA.read_bytes()).hexdigest() == dataset["resource"]["sha256"]
    )
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 14
    assert [float(r["pressure_gpa"]) for r in rows] == [
        0,
        0,
        42.7,
        42.8,
        47.4,
        55.3,
        55.8,
        60.2,
        66.1,
        66.5,
        67.7,
        75.7,
        84.5,
        91.2,
    ]
    assert [float(r["volume_a3"]) for r in rows] == [
        240.7,
        239.3,
        206.0,
        205.7,
        203.1,
        201.0,
        199.6,
        197.3,
        194.6,
        194.7,
        195.2,
        191.8,
        188.3,
        186.4,
    ]
    assert {float(r["temperature_k"]) for r in rows} == {300.0}
    assert float(rows[2]["b_angstrom_error"]) == 0.037
    assert float(rows[4]["c_angstrom"]) == 8.979
    assert [i for i, r in enumerate(rows) if float(r["volume_a3_error"]) == 0] == [
        10,
        13,
    ]
    assert all(c["role"] == "uncertainty" for c in dataset["columns"] if "of" in c)
    assert not any(c["name"] == "pressure_error" for c in dataset["columns"])


def test_ono_independent_table_reproduction_and_public_inversion():
    diagnostic = reproduce()
    assert diagnostic["all_rows"]["V0"] == pytest.approx(238.9, abs=0.9)
    assert diagnostic["all_rows"]["K0"] == pytest.approx(219.0, abs=6.0)
    assert diagnostic["published_rmse_gpa"] == pytest.approx(1.32624757635)
    assert diagnostic["high_pressure_only"]["K0"] > 240
    record = Material.from_eosmat(
        get_material_document("mgal2o4_cati2o4")
    ).get_eos_record(RECORD)
    volumes = np.genfromtxt(DATA, delimiter=",", names=True)["volume_a3"]
    assert record.pressure(volumes) == pytest.approx(
        bm2_pressure(volumes, 238.9, 219.0)
    )
    assert record.pressure(238.9) == pytest.approx(0.0, abs=1e-10)
    assert record.pressure(186.4) == pytest.approx(89.3705180868)
    for pressure in (0.0, 42.7, 60.0, 91.2):
        assert record.pressure(record.volume(pressure)) == pytest.approx(
            pressure, abs=1e-8
        )


def test_ono_audit_ledgers_and_direct_candidate_disposition():
    audit = json.loads(
        (ROOT / "peritheos/data/primary-source-audit.json").read_text(encoding="utf-8")
    )
    assert len([r for r in audit["records"] if r["record"] == RECORD]) == 1
    refits = json.loads(
        (ROOT / "docs/data/primary-eos-refits.json").read_text(encoding="utf-8")
    )
    refit = next(r for r in refits["records"] if r["record_identifier"] == RECORD)
    assert refit["status"] == "parity"
    assert refit["observations"] == 14
    report = (
        ROOT / "docs/literature-reproductions/ono-2006-mgal2o4-cati2o4.md"
    ).read_text(encoding="utf-8")
    assert DOI in report
    assert "Citation traces only" in report
    assert "No EOS" in report
    assert RECORD in (ROOT / "docs/material-eos-candidates.md").read_text(
        encoding="utf-8"
    )
