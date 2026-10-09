"""Primary-table, cell-basis and pressure-coordinate regression for KCl variants."""

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_kcl_variants import A3_PER_MOLAR, IDS, reproduce, rows

ROOT = Path(__file__).resolve().parents[1]


def test_walker_b1_retains_all_source_rows_and_real_sample_selection():
    data = rows("kcl-walker-2002-table1-pvt.csv")
    assert len(data) == 30
    assert sum(r["row_kind"] == "calibrant_spot_check" for r in data) == 6
    assert sum(r["included_in_fit"] == "1" for r in data) == 23
    anchor = data[0]
    assert anchor["row_kind"] == "derived_reference"
    assert anchor["pressure_kbar"] == anchor["b1_kcl_cell_volume_esd_a3"] == ""
    assert anchor["b1_kcl_cell_volume_a3"] == "249.53"
    last = data[-1]
    assert last["spectrum"] == "r57801"
    assert last["pressure_kbar"] == "2.10"
    assert last["included_in_fit"] == "0"
    record = get_eos_record(IDS[0])
    assert record.reference_volume == pytest.approx(249.080860076)
    assert record.reference_temperature == 296.15
    assert record.pressure(record.reference_volume, 296.15) == pytest.approx(
        0, abs=1e-12
    )
    # Abstract's directly reported product; this is not rounded K0*alpha0.
    assert record.pressure(
        record.reference_volume, 873.15, check_validity=False
    ) == pytest.approx(0.00195 * 577)
    # Independent off-reference measured Table 1 state, with propagated
    # coordinate/product-error tolerance (not exact source-fit parity).
    assert record.pressure(242.28, 873.15, check_validity=False) == pytest.approx(
        1.726, abs=0.08
    )


def test_official_tateno_pairing_has_correct_marker_and_sample_volume():
    raw = rows("kcl-tateno-2019-official-table-s1.csv")
    rounded = rows("kcl-tateno-2019-table-s1-pvt.csv")
    assert len(raw) == len(rounded) == 39
    for official, old in zip(raw, rounded):
        assert float(old["pressure_gpa"]) == pytest.approx(
            float(official["pressure_gpa"]), abs=0.051
        )
        assert float(old["kcl_unit_cell_volume_a3"]) == pytest.approx(
            float(official["kcl_unit_cell_volume_a3"]), abs=0.000051
        )
        assert float(old["platinum_unit_cell_volume_a3"]) == pytest.approx(
            float(official["platinum_unit_cell_volume_a3"]), abs=0.0051
        )
        # The article rounds a to four decimals and V to two decimals.
        lattice = float(old["platinum_lattice_a_angstrom"])
        assert lattice**3 == pytest.approx(
            float(old["platinum_unit_cell_volume_a3"]),
            abs=0.0051 + 3 * lattice**2 * 0.000051,
        )
    six = next(r for r in raw if r["source_excel_row"] == "24")
    assert float(six["platinum_unit_cell_volume_a3"]) == pytest.approx(59.1561974927066)
    assert float(six["kcl_unit_cell_volume_a3"]) == pytest.approx(45.04627509139811)
    corrections = json.loads(
        (ROOT / "peritheos/data/datasets/kcl_variant_sources/manifest.json").read_text()
    )["tateno_pt_pairing_corrections"]
    assert len(corrections) == 110
    assert all(r["column"].startswith("platinum_") for r in corrections)


def test_holmes_pressure_coordinate_is_separate_and_explicitly_qualified():
    result = reproduce()
    derived = result["holmes_pressure_coordinate"]
    assert len(derived["rows"]) == 39
    assert any(r["outside_holmes_thermal_approximation"] for r in derived["rows"])
    raw = rows("kcl-tateno-2019-official-table-s1.csv")
    # Direct native calibrant calculation confirms independent Eq. 11/12 replay.
    pt = get_eos_record("platinum_holmes_1989_vinet_1")
    for source, calculated in zip(raw, derived["rows"]):
        assert pt.pressure(
            float(source["platinum_unit_cell_volume_a3"]),
            float(source["temperature_k"]),
            check_validity=False,
        ) == pytest.approx(calculated["pressure_gpa"], abs=1e-10)
    assert (
        max(
            abs(r["pressure_gpa"] - float(s["pressure_gpa"]))
            for r, s in zip(derived["rows"], raw)
        )
        > 4
    )
    assert derived["calibration_parameter_errors"] is None
    for identifier in IDS[1:3]:
        assert "conditional Holmes" in result[identifier]["pressure_coordinate"]
        assert "T<2000 K" in result[identifier]["qualification"]
    assert (
        result[IDS[3]]["pressure_coordinate"] == "author-reported Sokolova Pt pressures"
    )


def test_complete_chidester_vinet_fit_and_independent_equation_reproduction():
    result = reproduce()
    c = result[IDS[4]]
    assert c["observations"] == 278
    assert (c["room_temperature_rows"], c["effective_temperature_rows"]) == (123, 155)
    doc = get_material_document("kcl")
    record = next(r for r in doc["eos_records"] if r["identifier"] == IDS[4])
    assert record["eos"]["parameters"]["V0"] == pytest.approx(34.3 * A3_PER_MOLAR)
    assert record["parameter_errors"]["V0"] == pytest.approx(0.5 * A3_PER_MOLAR)
    assert "q" not in record["thermal"]["fixed_parameters"]
    assert record["thermal"]["parameter_errors"]["q"] == 0.1
    cold = next(
        item
        for item in record["source_lineage"]
        if item["role"] == "room_temperature_observations"
    )
    assert cold["reference_calibration_record"] == "ruby_dorogokupets_oganov_2007"
    assert record["pressure_calibration"]["methods"][0]["reference_eos_record"] == (
        "platinum_dorogokupets_oganov_2007_vinet_4"
    )
    for name, fitted in c["fitted_parameters"].items():
        block = record if name in record["eos"]["parameters"] else record["thermal"]
        parameters = (
            record["eos"]["parameters"] if block is record else block["parameters"]
        )
        assert abs(fitted - parameters[name]) < block["parameter_errors"][name]
    for identifier in IDS:
        assert result[identifier]["native_max_difference_gpa"] < 1e-9
        assert result[identifier].get("quadrature_max_difference_gpa", 0) < 1e-10
    assert c["refit"]["rmse_gpa"] < 1.3
    assert result[IDS[0]]["fitted_parameters"]["K0"] == pytest.approx(17.68342185)
    assert (
        "strict coefficient/uncertainty parity is not established"
        in result[IDS[0]]["qualification"]
    )


@pytest.mark.parametrize("identifier", IDS)
def test_variants_remain_nondefault_and_execute_invert_and_round_trip(identifier):
    material_id = "kcl_b1" if identifier == IDS[0] else "kcl"
    document = get_material_document(material_id)
    record_doc = next(
        r for r in document["eos_records"] if r["identifier"] == identifier
    )
    assert record_doc["default"] is False
    material = Material.from_eosmat(document, record_identifiers=[identifier])
    record = material.eos_records[0]
    temp = 600 if material_id == "kcl_b1" else 1500
    volumes = (
        np.array([0.98, 0.96]) * record.reference_volume
        if material_id == "kcl_b1"
        else np.array([35, 30])
    )
    pressures = record.pressure(volumes, temp, check_validity=False)
    assert np.all(np.isfinite(pressures))
    assert record.volume(pressures, temp, check_validity=False) == pytest.approx(
        volumes, rel=1e-9
    )
    round_trip = Material.from_eosmat(material.to_eosmat()).eos_records[0]
    assert round_trip.pressure(volumes, temp, check_validity=False) == pytest.approx(
        pressures
    )
    assert round_trip.pressure_calibration == record.pressure_calibration


def test_defaults_and_every_variant_dataset_checksum():
    for material in ("kcl", "kcl_b1"):
        doc = get_material_document(material)
        assert get_eos_record(doc["eos_records"][0]["identifier"]).identifier not in IDS
        defaults = [r["identifier"] for r in doc["eos_records"] if r.get("default")]
        assert defaults == (
            ["kcl_b2_dewaele_2012_vinet_3"] if material == "kcl" else []
        )
        for dataset in doc["datasets"]:
            resource = dataset["resource"]
            payload = (ROOT / "peritheos/data" / resource["path"]).read_bytes()
            assert hashlib.sha256(payload).hexdigest() == resource["sha256"]
            parsed = list(csv.DictReader(payload.decode().splitlines()))
            assert parsed

    archive = ROOT / "peritheos/data/datasets/kcl_variant_sources"
    manifest = json.loads((archive / "manifest.json").read_text())
    for source in manifest["sources"]:
        if "filename" in source:
            assert (
                hashlib.sha256((archive / source["filename"]).read_bytes()).hexdigest()
                == source["sha256"]
            )
