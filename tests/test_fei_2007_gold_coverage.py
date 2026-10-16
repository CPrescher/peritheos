"""Fei Au coverage follows measured states, separate from EOS reference/outputs."""

import numpy as np
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_material_document
from scripts.reproduce_fei_2007_gold import observations, read_csv

RECORD = "gold_fei_2007_vinet_2"


def test_coverage_includes_cold_anchor_and_both_high_pressure_hot_states():
    document = get_material_document("gold")
    Draft202012Validator(eosmat_schema()).validate(document)
    raw = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    cold = read_csv("gold-dewaele-2004-table1-compression.csv")
    hot_rows, hot = observations()
    assert len(cold) == 37
    assert len(hot_rows) == 28
    cold_pressure = [float(r["ruby_pressure_revised_gpa"]) for r in cold]
    assert [min(cold_pressure), max(cold_pressure)] == [0.0, 93.6]
    assert hot["p"][-2:].tolist() == [106.5, 117.9]
    assert hot["t"][-2:].tolist() == [1340.0, 2330.0]
    measured_pressure = np.r_[cold_pressure, hot["p"]]
    measured_temperature = np.r_[np.full(len(cold), 298.0), hot["t"]]
    assert raw["experimental_pressure_range_gpa"] == [
        measured_pressure.min(),
        measured_pressure.max(),
    ]
    assert raw["experimental_temperature_range_k"] == [
        measured_temperature.min(),
        measured_temperature.max(),
    ]
    assert raw["pressure_range_status"] == "reported_exactly"
    assert raw["thermal"]["parameters"]["Tr"] == 300.0


def test_native_validity_retains_marginal_coverage_qualifications_on_export():
    material = Material.from_eosmat(
        get_material_document("gold"), record_identifiers=[RECORD]
    )
    record = material.eos_records[0]
    assert record.validity.pressure_gpa == (0.0, 117.9)
    assert record.validity.temperature_k == (298.0, 2330.0)
    notes = " ".join(record.validity.notes)
    assert "not a validated rectangular P-T domain" in notes
    assert "original author selection is not recovered" in notes
    assert "298 K" in notes and "300 K" in notes
    assert "pressure gap" in notes
    assert "calculated Figure 1 vertices" in notes
    assert "no explicit recommended numerical extrapolation range" in notes
    exported = material.to_eosmat()["eos_records"][0]
    assert exported["experimental_pressure_range_gpa"] == [0.0, 117.9]
    assert exported["experimental_temperature_range_k"] == [298.0, 2330.0]
    assert exported["validity"]["notes"] == list(record.validity.notes)
    assert exported["default"] is True
