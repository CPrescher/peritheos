"""Miozzi source transcription and explicitly unvalidated EOS candidates."""

import csv
import hashlib

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_eos_record, get_material_document
from scripts.reproduce_miozzi_2020_iron import (
    RECORDS,
    ROOT,
    ledger_outcome,
    observations,
    reproduce,
    source_pressure,
)


def test_complete_source_tables_and_missing_errors():
    doc = get_material_document("iron")
    source_datasets = [
        ds
        for ds in doc["datasets"]
        if ds["identifier"] in {"iron_miozzi_2020_he_pvt", "iron_miozzi_2020_mgo_pvt"}
    ]
    for ds, count in zip(source_datasets, (15, 116)):
        path = ROOT / "peritheos/data" / ds["resource"]["path"]
        assert hashlib.sha256(path.read_bytes()).hexdigest() == ds["resource"]["sha256"]
        with path.open(encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        assert len(rows) == count
        assert [int(r["source_row"]) for r in rows] == list(range(1, count + 1))
        assert ds["columns"][3]["role"] == "uncertainty"
        assert all(r["status"] == "unresolved" for r in ds["pressure_reductions"])
        if count == 116:
            assert [
                (r["source_row"], k) for r in rows for k, v in r.items() if v == ""
            ] == [
                ("1", "pressure_error_gpa"),
                ("27", "volume_error_a3"),
                ("76", "temperature_error_k"),
            ]
            assert rows[0]["mgo_volume_a3"] == "62.625"
            assert rows[-1]["volume_a3"] == "16.668"
            assert sum(r["pressure_error_gpa"] == "0" for r in rows) == 28
        else:
            assert rows[0]["volume_a3"] == "21.313"
            assert rows[-1]["pressure_gpa"] == "30.8"
            assert all(r["temperature_k"] == "298" for r in rows)
    data = observations()
    assert len(data["pressure_gpa"]) == 131
    assert np.sum(data["temperature_k"] <= 300) == 36
    assert max(data["temperature_k"]) == 3399
    assert max(data["pressure_gpa"]) == 127.9736


@pytest.mark.parametrize("identifier", RECORDS)
def test_published_parameters_independent_native_pressure_and_roundtrip(identifier):
    doc = get_material_document("iron")
    Draft202012Validator(eosmat_schema()).validate(doc)
    raw = next(r for r in doc["eos_records"] if r["identifier"] == identifier)
    assert raw["eos"]["parameters"] == dict(
        zip(("V0", "K0", "K0_prime"), RECORDS[identifier][:3])
    )
    assert raw["default"] is False
    assert raw["catalog_access"] == "explicit_selection"
    assert raw["scientific_validation"]["status"] == "not_reproduced"
    with pytest.raises(ValueError, match="not_reproduced"):
        Material.from_eosmat(doc, record_identifiers=[identifier])
    material = Material.from_eosmat(
        doc, record_identifiers=[identifier], require_primary_validation=False
    )
    exported = material.to_eosmat()
    assert (
        exported["eos_records"][0]["scientific_validation"]
        == raw["scientific_validation"]
    )
    source_datasets = [
        ds for ds in doc["datasets"] if ds["identifier"] in raw["fit_datasets"]
    ]
    for source, restored in zip(source_datasets, exported["datasets"]):
        for key in ("identifier", "resource", "columns", "uncertainty", "notes"):
            assert restored[key] == source[key]
        assert restored["used_by_eos_records"] == [identifier]
        assert restored["pressure_reductions"] == [
            r for r in source["pressure_reductions"] if r["eos_record"] == identifier
        ]
    record = get_eos_record(identifier)
    volumes = np.array([21.313, 19.269, 18.974, 16.668])
    temperatures = (
        np.array([300, 1759, 3399, 2386])
        if len(RECORDS[identifier]) == 5
        else np.full(4, 300)
    )
    independent = source_pressure(
        volumes, temperatures, RECORDS[identifier], vinet="vinet" in identifier
    )
    native = record.pressure(volumes, temperatures)
    np.testing.assert_allclose(native, independent, rtol=2e-11, atol=2e-10)
    np.testing.assert_allclose(record.volume(native, temperatures), volumes, rtol=1e-10)
    assert record.pressure(RECORDS[identifier][0], 300) == pytest.approx(0, abs=1e-10)


def test_source_reproduction_discrepancy_is_preserved():
    results = reproduce()
    assert results["iron_miozzi_2020_bm3"]["published_rmse_gpa"] == pytest.approx(
        3.7536167053
    )
    assert results["iron_miozzi_2020_bm3_mgd"]["published_rmse_gpa"] == pytest.approx(
        2.4988584007
    )
    assert results["iron_miozzi_2020_bm3_mgd"]["published_residual_range_gpa"][0] < -6.9
    raw = next(
        r
        for r in get_material_document("iron")["eos_records"]
        if r["identifier"] == "iron_miozzi_2020_bm3_mgd"
    )
    assert raw["fixed_parameters"] == []
    assert raw["thermal"]["parameters"]["n"] == 2
    communication = raw["scientific_validation"]["personal_communication"]
    assert communication["type"] == "personal_communication"
    assert communication["attribution"] == "Miozzi et al."
    assert communication["recorded_on"] == "2026-10-03"
    assert raw["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert raw["parameter_errors"]["V0"] == pytest.approx(0.06642156269)
    for r in get_material_document("iron")["eos_records"]:
        if r["identifier"] not in RECORDS:
            continue
        assert ledger_outcome(r)["status"] == "parity_not_achieved"
    for item in results.values():
        for fit in item["fits"].values():
            assert fit["solver_success"]
            assert len(fit["standard_errors"]) == len(fit["free_parameter_names"])
            assert np.all(np.array(fit["standard_errors"]) > 0)
    assert (
        results["iron_miozzi_2020_bm3"]["fits"]["effective_errors"]["row_count"] == 35
    )
