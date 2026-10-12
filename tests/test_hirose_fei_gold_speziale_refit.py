"""Recalculated pressure provenance and executable combined Au refit."""

import copy
import csv
import hashlib
import io
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    eosmat_schema,
    find_pressure_calibration_path,
    get_eos_record,
    get_eos_record_document,
    get_material,
    get_material_document,
)
from scripts.check_numerical_archive import check_csv
from scripts.refit_hirose_fei_gold_speziale import (
    CSV,
    DATA,
    DATASET_ID,
    MGO_ID,
    NAMES,
    RECORD_ID,
    REPORT,
    append_object,
    gold_pressure,
    ledger_outcome,
    reproduce,
)
from scripts.reproduce_fei_2007_gold import check_saved


@pytest.fixture(scope="module")
def result():
    return reproduce()


def test_exact_source_inputs_and_recalculated_pressure_coordinate(result):
    report, text = result
    rows = list(csv.DictReader(io.StringIO(text)))
    assert len(rows) == 47
    assert sum(r["fit_stage"] == "cold" for r in rows) == 9
    assert sum(r["source"] == "Fei2004" for r in rows) == 26
    assert (
        sum(r["source"] == "Hirose2008" and r["fit_stage"] == "thermal" for r in rows)
        == 12
    )
    originals = {}
    for source, name in (
        ("Hirose2008", "gold-hirose-2008-table1.csv"),
        ("Fei2004", "gold-fei-2004-table1.csv"),
    ):
        path = DATA / name
        assert (
            hashlib.sha256(path.read_bytes()).hexdigest()
            == report["source_hashes"][name]
        )
        with path.open() as stream:
            originals[source] = list(csv.DictReader(stream))
    for row in rows:
        original = originals[row["source"]][int(row["source_row"]) - 1]
        for key in (
            "volume_a3",
            "volume_uncertainty_a3",
            "mgo_volume_a3",
            "temperature_k",
        ):
            assert row[key] == original[key]
        assert row["reported_pressure_gpa"] == original["pressure_gpa"]
        assert (
            row["reported_pressure_width_gpa"] == original["pressure_gpa_uncertainty"]
        )
        expected = get_eos_record(MGO_ID).pressure(
            float(row["mgo_volume_a3"]), float(row["temperature_k"])
        )
        assert float(row["pressure_gpa"]) == pytest.approx(expected, abs=1e-9)
    assert min(float(r["reported_pressure_gpa"]) for r in rows) == 8.64
    assert any(r["volume_uncertainty_a3"] == "" for r in rows)
    assert report["fei_eos_work_paused"] is True


def test_registered_native_fallback_and_interchange_match_staged_fit(result):
    report, text = result
    record = get_eos_record(RECORD_ID)
    rows = list(csv.DictReader(io.StringIO(text)))
    v = np.array([float(r["volume_a3"]) for r in rows])
    t = np.array([float(r["temperature_k"]) for r in rows])
    expected = gold_pressure(
        v,
        t,
        report["thermal_stage"]["parameters"],
        report["cold_stage"]["parameters"][0],
    )
    np.testing.assert_allclose(record.pressure(v, t), expected, atol=1e-9)
    fallback = copy.copy(record.eos)
    del fallback._native
    np.testing.assert_allclose(fallback.pressure(v, t), expected, atol=1e-9)
    assert np.all(record.eos.bulk_modulus(v, t) > 0)
    np.testing.assert_allclose(record.eos.thermal_pressure(v, 300), 0, atol=1e-12)
    pressure = record.pressure([64, 51], [1473, 2200])
    np.testing.assert_allclose(
        record.volume(pressure, [1473, 2200]), [64, 51], atol=1e-8
    )
    selected = Material.from_eosmat(
        get_material_document("gold"), record_identifiers=[RECORD_ID]
    )
    restored = Material.from_eosmat(selected.to_eosmat()).get_eos_record(RECORD_ID)
    np.testing.assert_allclose(restored.pressure(v, t), expected, atol=1e-9)
    raw = get_eos_record_document(RECORD_ID)
    assert raw["record_kind"] == "refit"
    assert raw["default"] is False
    assert raw["parameter_covariance"] is None
    assert all(value is None for value in raw["thermal"]["parameter_errors"].values())
    assert report["thermal_stage"]["pressure_mae_gpa"] == pytest.approx(
        0.315847, abs=1e-6
    )
    assert raw["eos"]["parameters"]["K0_prime"] == pytest.approx(5.58115032)
    assert find_pressure_calibration_path(RECORD_ID, MGO_ID)
    assert get_eos_record("gold_fei_2007_vinet_2").is_default
    assert (
        get_eos_record_document("gold_hirose_2008_bm3_fit2")["thermal"]["parameters"][
            "beta1"
        ]
        == 1.03e-6
    )


def test_dataset_load_schema_archives_and_ledger_dispatch(result):
    report, table = result
    dataset = get_material("gold").get_dataset(DATASET_ID)
    pv = dataset.as_pressure_volume()
    assert len(pv.volume) == 47
    assert pv.pressure_uncertainty is None
    doc = get_material_document("gold")
    Draft202012Validator(eosmat_schema()).validate(doc)
    check_saved(json.loads(REPORT.read_text()), report)
    check_csv(CSV.read_text(), table, {"pressure_gpa"}, atol=1e-8)
    assert hashlib.sha256(CSV.read_bytes()).hexdigest() == report["derived_csv_sha256"]
    raw = get_eos_record_document(RECORD_ID)
    outcome = ledger_outcome(raw)
    assert outcome["observations"] == 47
    assert outcome["original_publication_reproduction_status"] == "not_reproduced"
    for name, value in zip(NAMES, report["thermal_stage"]["parameters"]):
        assert np.isclose(
            raw["thermal"]["parameters"][name], value, rtol=1e-4, atol=1e-16
        )


def test_targeted_registration_preserves_other_json_members_byte_for_byte():
    text = '{\n "eos_records": [\n  {"identifier":"keep", "value": 1.000}\n ]\n}\n'
    own = {"identifier": RECORD_ID, "value": 2}
    added = append_object(text, "eos_records", own)
    changed = append_object(added, "eos_records", {**own, "value": 3})
    assert '{"identifier":"keep", "value": 1.000}' in changed
    parsed = json.loads(changed)
    assert len(parsed["eos_records"]) == 2
    assert parsed["eos_records"][1]["value"] == 3
