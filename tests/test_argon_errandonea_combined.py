"""Combined-source diagnostics must preserve temperature and source anomalies."""

import csv
import hashlib
import json

import numpy as np
import pytest

from scripts.fit_argon_errandonea_combined import (
    INPUT,
    calculate,
    inverse,
    observations,
)
from scripts.reproduce_argon_errandonea_2006 import bm3


def test_registered_refit_preserves_defaults_and_selected_observations():
    from jsonschema import Draft202012Validator

    from peritheos import Material, get_material_document
    from scripts.fit_argon_errandonea_combined import PREFERRED, REPORT
    from scripts.register_argon_combined_refit import DATASET_ID, RECORD_ID
    from scripts.reproduce_argon_errandonea_2006 import ROOT

    doc = get_material_document("argon_fcc")
    Draft202012Validator(
        json.loads((ROOT / "peritheos/data/eosmat-v3.schema.json").read_text())
    ).validate(doc)
    record = next(r for r in doc["eos_records"] if r["identifier"] == RECORD_ID)
    assert record["record_kind"] == "refit"
    assert record["default"] is False
    assert [r["identifier"] for r in doc["eos_records"] if r.get("default")] == [
        "argon_fcc_dewaele_2021_vinet_mgd"
    ]
    assert (
        record["eos"]["parameters"]
        == json.loads(REPORT.read_text())["fits"][PREFERRED]["pressure"]["parameters"]
    )
    dataset = next(d for d in doc["datasets"] if d["identifier"] == DATASET_ID)
    path = ROOT / "peritheos/data" / dataset["resource"]["path"]
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest() == dataset["resource"]["sha256"]
    )
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 49
    assert sum(r["source"] == "ross_1986" for r in rows) == 41
    assert not any(r["source"] == "ross_1986" and r["source_row"] == "17" for r in rows)
    model = Material.from_eosmat(doc).get_eos_record(RECORD_ID)
    volumes = np.array([float(r["volume_a3"]) for r in rows])
    parameters = record["eos"]["parameters"]
    assert model.pressure(volumes, 300) == pytest.approx(
        bm3(volumes, parameters["V0"], parameters["K0"], parameters["K0_prime"]),
        abs=1e-10,
    )
    restored = Material.from_eosmat(
        Material.from_eosmat(doc).to_eosmat()
    ).get_eos_record(RECORD_ID)
    assert restored.pressure(volumes, 300) == pytest.approx(
        model.pressure(volumes, 300), abs=1e-10
    )


def test_combined_source_selection_and_units():
    source = json.loads(INPUT.read_text())
    ross = source["ross"]["rows"]
    assert len(ross) == 42
    for row in ross:
        assert float(row["pressure_gpa"]) == float(row["pressure_kbar"]) / 10
        assert float(row["volume_a3"]) == pytest.approx(
            float(row["molar_volume_cm3_mol"]) * 4e24 / 6.02214076e23, abs=1e-8
        )
    assert float(ross[16]["pressure_gpa"]) == 24.7
    rows = observations()
    assert len(rows) == 51
    assert {r["temperature_k"] for r in rows if r["source"] == "ross_1986"} == {298}
    assert rows[-1]["temperature_k"] is None
    assert rows[-1]["source"] == "anderson_swenson_1975"
    assert sum(r["included_in_preferred_fit"] for r in rows) == 49
    suspect = next(
        r for r in rows if r["source"] == "ross_1986" and r["source_row"] == 17
    )
    assert suspect["pressure_gpa"] == 24.7
    assert not suspect["included_in_preferred_fit"]
    assert "Suspected printed pressure error" in suspect["exclusion_reason"]


def test_fit_sensitivity_and_inversion():
    pressure = np.array([0, 1.6, 40, 114.0])
    assert bm3(inverse(pressure, [143, 6.5, 5.1])) == pytest.approx(pressure, abs=1e-9)
    report = calculate()
    fits = report["fits"]
    assert report["preferred_fit"]["subset"] == "room_temperature_without_ross_row17"
    assert fits["room_temperature_union"]["pressure"]["count"] == 50
    assert fits["room_temperature_without_ross_row17"]["pressure"]["count"] == 49
    assert fits["all_including_cryogenic_marker"]["pressure"]["count"] == 51
    assert fits["errandonea_only"]["pressure"]["active_bound_mask"][0] == 1
    rt = fits["room_temperature_union"]
    assert rt["pressure"]["parameters_inside_published_error_intervals"]
    assert rt["pressure"]["parameters"]["V0"] == pytest.approx(141.6932, abs=0.002)
    assert rt["volume"]["parameters"]["V0"] > 180
    assert all(
        result["solver_success"] for group in fits.values() for result in group.values()
    )
