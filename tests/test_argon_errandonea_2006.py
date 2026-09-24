"""Guard against phase/volume conflation and false reproduction claims."""

import csv
import hashlib
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, get_material_document
from scripts.reproduce_argon_errandonea_2006 import RECORD, ROOT, bm3, reproduce


def test_published_equation_native_and_volume_basis():
    doc = get_material_document("argon_fcc")
    schema = json.loads((ROOT / "peritheos/data/eosmat-v3.schema.json").read_text())
    Draft202012Validator(schema).validate(doc)
    model = Material.from_eosmat(doc).get_eos_record(RECORD)
    volumes = np.array([143.0, 90.0, 60.0, 49.0])
    assert model.pressure(volumes, 300) == pytest.approx(bm3(volumes), abs=1e-11)
    assert bm3(volumes / 4, v0=35.75) == pytest.approx(bm3(volumes), abs=1e-11)
    assert model.volume([34.9, 55, 114], 300) == pytest.approx(
        [67.04839488026091, 59.57824676453456, 48.83429160772853], abs=1e-8
    )
    restored = Material.from_eosmat(Material.from_eosmat(doc).to_eosmat())
    assert restored.get_eos_record(RECORD).pressure(volumes, 300) == pytest.approx(
        bm3(volumes), abs=1e-11
    )


def test_observation_provenance_and_no_hcp_fit_invention():
    doc = get_material_document("argon_fcc")
    record = next(r for r in doc["eos_records"] if r["identifier"] == RECORD)
    assert record["eos"]["parameters"] == {"V0": 143, "K0": 6.5, "K0_prime": 5.1}
    assert record["scientific_validation"]["reproduction_status"] == "not_reproduced"
    ds = next(d for d in doc["datasets"] if d["identifier"] in record["fit_datasets"])
    path = ROOT / "peritheos/data" / ds["resource"]["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ds["resource"]["sha256"]
    with path.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 8
    assert {r["phase"] for r in rows} == {"fcc"}
    audit = reproduce()
    assert audit["native_max_difference_gpa"] < 1e-11
    assert audit["published_partial_data_rmse_gpa"] < 3
    assert audit["diagnostic_equal_pressure_weight_refit"]["solver_success"]
    hcp = get_material_document("argon_hcp")
    assert (
        hcp["eos_records"][0]["scientific_validation"]["reproduction_status"]
        == "not_reproduced"
    )
    assert not any("errandonea" in r["identifier"] for r in hcp["eos_records"])
