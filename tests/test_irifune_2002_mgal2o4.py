import csv
import hashlib
import json
from pathlib import Path

import pytest

from peritheos import Material, get_material_document
from peritheos.eos.rt import BM3
from scripts.reproduce_irifune_2002_mgal2o4 import reproduce

ROOT = Path(__file__).parents[1]
RECORD = "mgal2o4_cafe2o4_irifune_2002_bm2_2"
DATASET = "mgal2o4_cafe2o4_irifune_2002_text_pv"
DOI = "10.1007/s00269-002-0275-1"


def test_irifune_record_preserves_constrained_source_and_existing_default():
    document = get_material_document("mgal2o4_cafe2o4")
    assert document["formula_units_per_cell"] == 4
    assert document["space_group"] == "Pbnm"
    assert document["eos_records"][0]["default"] is True
    material = Material.from_eosmat(document)
    stored = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert stored["reference"]["doi"] == DOI
    assert stored["default"] is False
    assert stored["eos"]["parameters"] == {"V0": 240.6, "K0": 213.0}
    assert stored["parameter_errors"] == {"V0": 0.3, "K0": 3.0}
    assert stored["fixed_parameters"] == ["V0"]
    assert stored["parameter_error_confidence"] is None
    assert stored["parameter_covariance"] is None
    assert stored["temperature_ref"] == pytest.approx(300.15)
    assert "thermal" not in stored
    roundtrip = material.to_eosmat()["eos_records"][1]
    assert {key: roundtrip[key] for key in stored} == stored
    executable = material.get_eos_record(RECORD)
    assert executable.pressure(240.6) == pytest.approx(0.0)
    bm3 = BM3(V0=240.6, K0=213.0, K0_prime=4.0)
    for volume in (216.1, 213.3, 213.1):
        assert executable.pressure(volume) == pytest.approx(bm3.pressure(volume))
    for pressure in (0.0, 27.4, 32.1, 34.7):
        volume = executable.volume(pressure, check_validity=True)
        assert executable.pressure(volume) == pytest.approx(pressure, abs=1e-8)
    with pytest.raises(ValueError):
        executable.volume(35.0, check_validity=True)


def test_irifune_complete_transcription_and_sparse_fit_selection():
    document = get_material_document("mgal2o4_cafe2o4")
    dataset = next(d for d in document["datasets"] if d["identifier"] == DATASET)
    assert dataset["used_by_eos_records"] == [RECORD]
    path = ROOT / "peritheos/data" / dataset["resource"]["path"]
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest() == dataset["resource"]["sha256"]
    )
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 6
    assert [float(r["volume_a3"]) for r in rows] == [
        216.1,
        213.3,
        213.1,
        240.7,
        240.0,
        241.0,
    ]
    assert [float(r["volume_error_a3"]) for r in rows] == [0.3, 0.4, 0.7, 0.3, 0.2, 0.3]
    assert [float(r["pressure_gpa"]) for r in rows] == [27.4, 32.1, 34.7, 0, 0, 0]
    assert [r["fit_included"] for r in rows] == ["1"] * 3 + ["0"] * 3


def test_irifune_independent_reproduction_and_ledger():
    result = reproduce()
    assert result["compressed_observations"] == 3
    assert result["free_parameters"] == ["K0"]
    assert result["fixed_v0_pressure_refit_k0_gpa"] == pytest.approx(213.86541623)
    assert result["fixed_v0_normalized_volume_refit_k0_gpa"] == pytest.approx(
        213.65992611
    )
    assert result["published_pressure_rmse_gpa"] == pytest.approx(1.188482859)
    executable = Material.from_eosmat(
        get_material_document("mgal2o4_cafe2o4")
    ).get_eos_record(RECORD)
    for volume, expected in zip(
        (216.1, 213.3, 213.1), result["published_pressures_gpa"]
    ):
        assert executable.pressure(volume) == pytest.approx(expected, abs=1e-10)
    # 2 GPa envelopes the disclosed scatter, not an invented measurement error.
    for volume, observed in zip((216.1, 213.3, 213.1), (27.4, 32.1, 34.7)):
        assert executable.pressure(volume) == pytest.approx(observed, abs=2.0)
    ledger = json.loads(
        (ROOT / "docs/data/primary-eos-refits.json").read_text(encoding="utf-8")
    )
    item = next(r for r in ledger["records"] if r["record_identifier"] == RECORD)
    assert item["status"] == "parity"
    assert item["observations"] == 3
    assert item["fixed_parameters"] == ["V0"]
    assert "diagnostic" in item["qualification"]
    candidates = json.loads(
        (ROOT / "docs/data/litcurate-eos-candidates.json").read_text(encoding="utf-8")
    )
    assert DOI not in json.dumps(candidates)
