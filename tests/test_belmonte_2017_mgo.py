import hashlib
import importlib.util
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]
IDENTIFIER = "mgo_belmonte_2017_b3lyp_static_bm3"


def _load_reproduction():
    path = ROOT / "scripts/reproduce_belmonte_2017_mgo.py"
    spec = importlib.util.spec_from_file_location("reproduce_belmonte_2017_mgo", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_belmonte_record_and_table4_reproduction():
    document = get_material_document("mgo")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == IDENTIFIER
    )
    assert record["reference"]["doi"] == "10.3390/min7100183"
    assert record["eos"] == {
        "type": "BM3",
        "model": "birch_murnaghan_3",
        "parameters": {"V0": 75.762, "K0": 167.01, "K0_prime": 3.95},
    }
    assert record["volume_basis"] == {
        "kind": "formula_units",
        "formula_units": 4.0,
        "molar_mass_g_mol": 40.304,
    }
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    assert eos.pressure(75.762) == pytest.approx(0.0, abs=1e-12)

    result = _load_reproduction().reproduce()
    assert result["observations"] == 8
    assert result["max_abs_residual_gpa"] < 0.005


def test_belmonte_dataset_hash_and_schema():
    path = ROOT / "peritheos/data/datasets/mgo-belmonte-2017-table4-300k.csv"
    assert (
        hashlib.sha256(path.read_bytes()).hexdigest()
        == "85ba8923a840aec8eacf084807dc62a5afacfddbf0c04964569d01e408c265b7"
    )
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    assert (
        list(Draft202012Validator(schema).iter_errors(get_material_document("mgo")))
        == []
    )


def test_belmonte_audit_enumerates_all_seven_candidates():
    text = (ROOT / "docs/literature-reproductions/belmonte-2017-mgo.md").read_text(
        encoding="utf-8"
    )
    identifiers = {
        "litcurate_4800133db6bddd4b",
        "litcurate_aa2c82a02036dbfa",
        "litcurate_516d0cf64b5f460d",
        "litcurate_3f43b1de0ed65a6f",
        "litcurate_fd66dc4be752b6ef",
        "litcurate_841c8a89f6b4ed48",
        "litcurate_a29f90dd6872f077",
    }
    assert all(identifier in text for identifier in identifiers)
