import importlib.util
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_deng_record_reproduction_and_schema():
    path = ROOT / "scripts/reproduce_deng_2006_bridgmanite.py"
    spec = importlib.util.spec_from_file_location("reproduce_deng_2006", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert module.reproduce()["pressure_at_0.8_v0_gpa"] > 80.0
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    assert (
        list(
            Draft202012Validator(schema).iter_errors(
                get_material_document("bridgmanite")
            )
        )
        == []
    )


def test_deng_audit_enumerates_all_candidates():
    text = (ROOT / "docs/literature-reproductions/deng-2006-bridgmanite.md").read_text(
        encoding="utf-8"
    )
    identifiers = {
        "litcurate_e35eb6e4756f8473",
        "litcurate_d61fdcb2a93f7b6f",
        "litcurate_1402a7cf80b58a03",
        "litcurate_dbecb6f0878cef2b",
        "litcurate_df4f6d2ea24654be",
        "litcurate_fc6bb3f5f0b1a72c",
        "litcurate_d97945919d0deabc",
    }
    assert all(identifier in text for identifier in identifiers)
