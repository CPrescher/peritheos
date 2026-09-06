import importlib.util
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_yang_reproduction_and_schema():
    path = ROOT / "scripts/reproduce_yang_2015_ferropericlase.py"
    spec = importlib.util.spec_from_file_location("reproduce_yang_2015", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert module.reproduce()["record"] == "mg092fe008o_yang_2015_hs_bm3_reference"
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    assert (
        list(
            Draft202012Validator(schema).iter_errors(
                get_material_document("mg092fe008o")
            )
        )
        == []
    )


def test_yang_audit_enumerates_all_candidates():
    text = (
        ROOT / "docs/literature-reproductions/yang-2015-ferropericlase.md"
    ).read_text(encoding="utf-8")
    assert all(
        identifier in text
        for identifier in {
            "litcurate_681529622576473b",
            "litcurate_0cb0d7a43e4cda51",
        }
    )
