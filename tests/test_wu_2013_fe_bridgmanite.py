import importlib.util
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_wu_reproduction_and_schemas():
    path = ROOT / "scripts/reproduce_wu_2013_fe_bridgmanite.py"
    spec = importlib.util.spec_from_file_location("reproduce_wu_2013", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert len(module.reproduce()["records"]) == 2
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    validator = Draft202012Validator(schema)
    for material in ("bridgmanite", "mg075fe025sio3_bridgmanite"):
        assert list(validator.iter_errors(get_material_document(material))) == []


def test_wu_audit_enumerates_all_candidates():
    text = (ROOT / "docs/literature-reproductions/wu-2013-fe-bridgmanite.md").read_text(
        encoding="utf-8"
    )
    assert all(
        identifier in text
        for identifier in {
            "litcurate_1fd4855251635f0a",
            "litcurate_eb3b372892561702",
            "litcurate_2edf2b54fb35f7a2",
            "litcurate_8d4961bd7435eb6b",
        }
    )
