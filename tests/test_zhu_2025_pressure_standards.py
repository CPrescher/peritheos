import importlib.util
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_zhu_reproduction_and_schemas():
    path = ROOT / "scripts/reproduce_zhu_2025_pressure_standards.py"
    spec = importlib.util.spec_from_file_location("reproduce_zhu_2025", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert len(module.reproduce()["records"]) == 3
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    validator = Draft202012Validator(schema)
    for material in ("platinum", "gold", "mgo"):
        assert list(validator.iter_errors(get_material_document(material))) == []


def test_zhu_audit_enumerates_all_candidates():
    text = (
        ROOT / "docs/literature-reproductions/zhu-2025-pressure-standards.md"
    ).read_text(encoding="utf-8")
    assert "litcurate_b1a3d7a39cb0fae8" in text
    assert "litcurate_fdde71245d253100" in text
