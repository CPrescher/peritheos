import importlib.util
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_kiefer_reproduction_and_schema():
    path = ROOT / "scripts/reproduce_kiefer_2002_fe_bridgmanite.py"
    spec = importlib.util.spec_from_file_location("reproduce_kiefer_2002", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert set(result["records"]) == {
        "bridgmanite_kiefer_2002_gga_bm3",
        "mg075fe025sio3_bridgmanite_kiefer_2002_gga_bm3",
    }
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    validator = Draft202012Validator(schema)
    for material in ("bridgmanite", "mg075fe025sio3_bridgmanite"):
        assert list(validator.iter_errors(get_material_document(material))) == []


def test_kiefer_audit_enumerates_all_candidates():
    text = (
        ROOT / "docs/literature-reproductions/kiefer-2002-fe-bridgmanite.md"
    ).read_text(encoding="utf-8")
    identifiers = {
        "litcurate_6fc87b3bb31010b7",
        "litcurate_2add53b06987f13d",
        "litcurate_f948de0084755234",
        "litcurate_ff5b474ba6734f85",
        "litcurate_8774ba014956f30e",
        "litcurate_75b3b031cfd6ccd6",
        "litcurate_aa7e18d94a2e00f1",
        "litcurate_2c92df1a000e5fdc",
        "litcurate_4e5c23611ff4a5ae",
    }
    assert all(identifier in text for identifier in identifiers)
