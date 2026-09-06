import importlib.util
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_karki_reproduction_and_schema():
    path = ROOT / "scripts/reproduce_karki_1997_mgo.py"
    spec = importlib.util.spec_from_file_location("reproduce_karki_1997", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    result = module.reproduce()
    assert len(result["records"]) == 2
    assert all(
        row["max_abs_equation_residual_gpa"] < 2e-11
        for row in result["records"].values()
    )
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    assert (
        list(Draft202012Validator(schema).iter_errors(get_material_document("mgo")))
        == []
    )


def test_karki_audit_enumerates_all_candidates():
    text = (ROOT / "docs/literature-reproductions/karki-1997-mgo.md").read_text(
        encoding="utf-8"
    )
    identifiers = {
        "litcurate_c2d9fd33177eb390",
        "litcurate_df447e7b3becd6df",
        "litcurate_e6b89c003860499c",
        "litcurate_63c5aca81ea6f807",
        "litcurate_25c45c2113a30788",
        "litcurate_4e7138b29be1cf65",
        "litcurate_02a07b79e480ba9e",
        "litcurate_25ffc08073a3042c",
        "litcurate_64446b02fc337884",
        "litcurate_6a40e1d5485ae1cb",
        "litcurate_399535c9578d02c7",
    }
    assert all(identifier in text for identifier in identifiers)
