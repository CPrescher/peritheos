import importlib.util
import json
from pathlib import Path

from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def test_fu_2023_reproduction_and_schemas():
    path = ROOT / "scripts/reproduce_fu_2023_eos.py"
    spec = importlib.util.spec_from_file_location("reproduce_fu_2023", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    assert len(module.reproduce()["records"]) == 3
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    validator = Draft202012Validator(schema)
    for material in ("mg088fe010al014si090o3_bridgmanite", "ca_perovskite"):
        assert list(validator.iter_errors(get_material_document(material))) == []


def test_fu_2023_audit_enumerates_all_candidates():
    text = (ROOT / "docs/literature-reproductions/fu-2023-eos.md").read_text(
        encoding="utf-8"
    )
    identifiers = {
        "litcurate_5661c5f835e32b74",
        "litcurate_efb817fde2eec8ae",
        "litcurate_1347395c4b56d29a",
        "litcurate_a70856250f14a2ee",
        "litcurate_7a71575938f3f9fb",
        "litcurate_a2591f9f56f627a0",
        "litcurate_438c7d5970e1cd0f",
        "litcurate_a6ab9ce1494d1f92",
        "litcurate_b601d5b21bd6f007",
        "litcurate_a3d1728dde300a01",
        "litcurate_6f0c0a415873bd49",
        "litcurate_dc5880a87d2a1522",
        "litcurate_f150ba8e2eb3c5db",
        "litcurate_7e1c43165d898ecb",
        "litcurate_57faf342bab1540f",
        "litcurate_c570565a82d28af4",
        "litcurate_62c01575e6a31000",
        "litcurate_8f536afdfadf777c",
        "litcurate_b46a77f6e6173d23",
    }
    assert len(identifiers) == 19
    assert all(identifier in text for identifier in identifiers)
