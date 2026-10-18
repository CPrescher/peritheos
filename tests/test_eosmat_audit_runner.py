"""Historical scientific scripts keep their fingerprints and can edit v4 safely."""

import json
import sys
from types import ModuleType

import pytest

from peritheos import get_material_document, load_eosmat
from scripts.eosmat_io import write_material
from scripts.run_eosmat_audit import run_audit


@pytest.mark.parametrize("outcome", ["check", "update", "failure", "concurrent"])
def test_audit_uses_temporary_legacy_file_and_restores_native_state(tmp_path, outcome):
    path = tmp_path / "gold.eosmat"
    write_material(path, get_material_document("gold"))
    original = path.read_bytes()
    module = ModuleType("historical_audit")
    module.MATERIAL = path
    argv = sys.argv

    def main():
        assert module.MATERIAL != path
        document = json.loads(module.MATERIAL.read_text())
        assert document["format_version"] == 3
        assert sys.argv == ["historical_audit", "--check"]
        if outcome == "check":
            return
        document["audit_extension"] = {"ordered": [3, 1, 2]}
        module.MATERIAL.write_text(json.dumps(document))
        if outcome == "failure":
            raise RuntimeError("Audit failed")
        if outcome == "concurrent":
            path.write_bytes(original + b"\n")

    module.main = main
    if outcome in {"failure", "concurrent"}:
        with pytest.raises(RuntimeError):
            run_audit(module, "MATERIAL", ["--check"])
    else:
        run_audit(module, "MATERIAL", ["--check"])
    assert module.MATERIAL == path
    assert sys.argv is argv
    if outcome == "update":
        updated = load_eosmat(path)
        assert updated["format_version"] == 4
        assert updated.pop("audit_extension") == {"ordered": [3, 1, 2]}
        assert updated == json.loads(original)
    else:
        assert path.read_bytes() == original + (
            b"\n" if outcome == "concurrent" else b""
        )
