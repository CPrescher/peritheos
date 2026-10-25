"""The assessment release preserves the complete published scientific baseline."""

import copy
import hashlib
import json
from pathlib import Path

from peritheos import load_eosmat

ROOT = Path(__file__).resolve().parents[1]


def test_all_published_material_sources_and_resources_are_unchanged():
    proof = json.loads(
        (ROOT / "docs/data/assessment-release-source-bindings.json").read_text()
    )
    assert proof["base_commit"] == "239b47be50e1d65064b118d5cae76e055a3e7ae6"
    assert len(proof["source_documents"]) == 233
    for relative, expected in proof["source_documents"].items():
        document = copy.deepcopy(load_eosmat(ROOT / relative))
        for record in document["eos_records"]:
            record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)
        actual = hashlib.sha256(
            json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        assert actual == expected, relative
    for relative, expected in proof["source_resources"].items():
        assert hashlib.sha256((ROOT / relative).read_bytes()).hexdigest() == expected, (
            relative
        )
