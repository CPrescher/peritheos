"""Frozen release sources and explicitly reviewed integration fingerprints."""

import copy
import hashlib
import json
from pathlib import Path

from peritheos import load_eosmat
from scripts.merge_source_bindings import accepts_merge_binding

ROOT = Path(__file__).resolve().parents[1]


def test_sources_match_frozen_release_or_reviewed_merge_fingerprints():
    proof = json.loads(
        (ROOT / "docs/data/assessment-release-source-bindings.json").read_text()
    )
    assert proof["base_commit"] == "239b47be50e1d65064b118d5cae76e055a3e7ae6"
    assert len(proof["source_documents"]) == 233
    for relative, expected in proof["source_documents"].items():
        document = copy.deepcopy(load_eosmat(ROOT / relative))
        if document["identifier"] == "cementite":
            # Version 0.17.2 adds this one reviewed source and two datasets.
            # Project only those additions out before checking the frozen
            # 0.17.0 source baseline; reviewed integration differences require
            # explicit original-to-merged fingerprint pairs below.
            added_record = "cementite_prescher_2012_full_range_bm3"
            added_datasets = {
                "cementite_prescher_2012_table_s1",
                "cementite_prescher_2012_figure1_digitized",
            }
            assert (
                sum(r["identifier"] == added_record for r in document["eos_records"])
                == 1
            )
            assert added_datasets <= {d["identifier"] for d in document["datasets"]}
            document["eos_records"] = [
                r for r in document["eos_records"] if r["identifier"] != added_record
            ]
            document["datasets"] = [
                d for d in document["datasets"] if d["identifier"] not in added_datasets
            ]
        for record in document["eos_records"]:
            record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)
        actual = hashlib.sha256(
            json.dumps(document, sort_keys=True, separators=(",", ":")).encode()
        ).hexdigest()
        assert actual == expected or accepts_merge_binding(
            "release_document", relative, expected, actual
        ), relative
    for relative, expected in proof["source_resources"].items():
        actual = hashlib.sha256((ROOT / relative).read_bytes()).hexdigest()
        assert actual == expected or accepts_merge_binding(
            "release_resource", relative, expected, actual
        ), relative
