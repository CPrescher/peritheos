"""Later calculated declarations must not invalidate or replace frozen sources."""

import copy
import hashlib
import json

from peritheos import load_eosmat
from peritheos_build import ROOT
from scripts.merge_source_bindings import accepts_merge_binding


def fingerprint(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def test_removing_only_reviewed_additions_recovers_accepted_predecessors():
    proof = json.loads(
        (ROOT / "docs/data/dataset-addition-source-bindings.json").read_bytes()
    )
    assert (
        sum(len(e["appended_declarations"]) for e in proof["materials"].values()) == 7
    )
    for entry in proof["materials"].values():
        path = entry["path"]
        document = copy.deepcopy(load_eosmat(ROOT / path))
        for record in document["eos_records"]:
            record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)
        binding = proof["bindings"]["release_document"][path]
        actual = fingerprint(document)
        assert actual == binding["current_sha256"]
        assert accepts_merge_binding(
            "release_document", path, binding["historical_sha256"], actual
        )
        assert not accepts_merge_binding("release_document", path, "0" * 64, actual)
        assert not accepts_merge_binding(
            "unreviewed", path, binding["historical_sha256"], actual
        )
        for dataset in document["datasets"]:
            if dataset["identifier"] in entry["appended_declarations"]:
                assert (
                    fingerprint(dataset)
                    == entry["appended_declarations"][dataset["identifier"]]
                )
        document["datasets"] = [
            d
            for d in document["datasets"]
            if d["identifier"] not in entry["appended_declarations"]
        ]
        previous = fingerprint(document)
        assert previous == entry["previous_reviewed_sha256"]
        assert previous == binding["historical_sha256"] or accepts_merge_binding(
            "release_document", path, binding["historical_sha256"], previous
        )
        document["eos_records"][0]["parameters"]["K0"]["value"] += 1
        assert not accepts_merge_binding(
            "release_document",
            path,
            binding["historical_sha256"],
            fingerprint(document),
        )
