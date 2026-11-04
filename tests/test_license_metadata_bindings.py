"""License-only pins must recover the previous complete scientific declaration."""

import copy
import hashlib
import json
from pathlib import Path

import pytest

from peritheos import load_eosmat
from scripts.merge_source_bindings import accepts_merge_binding

ROOT = Path(__file__).resolve().parents[1]
PROOF = json.loads(
    (ROOT / "docs/data/license-metadata-source-bindings.json").read_text(
        encoding="utf-8"
    )
)


def semantic(value):
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    ).hexdigest()


def verify_license_only_change(document, review):
    restored = copy.deepcopy(document)
    for record in restored["eos_records"]:
        record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)
    assert semantic(restored) == review["licensed_document_sha256"]
    datasets = {d["identifier"]: d for d in restored["datasets"]}
    for edit in review["dataset_edits"]:
        dataset = datasets[edit["dataset_id"]]
        assert {k: dataset[k] for k in PROOF["license_fields"] if k in dataset} == edit[
            "after"
        ]
        for key in PROOF["license_fields"]:
            dataset.pop(key, None)
        dataset.update(edit["before"])
    # Every coefficient, observation, source hash and other field is still bound
    # by the previous whole-document digest, not just the new licensed digest.
    assert semantic(restored) == review["previous_document_sha256"]


@pytest.mark.parametrize("material", sorted(PROOF["materials"]))
def test_only_the_exact_reviewed_license_fields_changed(material):
    review = PROOF["materials"][material]
    binding = PROOF["bindings"]["release_document"][review["path"]]
    assert review["previous_document_sha256"] == binding["previous_reviewed_sha256"]
    assert review["licensed_document_sha256"] == binding["licensed_sha256"]
    verify_license_only_change(load_eosmat(ROOT / review["path"]), review)


@pytest.mark.parametrize("material", sorted(PROOF["materials"]))
def test_scientific_change_cannot_be_hidden_by_repinning_the_new_digest(material):
    review = copy.deepcopy(PROOF["materials"][material])
    document = load_eosmat(ROOT / review["path"])
    document["eos_records"][0]["parameters"]["V0"]["value"] += 1e-12
    stripped = copy.deepcopy(document)
    for record in stripped["eos_records"]:
        record.get("evidence", {}).get("reproduction", {}).pop("assessment", None)
    review["licensed_document_sha256"] = semantic(stripped)
    with pytest.raises(AssertionError):
        verify_license_only_change(document, review)


def test_compatibility_pins_require_the_exact_historical_digest_and_context():
    for context, items in PROOF["bindings"].items():
        for identifier, entry in items.items():
            old, new = entry["historical_sha256"], entry["licensed_sha256"]
            assert accepts_merge_binding(context, identifier, old, new)
            assert not accepts_merge_binding(context, identifier, "0" * 64, new)
            assert not accepts_merge_binding(context, identifier, old, "0" * 64)
            assert not accepts_merge_binding("unreviewed", identifier, old, new)
