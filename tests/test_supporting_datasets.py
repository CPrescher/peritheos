"""Observation-only studies must survive without invented executable records."""

import copy

import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    eosmat_schema,
    get_material_document,
    validate_eosmat_document,
)


def test_standalone_dataset_survives_material_and_record_selection():
    document = get_material_document("b4c")
    supporting = copy.deepcopy(document["datasets"][0])
    supporting["identifier"] = "supporting_observations"
    supporting["used_by_eos_records"] = []
    document["datasets"].append(supporting)
    assert not list(Draft202012Validator(eosmat_schema()).iter_errors(document))
    validate_eosmat_document(document)
    selected = document["eos_records"][0]["identifier"]
    material = Material.from_eosmat(document, record_identifiers=[selected])
    dataset = material.get_dataset("supporting_observations")
    assert dataset.used_by_eos_records == ()
    assert len(dataset) == len(supporting["rows"])
    restored = Material.from_eosmat(material.to_eosmat())
    assert restored.get_dataset("supporting_observations").used_by_eos_records == ()
    observation_only = copy.deepcopy(document)
    observation_only["eos_records"] = []
    observation_only["datasets"] = [supporting]
    validate_eosmat_document(observation_only)
    assert not list(Draft202012Validator(eosmat_schema()).iter_errors(observation_only))
    # Observation-only documents remain non-executable by design.
    with pytest.raises(ValueError, match="at least one EOS record"):
        Material.from_eosmat(observation_only)
    # Allowing no links must not allow malformed or dangling links.
    for invalid in (None, "", ["missing_record"]):
        document["datasets"][-1]["used_by_eos_records"] = invalid
        with pytest.raises(ValueError):
            validate_eosmat_document(document)
