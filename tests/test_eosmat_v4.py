"""Native v4 files are the scientific record, with complete legacy migration."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    EosmatError,
    Material,
    convert_eosmat,
    create_eosmat,
    eosmat_schema,
    get_material_document,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
)
from peritheos.__main__ import main
from peritheos.observations import observation_material_fingerprint

ROOT = Path(__file__).resolve().parents[1]


def semantic(document):
    """Only optional empty metadata maps/lists have equivalent absent spellings."""
    if isinstance(document, dict):
        return {
            k: semantic(v)
            for k, v in document.items()
            if not (k in {"parameter_errors", "fixed_parameters"} and not v)
        }
    if isinstance(document, list):
        return [semantic(v) for v in document]
    return document


def test_entire_catalog_migrates_with_every_scientific_field_and_source_array_intact():
    validator = Draft202012Validator(eosmat_schema(4))
    paths = sorted((ROOT / "peritheos/data/materials").glob("*.eosmat"))
    assert len(paths) == 232
    for path in paths:
        original = json.loads(path.read_text())
        before = copy.deepcopy(original)
        native = convert_eosmat(original)
        validator.validate(native)
        restored = convert_eosmat(native, version=3)
        assert semantic(restored) == semantic(original), path.name
        assert observation_material_fingerprint(
            original
        ) == observation_material_fingerprint(native)
        assert native.get("datasets") == original.get("datasets"), path.name
        assert original == before


@pytest.mark.parametrize("identifier", ["gold", "mgo", "iron", "diamond", "nacl_b2"])
def test_native_loading_keeps_numerical_behavior_uncertainty_defaults_and_restrictions(
    identifier,
):
    original = get_material_document(identifier)
    native = convert_eosmat(original)
    old = Material.from_eosmat(original, require_primary_validation=False)
    new = Material.from_eosmat(native, require_primary_validation=False)
    assert [r.identifier for r in old.eos_records] == [
        r.identifier for r in new.eos_records
    ]
    assert new.to_eosmat()["format_version"] == 4
    assert new.to_eosmat(format_version=3)["format_version"] == 3
    assert (
        old.default_equilibrium_record().identifier
        == new.default_equilibrium_record().identifier
    )
    for before, after in zip(old.eos_records, new.eos_records):
        assert before.parameter_errors == after.parameter_errors
        assert before.eosmat_metadata.get(
            "fixed_parameters"
        ) == after.eosmat_metadata.get("fixed_parameters")
        assert before.scientific_validation_status == after.scientific_validation_status
        assert before.reference_volume == after.reference_volume
        assert before.reference_temperature == after.reference_temperature
        # The same evaluation and domain errors must occur for qualified models.
        volumes = np.asarray([0.99, 1.0]) * before.reference_volume
        try:
            expected = before.pressure(volumes, before.reference_temperature)
        except (ValueError, ArithmeticError, NotImplementedError) as error:
            with pytest.raises(type(error)):
                after.pressure(volumes, after.reference_temperature)
        else:
            np.testing.assert_array_equal(
                after.pressure(volumes, after.reference_temperature), expected
            )


def test_native_round_trip_has_no_hidden_generated_file_and_preserves_extensions(
    tmp_path,
):
    document = create_eosmat(
        "test",
        "X",
        "BM3",
        {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0},
        temperature_ref=300,
    )
    document["extension"] = {"tokens": [None, -0.0, 1.2345678901234567]}
    record = document["eos_records"][0]
    record["extension"] = {"original_order": [3, 1, 2]}
    record["parameters"]["K0"] = {"value": 160.0, "error": 0.0, "fixed": True}
    path = tmp_path / "native.eosmat"
    save_eosmat(path, document)
    assert load_eosmat(path) == document
    assert list(tmp_path.iterdir()) == [path]
    assert "format" not in document and "units" not in document
    assert "eos" not in record and "parameter_errors" not in record
    with pytest.raises(ValueError, match="pending_primary_source_check"):
        Material.from_eosmat(document)
    assert (
        convert_eosmat(document, version=3)["eos_records"][0]["parameter_errors"]["K0"]
        == 0
    )


@pytest.mark.parametrize(
    "change",
    [
        lambda d: d.update(format="peritheos.material"),
        lambda d: d["eos_records"][0].update(eos={}),
        lambda d: d["eos_records"][0].update(model="Unknown"),
        lambda d: d["eos_records"][0]["parameters"].update(
            V0={"value": 10.0, "fixed": 1}
        ),
        lambda d: d["eos_records"][0]["parameters"].update(K0={"fixed": False}),
        lambda d: d["eos_records"][0].update(equation_options={"model": "BM2"}),
        lambda d: d["eos_records"][0].update(evidence={"identifier": "other"}),
        lambda d: d["eos_records"][0].update(evidence={"model": "BM2"}),
        lambda d: d["eos_records"][0].update(evidence={"default_for": "equilibrium"}),
        lambda d: d["eos_records"][0].update(evidence={"thermal": {}}),
        lambda d: d.update(units={"volume": "cm^3/mol"}),
    ],
)
def test_ambiguous_or_conflicting_native_fields_fail_closed(change):
    document = create_eosmat(
        "test",
        "X",
        "BM3",
        {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0},
        temperature_ref=300,
    )
    change(document)
    assert not Draft202012Validator(eosmat_schema(4)).is_valid(document)
    with pytest.raises(EosmatError):
        validate_eosmat_document(document)


def test_migration_command_writes_one_native_file_and_can_export_legacy(tmp_path):
    source = tmp_path / "source.eosmat"
    native = tmp_path / "native.eosmat"
    legacy = tmp_path / "legacy.eosmat"
    original = get_material_document("gold")
    save_eosmat(source, original)
    assert main(["convert", str(source), str(native)]) == 0
    assert load_eosmat(native)["format_version"] == 4
    assert main(["convert", str(native), str(legacy), "--format-version", "3"]) == 0
    assert semantic(load_eosmat(legacy)) == semantic(original)
    with pytest.raises(SystemExit):
        main(["convert", str(source), str(native)])


def test_existing_v3_extension_collision_is_reported_instead_of_discarded():
    document = get_material_document("gold")
    document["eos_records"][0]["model"] = {"unrecognized_extension": 42}
    with pytest.raises(EosmatError, match="collides"):
        convert_eosmat(document)
