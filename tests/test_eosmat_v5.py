"""V5 parameter-unit uncertainties, physical preservation and old-file readers."""

import copy
import json
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    convert_eosmat,
    create_eosmat,
    eosmat_schema,
    load_eosmat,
    save_eosmat,
    validate_eosmat_document,
)
from peritheos.errors import EosmatError
from scripts.migrate_eosmat_v5 import prior_material_bytes, restore_document

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = ROOT / "crates/peritheos/tests/data/eosmat-v5-uncertainty.json"
SCALE = 6.02214076e23 * 1e-25 / 4


def fixture():
    return json.loads(FIXTURE.read_bytes())


def load(document):
    return Material.from_eosmat(document, require_primary_validation=False)


def test_native_rule_is_simple_and_matches_an_independent_pressure_uncertainty():
    document = fixture()
    Draft202012Validator(eosmat_schema(5)).validate(document)
    record = load(document).eos_records[0]
    assert record.parameter_errors["rt_eos.V0"] == pytest.approx(0.2 * SCALE)
    # At V0,Tref: dP/dV0 = K0/V0 = 1 GPa/Å³; all other derivatives vanish.
    assert record.pressure_with_uncertainty(100, 300).standard_error == pytest.approx(
        0.2, rel=1e-7
    )
    expected = np.array(
        [
            [0.04 * SCALE**2, 0.1 * SCALE, 0.002 * SCALE],
            [0.1 * SCALE, 4, 0.01],
            [0.002 * SCALE, 0.01, 0.01],
        ]
    )
    np.testing.assert_allclose(record.parameter_covariance, expected, rtol=1e-15)
    exported = load(document).to_eosmat()
    raw = exported["eos_records"][0]
    assert raw["parameters"]["V0"]["error"] == 0.2
    assert (
        raw["parameter_covariance"]
        == document["eos_records"][0]["parameter_covariance"]
    )
    assert "parameter_error_volume_basis" not in raw
    assert "volume_basis" not in raw["parameter_covariance"]


@pytest.mark.parametrize("version", [3, 4])
def test_legacy_export_is_old_schema_compatible_and_preserves_predictions(version):
    document = fixture()
    legacy = convert_eosmat(document, version=version)
    raw = convert_eosmat(legacy, version=3, legacy_uncertainty_units=False)[
        "eos_records"
    ][0]
    assert raw["parameter_errors"]["V0"] == pytest.approx(0.2 * SCALE)
    assert raw["parameter_covariance"]["matrix"][0][0] == pytest.approx(0.04 * SCALE**2)
    assert "parameter_error_volume_basis" not in raw
    assert "volume_basis" not in raw["parameter_covariance"]
    schema = eosmat_schema(version)
    schema["$defs"]["eos_record"]["properties"].pop("parameter_error_volume_basis")
    schema["$defs"]["parameter_covariance"]["oneOf"][1]["properties"].pop(
        "volume_basis"
    )
    Draft202012Validator(schema).validate(legacy)
    # Re-saving an unmarked old file must not add unit fields that its original
    # reader's closed schema does not understand.
    Draft202012Validator(schema).validate(load(legacy).to_eosmat())
    first, second = load(document).eos_records[0], load(legacy).eos_records[0]
    for volume, temperature in ((100, 300), (90, 1200)):
        a = first.pressure_with_uncertainty(volume, temperature)
        b = second.pressure_with_uncertainty(volume, temperature)
        assert b.value == a.value
        assert b.standard_error == pytest.approx(a.standard_error, rel=1e-12)


@pytest.mark.parametrize("version", [3, 4])
def test_unmarked_old_files_keep_model_units_and_migrate_explicitly(version, tmp_path):
    legacy = convert_eosmat(fixture(), version=version)
    original = copy.deepcopy(legacy)
    path = tmp_path / "old.eosmat"
    save_eosmat(path, legacy)
    assert load_eosmat(path) == original
    native = convert_eosmat(legacy, version=5)
    assert legacy == original
    assert native["format_version"] == 5
    raw = native["eos_records"][0]
    assert raw["parameters"]["V0"]["error"] == pytest.approx(0.2)
    np.testing.assert_allclose(
        raw["parameter_covariance"]["matrix"],
        fixture()["eos_records"][0]["parameter_covariance"]["matrix"],
        rtol=1e-15,
    )
    np.testing.assert_allclose(
        load(native).eos_records[0].parameter_covariance,
        load(legacy).eos_records[0].parameter_covariance,
        rtol=1e-15,
    )


def test_mixed_v4_annotations_remain_readable_and_normalize_without_guessing():
    old = convert_eosmat(fixture(), version=4)
    old["parameter_error_volume_basis"] = "ignored extension at material scope"
    raw = old["eos_records"][0]
    raw["parameters"]["V0"]["error"] = 0.2
    raw["parameter_error_volume_basis"] = "public"
    raw["parameter_covariance"]["volume_basis"] = "model"
    migrated = convert_eosmat(old)
    assert migrated["eos_records"][0]["parameters"]["V0"]["error"] == 0.2
    assert (
        migrated["parameter_error_volume_basis"] == old["parameter_error_volume_basis"]
    )
    np.testing.assert_allclose(
        load(old).eos_records[0].parameter_covariance,
        load(migrated).eos_records[0].parameter_covariance,
        rtol=1e-15,
    )


@pytest.mark.parametrize("covariance", [False, True])
def test_v5_rejects_legacy_unit_switches(covariance):
    document = fixture()
    raw = document["eos_records"][0]
    if covariance:
        raw["parameter_covariance"]["volume_basis"] = "public"
    else:
        raw["parameter_error_volume_basis"] = "public"
    with pytest.raises(EosmatError, match="always use parameter units"):
        validate_eosmat_document(document)
    assert not Draft202012Validator(eosmat_schema(5)).is_valid(document)


def test_zero_null_and_fixed_annotations_and_extensions_survive_migration():
    document = fixture()
    record = document["eos_records"][0]
    record["parameters"]["K0"]["error"] = 0
    record["extension"] = {"tokens": [None, -0.0, "0.200"]}
    native = convert_eosmat(document)
    assert native == document
    old = convert_eosmat(native, version=4)
    restored = convert_eosmat(old)
    assert restored["eos_records"][0]["parameters"]["K0"] == record["parameters"]["K0"]
    assert (
        restored["eos_records"][0]["parameters"]["K0_prime"]
        == record["parameters"]["K0_prime"]
    )
    assert restored["eos_records"][0]["extension"] == record["extension"]


def test_catalog_migration_preserves_every_predecessor_and_runtime_uncertainty():
    receipt = json.loads((ROOT / "docs/data/eosmat-v5-migration.json").read_bytes())
    for relative, binding in receipt["materials"].items():
        current = load_eosmat(ROOT / relative)
        assert current["format_version"] == 5
        previous = json.loads(prior_material_bytes(Path(relative).stem))
        assert restore_document(current, binding) == previous
        assert current.get("datasets") == previous.get("datasets")
        if not current["eos_records"]:
            continue
        old, new = load(previous), load(current)
        for a, b in zip(old.eos_records, new.eos_records):
            assert a.eos.parameter_values() == b.eos.parameter_values()
            assert (
                a.parameter_errors.keys() == b.parameter_errors.keys()
                if a.parameter_errors is not None
                else b.parameter_errors is None
            )
            if a.parameter_errors:
                np.testing.assert_allclose(
                    list(a.parameter_errors.values()),
                    list(b.parameter_errors.values()),
                    rtol=3e-15,
                    atol=0,
                )
            if a.parameter_covariance is not None:
                np.testing.assert_allclose(
                    a.parameter_covariance, b.parameter_covariance, rtol=3e-15, atol=0
                )
            assert a.scientific_validation_status == b.scientific_validation_status


def test_creation_defaults_to_v5_and_legacy_outputs_keep_working_units():
    for version in (3, 4, 5):
        document = create_eosmat(
            "example",
            "MgO",
            "BM3",
            {"V0": 100, "K0": 100, "K0_prime": 4},
            temperature_ref=300,
            parameter_errors={"V0": 0.2},
            formula_units_per_cell=4,
            thermal={
                "model": "MieGruneisenDebye",
                "parameters": {"Tr": 300, "theta0": 760, "gamma0": 1.5, "q": 1, "n": 2},
            },
            format_version=version,
        )
        assert document["format_version"] == version
        assert load(document).eos_records[0].parameter_errors[
            "rt_eos.V0"
        ] == pytest.approx(0.2 * SCALE)


def test_legacy_v2_still_loads_and_can_be_exported_through_material(tmp_path):
    legacy = convert_eosmat(fixture(), version=3)
    legacy.pop("format")
    legacy["format_version"] = 2
    path = tmp_path / "legacy-v2.eosmat"
    save_eosmat(path, legacy)
    assert load_eosmat(path) == legacy
    material = load(legacy)
    migrated = material.to_eosmat(format_version=5)
    assert load(migrated).eos_records[0].pressure_with_uncertainty(
        100, 300
    ).standard_error == pytest.approx(0.2, rel=1e-7)


@pytest.mark.parametrize("z", [0, -1, None, "four", True])
def test_v5_rejects_unusable_volume_normalization_with_typed_errors(z):
    document = fixture()
    document["formula_units_per_cell"] = z
    with pytest.raises(EosmatError):
        validate_eosmat_document(document)
