"""Uncertainty units, source preservation and legacy EOSMAT compatibility."""

import hashlib
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
    get_eos_record,
    get_material_document,
    validate_eosmat_document,
)
from peritheos.errors import EosmatError
from peritheos.materials import MGO_TANGE_2009
from peritheos_build import digest

ROOT = Path(__file__).resolve().parents[1]


def thermal_document():
    return create_eosmat(
        "uncertainty_example",
        "MgO",
        "BM3",
        {"V0": 74.698, "K0": 160.63, "K0_prime": 4.367},
        temperature_ref=300,
        parameter_errors={"V0": 0.016},
        formula_units_per_cell=4,
        thermal={
            "model": "MieGruneisenDebye",
            "parameters": {"Tr": 300, "theta0": 761, "gamma0": 1.442, "q": 1, "n": 2},
        },
    )


def test_catalog_mgo_uncertainty_matches_closed_form_and_compatibility_record():
    record = get_eos_record("mgo_b1_tange_2009_vinet")
    result = record.pressure_with_uncertainty(74.698, 300)
    # At the reference state only V0 contributes: dP/dV0 = K0/V0.
    expected = 160.63 * 0.016 / 74.698
    assert result.standard_error == pytest.approx(expected, rel=1e-7)
    assert result.standard_error == pytest.approx(
        MGO_TANGE_2009.pressure_with_uncertainty(74.698, 300).standard_error,
        rel=1e-7,
    )


def test_authored_public_volume_error_and_legacy_export_preserve_relative_error():
    document = thermal_document()
    material = Material.from_eosmat(document, require_primary_validation=False)
    record = material.eos_records[0]
    assert record.parameter_errors["rt_eos.V0"] / record.eos.rt_eos.V0 == pytest.approx(
        0.016 / 74.698
    )
    for version in (3, 4, 5):
        exported = material.to_eosmat(format_version=version)
        row = convert_eosmat(exported, version=3, legacy_uncertainty_units=False)[
            "eos_records"
        ][0]
        if version in (3, 4):
            assert "parameter_error_volume_basis" not in row
            assert row["parameter_errors"]["V0"] == pytest.approx(
                0.016 * record.volume_scale
            )
        else:
            assert row["parameter_error_volume_basis"] == "public"
            assert row["parameter_errors"]["V0"] == pytest.approx(0.016)
        reloaded = Material.from_eosmat(exported, require_primary_validation=False)
        assert reloaded.eos_records[0].parameter_errors == pytest.approx(
            record.parameter_errors
        )
        assert reloaded.eos_records[0].pressure_with_uncertainty(
            70, 1000
        ).standard_error == pytest.approx(
            record.pressure_with_uncertainty(70, 1000).standard_error
        )


@pytest.mark.parametrize("basis", ["public", "model"])
def test_covariance_volume_rows_columns_and_cross_terms_scale_explicitly(basis):
    document = convert_eosmat(
        thermal_document(), version=4, legacy_uncertainty_units=False
    )
    raw = document["eos_records"][0]
    raw["parameter_covariance"] = {
        "parameter_order": ["rt_eos.V0", "rt_eos.K0"],
        "volume_basis": basis,
        "matrix": [[0.04, 0.1], [0.1, 4.0]],
    }
    material = Material.from_eosmat(document, require_primary_validation=False)
    record = material.eos_records[0]
    scale = record.volume_scale if basis == "public" else 1.0
    expected = [[0.04 * scale**2, 0.1 * scale], [0.1 * scale, 4.0]]
    np.testing.assert_allclose(record.parameter_covariance, expected)
    for version in (3, 4, 5):
        output = material.to_eosmat(format_version=version)
        reloaded = Material.from_eosmat(output, require_primary_validation=False)
        np.testing.assert_allclose(
            reloaded.eos_records[0].parameter_covariance, expected
        )


def test_unmarked_historical_volume_uncertainties_keep_model_units():
    document = convert_eosmat(
        thermal_document(), version=4, legacy_uncertainty_units=False
    )
    document["eos_records"][0].pop("parameter_error_volume_basis")
    record = Material.from_eosmat(
        document, require_primary_validation=False
    ).eos_records[0]
    assert record.parameter_errors["rt_eos.V0"] == 0.016


def test_v3_uncertainty_exports_validate_with_the_previous_closed_schema():
    schema = eosmat_schema(3)
    schema["$defs"]["eos_record"]["properties"].pop("parameter_error_volume_basis")
    schema["$defs"]["parameter_covariance"]["oneOf"][1]["properties"].pop(
        "volume_basis"
    )
    document = thermal_document()
    document["eos_records"][0]["parameter_covariance"] = {
        "parameter_order": ["rt_eos.V0"],
        "matrix": [[0.02**2]],
    }
    material = Material.from_eosmat(document, require_primary_validation=False)
    legacy = material.to_eosmat(format_version=3)
    Draft202012Validator(schema).validate(legacy)
    row = legacy["eos_records"][0]
    assert "parameter_error_volume_basis" not in row
    assert "volume_basis" not in row["parameter_covariance"]
    direct = create_eosmat(
        "legacy_thermal",
        "MgO",
        "BM3",
        {"V0": 74.7, "K0": 160, "K0_prime": 4},
        temperature_ref=300,
        parameter_errors={"V0": 0.02},
        formula_units_per_cell=4,
        thermal={
            "model": "MieGruneisenDebye",
            "parameters": {"Tr": 300, "theta0": 760, "gamma0": 1.5, "q": 1, "n": 2},
        },
        format_version=3,
    )
    Draft202012Validator(schema).validate(direct)
    assert direct["eos_records"][0]["parameter_errors"]["V0"] == pytest.approx(
        0.02 * material.eos_records[0].volume_scale
    )


def test_historical_recovery_rejects_changes_to_reviewed_source_numbers(
    tmp_path, monkeypatch
):
    from scripts import uncertainty_volume_basis_review as review

    relative = "peritheos/data/materials/mgo.eosmat"
    document = json.loads((ROOT / relative).read_bytes())
    document["eos_records"][0]["parameters"]["V0"] = 123
    path = tmp_path / relative
    path.parent.mkdir(parents=True)
    path.write_text(json.dumps(document))
    monkeypatch.setattr(review, "ROOT", tmp_path)
    with pytest.raises(ValueError, match="source changed"):
        review.prior_material_bytes("mgo")


def test_catalog_volume_unit_annotations_preserve_every_original_source_value():
    receipt = json.loads(
        (ROOT / "docs/data/uncertainty-volume-basis-review.json").read_text()
    )
    for relative, binding in receipt["materials"].items():
        from scripts.migrate_eosmat_v5 import prior_material_bytes

        data = prior_material_bytes(Path(relative).stem)
        assert hashlib.sha256(data).hexdigest() == binding["after_sha256"]
        document = json.loads(data)
        for row in document["eos_records"]:
            row.pop("parameter_error_volume_basis", None)
            if isinstance(row.get("parameter_covariance"), dict):
                row["parameter_covariance"].pop("volume_basis", None)
        assert digest(document) == binding["preserved_document_sha256"]
        native = get_material_document(Path(relative).stem, format_version=5)
        material = Material.from_eosmat(native, require_primary_validation=False)
        before = {
            r["identifier"]: r for r in convert_eosmat(native, version=3)["eos_records"]
        }
        after = {
            r["identifier"]: r
            for r in convert_eosmat(material.to_eosmat(), version=3)["eos_records"]
        }
        for item in receipt["records"]:
            if item["material"] != material.identifier:
                continue
            identifier = item["identifier"]
            for field in ("parameter_errors", "parameter_covariance"):
                assert before[identifier].get(field) == after[identifier].get(field)
    for item in receipt["records"]:
        document = get_material_document(item["material"])
        material = Material.from_eosmat(
            document,
            record_identifiers=[item["identifier"]],
            require_primary_validation=False,
        )
        record = material.eos_records[0]
        if "original_V0_error" in item:
            scale = (
                record.volume_scale if item["error_volume_basis"] == "public" else 1.0
            )
            assert record.parameter_errors["rt_eos.V0"] == pytest.approx(
                item["original_V0_error"] * scale
            )


@pytest.mark.parametrize("version", [3, 4, 5])
@pytest.mark.parametrize("thermal", [False, True])
def test_negative_parameter_errors_are_rejected_by_reader_and_schema(version, thermal):
    document = thermal_document()
    component = document["eos_records"][0]
    if thermal:
        component = component["thermal"]
        component["parameters"]["gamma0"] = {"value": 1.442, "error": -0.01}
    else:
        component["parameters"]["V0"]["error"] = -0.01
    # Construct v3 without the public converter's validation, to test the reader.
    if version in (3, 4):
        from peritheos._eosmat_v4 import to_v3, to_v4

        document = to_v3(document)
        if version == 4:
            document = to_v4(document)
    with pytest.raises(EosmatError, match="nonnegative"):
        validate_eosmat_document(document)
    assert list(Draft202012Validator(eosmat_schema(version)).iter_errors(document))


@pytest.mark.parametrize("covariance", [False, True])
@pytest.mark.parametrize("bad_basis", ["unknown", None, {}, []])
def test_unknown_volume_bases_are_rejected(covariance, bad_basis):
    document = convert_eosmat(
        thermal_document(), version=4, legacy_uncertainty_units=False
    )
    raw = document["eos_records"][0]
    if covariance:
        raw["parameter_covariance"] = {
            "parameter_order": ["rt_eos.V0"],
            "matrix": [[1.0]],
            "volume_basis": bad_basis,
        }
    else:
        raw["parameter_error_volume_basis"] = bad_basis
    with pytest.raises(EosmatError, match="must be 'public' or 'model'"):
        validate_eosmat_document(document)
