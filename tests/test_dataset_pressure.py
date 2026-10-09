"""Dataset pressure contracts, independent numerical checks, and interchange."""

import copy
import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from peritheos import (
    DatasetError,
    EosmatError,
    Material,
    eosmat_schema,
    get_material_document,
    load_eosmat,
    resolve_dataset_pressure,
    save_eosmat,
    validate_eosmat_document,
)
from scripts.audit_dewaele_2019_static_dac import SOURCES, _fit_pressure, _published

ROOT = Path(__file__).resolve().parents[1]
FIXTURE = json.loads(
    (ROOT / "crates/peritheos/tests/data/dataset-pressure.json").read_text()
)
METALS = {
    "aluminum": 40,
    "copper": 42,
    "gold": 37,
    "platinum": 36,
    "tantalum": 36,
    "tungsten": 42,
}


def table(document):
    return next(
        d for d in document["datasets"] if "dewaele_2004_table1" in d["identifier"]
    )


@pytest.mark.parametrize("metal,count", METALS.items())
def test_migrated_reductions_preserve_rows_scales_and_exports(metal, count, tmp_path):
    document = get_material_document(metal)
    dataset = table(document)
    assert sum(r["status"] != "unresolved" for r in dataset["pressure_reductions"]) == (
        2 if metal == "gold" else 4
    )
    assert [
        c["pressure_scale"]["calibration_record"] for c in dataset["columns"][:2]
    ] == ["ruby_mao_1986", "ruby_dewaele_2004"]
    jsonschema.validate(document, eosmat_schema())
    assert (
        Material.from_eosmat(document, require_primary_validation=False).to_eosmat()[
            "datasets"
        ]
        == document["datasets"]
    )
    path = tmp_path / "material.eosmat"
    save_eosmat(path, document)
    assert load_eosmat(path) == document
    for reduction in dataset["pressure_reductions"]:
        selected = Material.from_eosmat(
            document, record_identifiers=(reduction["eos_record"],)
        ).to_eosmat()
        assert table(selected)["pressure_reductions"] == [reduction]
        jsonschema.validate(selected, eosmat_schema())
        if reduction["status"] == "unresolved":
            continue
        result = resolve_dataset_pressure(
            document, reduction["eos_record"], dataset["identifier"]
        )
        selected_pressure = resolve_dataset_pressure(
            selected, reduction["eos_record"], dataset["identifier"]
        )
        np.testing.assert_array_equal(
            selected_pressure.pressure_gpa, result.pressure_gpa
        )
        assert len(result.pressure_gpa) == count
        np.testing.assert_array_equal(result.row_indices, np.arange(count))
        assert result.dataset.checksum_verified
        np.testing.assert_array_equal(
            result.reported_pressure, result.dataset[reduction["pressure_column"]]
        )
        assert not result.pressure_gpa.flags.writeable
        assert not result.row_indices.flags.writeable
        if reduction["status"] == "as_reported":
            np.testing.assert_array_equal(result.pressure_gpa, result.reported_pressure)
        else:
            # Independent equation from the primary-source audit, all rows.
            ratio = (1 + 7.665 * result.reported_pressure / 1904) ** (1 / 7.665)
            expected = 1884 * (ratio - 1) * (1 + 5.5 * (ratio - 1))
            np.testing.assert_allclose(result.pressure_gpa, expected, rtol=2e-14)
            assert not np.allclose(
                result.pressure_gpa, result.dataset["ruby_pressure_revised_gpa"]
            )


@pytest.mark.parametrize("metal", set(METALS) - {"gold"})
@pytest.mark.parametrize("scale", ["mao", "dor"])
def test_complete_2019_reductions_reproduce_existing_audit(metal, scale):
    document = get_material_document(metal)
    spec = next(s for s in SOURCES if s.material_file == f"{metal}.eosmat")
    identifier, published, errors = _published(spec, scale)
    result = resolve_dataset_pressure(
        document, identifier, table(document)["identifier"]
    )
    fit, _, _ = _fit_pressure(
        result.dataset["atomic_volume_a3"], result.pressure_gpa, published
    )
    audit = json.loads(
        (ROOT / "docs/data/dewaele-2019-static-dac-refit.json").read_text()
    )
    expected = audit["row_level_refits"][identifier][
        "unweighted_pressure_residual_fit"
    ]["parameters"]
    np.testing.assert_allclose(fit, expected, atol=3e-5, rtol=1e-6)
    assert np.all(abs(fit - published) <= errors)


def test_aluminum_validity_temperatures_and_uncertainties_remain_distinct():
    document = get_material_document("aluminum")
    result = resolve_dataset_pressure(
        document, "aluminum_dewaele_2019_dor_vinet", table(document)["identifier"]
    )
    assert len(result.pressure_gpa) == 40
    assert result.source_calibration["identifier"] == "ruby_mao_1986"
    assert result.target_calibration["identifier"] == "ruby_dorogokupets_oganov_2007"
    assert result.source_calibration["validity"]["temperature_k"] == [298, 298]
    assert result.target_calibration["validity"]["temperature_k"] == [298, 298]
    assert result.eos_reference_temperature_k == 300
    assert not result.dataset.find_columns(quantity="temperature")
    assert result.uncertainty_treatment == "not_propagated"
    assert result.source_pressure_extrapolated.sum() == 12
    assert not result.target_pressure_extrapolated.any()
    assert result.pressure_gpa[-1] == pytest.approx(155.29074838, abs=1e-8)
    assert result.dataset["atomic_volume_uncertainty_a3"][0] == 0.01
    assert result.dataset.get_column("atomic_volume_uncertainty_a3").role == "bound"
    with pytest.raises(DatasetError, match="validity"):
        resolve_dataset_pressure(
            document, result.eos_record, result.dataset.identifier, check_validity=True
        )


def test_printed_copper_anomaly_is_not_recomputed():
    document = get_material_document("copper")
    result = resolve_dataset_pressure(
        document, "copper_dewaele_2004_vinet_1", table(document)["identifier"]
    )
    index = np.flatnonzero(result.dataset["ruby_pressure_classical_gpa"] == 37)[0]
    assert result.pressure_gpa[index] == 37.1


@pytest.mark.parametrize(
    "record",
    ["platinum_dorogokupets_oganov_2007_vinet_4", "platinum_fei_2007_vinet_300k"],
)
def test_reference_model_links_do_not_authorize_reduction(record):
    document = get_material_document("platinum")
    with pytest.raises(DatasetError, match="Unresolved"):
        resolve_dataset_pressure(document, record, table(document)["identifier"])


@pytest.mark.parametrize("case", FIXTURE["invalid_cases"], ids=lambda c: c["name"])
def test_invalid_metadata_shared_with_rust(case):
    document = copy.deepcopy(FIXTURE["document"])
    parts = case["path"].strip("/").split("/")
    target = document
    for part in parts[:-1]:
        target = target[int(part)] if isinstance(target, list) else target[part]
    target[int(parts[-1]) if isinstance(target, list) else parts[-1]] = case["value"]
    with pytest.raises(EosmatError):
        validate_eosmat_document(document)


def test_generic_names_units_and_backward_compatibility():
    document = copy.deepcopy(FIXTURE["document"])
    result = resolve_dataset_pressure(document, "test_eos", "observations")
    assert len(result.pressure_gpa) == 2
    document["datasets"][0]["columns"][0]["unit"] = "MPa"
    document["datasets"][0]["rows"] = [[1000, 9], [100000, 7]]
    converted = resolve_dataset_pressure(document, "test_eos", "observations")
    np.testing.assert_array_equal(converted.pressure_gpa, result.pressure_gpa)
    assert converted.reported_pressure_unit == "MPa"
    np.testing.assert_array_equal(converted.reported_pressure, [1000, 100000])
    del document["datasets"][0]["pressure_reductions"]
    del document["datasets"][0]["columns"][0]["pressure_scale"]
    validate_eosmat_document(document)
    jsonschema.validate(document, eosmat_schema())
    with pytest.raises(DatasetError, match="No explicit reduction"):
        resolve_dataset_pressure(document, "test_eos", "observations")
    with pytest.raises(DatasetError, match="Unknown dataset"):
        resolve_dataset_pressure(document, "test_eos", "missing")


@pytest.mark.parametrize("value", [None, -1])
def test_unusable_source_pressures_are_not_dropped(value):
    document = copy.deepcopy(FIXTURE["document"])
    document["datasets"][0]["rows"][0][0] = value
    with pytest.raises(DatasetError, match="finite non-negative"):
        resolve_dataset_pressure(document, "test_eos", "observations")


@pytest.mark.parametrize("status", ["as_reported", "transformed", "unresolved"])
def test_schema_variants(status):
    document = copy.deepcopy(FIXTURE["document"])
    reduction = document["datasets"][0]["pressure_reductions"][0]
    if status == "as_reported":
        document["datasets"][0]["columns"][0]["pressure_scale"][
            "calibration_record"
        ] = reduction["target_calibration_record"]
        del reduction["convention"]
    if status == "unresolved":
        for key in [
            "convention",
            "pressure_column",
            "target_calibration_record",
            "row_scope",
        ]:
            del reduction[key]
        reduction["notes"] = "Source does not establish the coordinate."
    reduction["status"] = status
    jsonschema.validate(document, eosmat_schema())
    validate_eosmat_document(document)
    reduction["convention"] = "guess"
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(document, eosmat_schema())


def test_calibration_kind_must_match_eos_method():
    from peritheos import list_diamond_raman_calibrations

    document = copy.deepcopy(FIXTURE["document"])
    diamond = list_diamond_raman_calibrations()[0]
    reduction = document["datasets"][0]["pressure_reductions"][0]
    document["datasets"][0]["columns"][0]["pressure_scale"]["calibration_record"] = (
        diamond
    )
    reduction.update(status="as_reported", target_calibration_record=diamond)
    del reduction["convention"]
    document["eos_records"][0]["pressure_calibration"]["methods"][0][
        "reference_calibration_record"
    ] = diamond
    with pytest.raises(EosmatError, match="conflicts"):
        validate_eosmat_document(document)
