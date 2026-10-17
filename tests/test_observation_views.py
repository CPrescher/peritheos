"""Observation contracts preserve source evidence while making coordinates explicit."""

import copy
import importlib.util
import json
import subprocess
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from peritheos import (
    DatasetError,
    EosmatError,
    get_material_document,
    observation_registry,
    plan_dataset_observations,
    resolve_dataset_observations,
    validate_eosmat_document,
)
from peritheos.observations import (
    CELL_UNIT,
    observation_material_fingerprint,
    observation_source_fingerprint,
    validate_observation_view,
)

ROOT = Path(__file__).resolve().parents[1]
SCHEMA = json.loads(
    (ROOT / "peritheos/data/observation-view-v1.schema.json").read_text()
)
EMPTY_REGISTRY = {"datasets": []}


def test_entire_packaged_registry_is_bound_valid_unambiguous_and_resolves():
    registry = observation_registry()
    pairs = set()
    assert len(registry["datasets"]) >= 300
    for entry in registry["datasets"]:
        document = get_material_document(entry["material_id"])
        dataset = next(
            d for d in document["datasets"] if d["identifier"] == entry["dataset_id"]
        )
        assert observation_source_fingerprint(dataset) == entry["declaration_sha256"]
        assert (
            observation_material_fingerprint(document)
            == entry["material_identity_sha256"]
        )
        for view in entry["views"]:
            jsonschema.validate(view, SCHEMA)
            validate_observation_view(
                view, dataset, {r["identifier"] for r in document["eos_records"]}
            )
            for record in view["eos_records"]:
                pair = (entry["material_id"], entry["dataset_id"], record)
                assert pair not in pairs
                pairs.add(pair)
            # Resolve each distinct coordinate interpretation through real packaged data.
            data = resolve_dataset_observations(
                document, view["eos_records"][0], entry["dataset_id"], registry=registry
            )
            assert len(data.row_indices) > 0, pair
            assert np.isfinite(data.pressure_gpa).all()
            assert (data.volume_a3 > 0).all()
            assert not data.pressure_gpa.flags.writeable


@pytest.mark.parametrize(
    "material,record,count,fit_count",
    [
        ("nacl_b1", "nacl_b1_shen_2026_vinet_1", 67, 55),
        ("nacl_b2", "nacl_b2_shen_2026_vinet_2", 131, 125),
    ],
)
def test_shen_nacl_uses_cu_with_exact_run_pairs_and_unknown_sample_temperature(
    material, record, count, fit_count
):
    document = get_material_document(material)
    original = copy.deepcopy(document)
    data = resolve_dataset_observations(
        document, record, "shen_smith_2026_table_s1_simultaneous_volumes"
    )
    assert len(data.row_indices) == count
    assert data.selected_for_fit.sum() == fit_count
    assert np.isnan(data.temperature_k).all()
    assert np.isnan(data.pressure_uncertainty).all()
    assert data.volume_uncertainty_role == "standard_error"
    assert data.provenance["calibration"]["temperature_k"] == 298
    for i, j in zip(data.row_indices, data.calibrant_row_indices):
        assert data.dataset["experiment"][i] == data.dataset["experiment"][j]
        assert data.dataset["run_number"][i] == data.dataset["run_number"][j]
        assert data.dataset["phase"][j] == "Cu"
    if material == "nacl_b2":
        # Independent numerical checkpoint retained from the full Cu-record replay.
        np.testing.assert_allclose(
            data.pressure_gpa[:3], [32.93748324, 33.32140952, 33.66751048], atol=5e-9
        )
    assert document == original


def fixture():
    document = get_material_document("coesite")
    record = document["eos_records"][0]["identifier"]
    dataset = {
        "identifier": "fixture",
        "kind": "pressure_volume",
        "reference": "isolated fixture",
        "source_location": "fixture",
        "used_by_eos_records": [record],
        "columns": [
            {"name": "p", "quantity": "pressure", "unit": "kbar", "role": "value"},
            {
                "name": "v",
                "quantity": "unit_cell_volume",
                "unit": CELL_UNIT,
                "role": "value",
            },
            {"name": "t", "quantity": "temperature", "unit": "degC", "role": "value"},
            {
                "name": "dv",
                "quantity": "unit_cell_volume",
                "unit": CELL_UNIT,
                "role": "standard_error",
                "of": "v",
            },
            {
                "name": "branch",
                "quantity": "phase_identifier",
                "unit": "dimensionless",
                "role": "flag",
            },
        ],
        "rows": [
            [10, 100, None, 0.2, 1],
            [20, 90, 20, 0.3, 1],
            [30, 80, 30, 0.4, 2],
            [None, 70, 40, 0.5, 1],
        ],
    }
    view = {
        "identifier": "fixture_view",
        "eos_records": [record],
        "evidence_kind": "experimental",
        "notes": "The source explicitly identifies branch 1.",
        "phase_reviewed": True,
        "row_selection": [{"column": "branch", "values": [1]}],
        "coordinates": {
            "pressure": {"kind": "column", "column": "p"},
            "volume": {"kind": "column", "column": "v", "unit": CELL_UNIT},
            "temperature": {"column": "t"},
        },
        "uncertainties": {"volume": "dv"},
    }
    dataset["observation_views"] = [view]
    document["datasets"] = [dataset]
    return document, record, dataset, view


def test_phase_filter_unit_conversion_and_missing_temperature_do_not_drop_pv():
    document, record, _, _ = fixture()
    data = resolve_dataset_observations(
        document, record, "fixture", registry=EMPTY_REGISTRY
    )
    np.testing.assert_equal(data.row_indices, [0, 1])
    np.testing.assert_allclose(data.pressure_gpa, [1, 2])
    np.testing.assert_allclose(data.temperature_k, [np.nan, 293.15])
    np.testing.assert_allclose(data.volume_uncertainty, [0.2, 0.3])
    assert data.volume_uncertainty_role == "standard_error"
    assert data.omissions == ({"row_index": 3, "reason": "invalid_pressure_or_volume"},)


def test_no_fallback_to_selected_eos_v0_or_automatic_phase_selection():
    document, record, dataset, _ = fixture()
    del dataset["observation_views"]
    assert (
        plan_dataset_observations(document, record, "fixture", registry=EMPTY_REGISTRY)[
            "reason"
        ]
        == "sample_or_phase_selection_required"
    )
    dataset["columns"].pop()
    for row in dataset["rows"]:
        row.pop()
    dataset["columns"][1].update(quantity="volume_ratio", unit="dimensionless")
    assert (
        plan_dataset_observations(document, record, "fixture", registry=EMPTY_REGISTRY)[
            "reason"
        ]
        == "coordinate_binding_required"
    )


def test_stale_source_or_material_binding_fails_closed():
    document = get_material_document("nacl_b1")
    document["formula_units_per_cell"] = 8
    with pytest.raises(DatasetError, match="source binding changed"):
        resolve_dataset_observations(
            document,
            "nacl_b1_shen_2026_vinet_1",
            "shen_smith_2026_table_s1_simultaneous_volumes",
        )


@pytest.mark.parametrize(
    "mutation",
    [
        lambda v: v.update(evidence_kind="observations"),
        lambda v: v["coordinates"]["volume"].update(factor=0),
        lambda v: v["coordinates"]["volume"].update(
            factor=4
        ),  # Missing conversion reference.
        lambda v: v["coordinates"]["pressure"].update(column="missing"),
        lambda v: v["uncertainties"].update(pressure="dv"),
    ],
)
def test_embedded_schema_validator_rejects_bad_bindings(mutation):
    document, _, _, view = fixture()
    mutation(view)
    with pytest.raises(EosmatError):
        validate_eosmat_document(document)


def test_lattice_and_reciprocal_density_conversion_preserve_error_meaning():
    document, record, dataset, view = fixture()
    volume = view["coordinates"]["volume"]
    volume.update(
        factor=1000, reciprocal=True, reference="Explicit source mass/basis conversion"
    )
    data = resolve_dataset_observations(
        document, record, "fixture", registry=EMPTY_REGISTRY
    )
    np.testing.assert_allclose(data.volume_a3, [10, 1000 / 90])
    np.testing.assert_allclose(data.volume_uncertainty, [0.02, 1000 * 0.3 / 90**2])
    dataset["columns"][1].update(quantity="lattice_parameter_a", unit="angstrom")
    del view["uncertainties"]
    view["coordinates"]["volume"] = {
        "kind": "lattice",
        "unit": CELL_UNIT,
        "columns": {"a": "v"},
        "constraints": {"b": "a", "c": "a", "alpha": 90, "beta": 90, "gamma": 90},
        "reference": "Explicit cubic conventional setting",
    }
    data = resolve_dataset_observations(
        document, record, "fixture", registry=EMPTY_REGISTRY
    )
    np.testing.assert_allclose(data.volume_a3, [100**3, 90**3])
    assert np.isnan(data.volume_uncertainty).all()


def test_duplicate_calibrants_are_omitted_without_averaging():
    document = get_material_document("nacl_b1")
    registry = observation_registry()
    entry = next(
        e
        for e in registry["datasets"]
        if e["material_id"] == "nacl_b1" and "shen" in e["dataset_id"]
    )
    dataset = next(
        d for d in document["datasets"] if d["identifier"] == entry["dataset_id"]
    )
    view = entry["views"][0]
    # Long-table join on one deliberately nonunique key is still rejected per row.
    view["coordinates"]["pressure"]["pair_by"] = ["experiment"]
    dataset["observation_views"] = [view]
    data = resolve_dataset_observations(
        document, view["eos_records"][0], dataset["identifier"], registry=EMPTY_REGISTRY
    )
    assert not len(data.row_indices)
    assert len(data.omissions) == 67
    assert {x["reason"] for x in data.omissions} == {"ambiguous_calibrant"}


def test_immutable_export_includes_embedded_views_without_a_registry_entry(
    tmp_path, monkeypatch
):
    document, _, _, _ = fixture()
    material = tmp_path / "peritheos/data/materials/coesite.eosmat"
    material.parent.mkdir(parents=True)
    material.write_text(json.dumps(document))
    registry = tmp_path / "peritheos/data/observation-views.json"
    registry.write_text(
        json.dumps(
            {
                "format": "peritheos.observation-views",
                "format_version": 1,
                "datasets": [],
            }
        )
    )
    subprocess.run(["git", "init", "-q", str(tmp_path)], check=True)
    subprocess.run(["git", "-C", str(tmp_path), "add", "."], check=True)
    subprocess.run(
        [
            "git",
            "-C",
            str(tmp_path),
            "-c",
            "user.name=Fixture",
            "-c",
            "user.email=fixture@example.invalid",
            "commit",
            "-qm",
            "isolated observation fixture",
        ],
        check=True,
    )
    spec = importlib.util.spec_from_file_location(
        "observation_export_fixture", ROOT / "scripts/export_observation_views.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    monkeypatch.setattr(module, "ROOT", tmp_path)
    result = module.export_observation_views("HEAD")
    assert len(result["datasets"]) == 1
    assert (
        result["datasets"][0]["views"] == document["datasets"][0]["observation_views"]
    )
    assert not result["incompatible_datasets"]
    assert result["view_registry_origin"].startswith("immutable Peritheos commit")
