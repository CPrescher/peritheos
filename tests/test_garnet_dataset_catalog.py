"""Recovered inputs stay usable without inventing source errors or fit parity."""

import copy
import csv
import hashlib
import io
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_material
from peritheos.datasets import Dataset
from peritheos.errors import DatasetError
from peritheos_build import stage_core, stage_datasets
from scripts.garnet_dataset_bindings import RECEIPT, prior_dataset_bytes
from scripts.garnet_dataset_catalog import AUDIT, PYROPE_IDS, ROOT, tables
from scripts.scientific_source_snapshot import (
    assert_snapshot_lineage,
    scientific_snapshot,
)


def test_each_declaration_and_resource_matches_the_recovered_source():
    for material, entries in tables().items():
        actual = get_material(material)
        for declaration, payload, _ in entries:
            dataset = actual.get_dataset(declaration["identifier"])
            assert dataset.checksum_verified
            assert (
                len(dataset) == len(list(csv.reader(payload.decode().splitlines()))) - 1
            )
            assert all(not array.flags.writeable for array in dataset.data.values())
            path = ROOT / "peritheos/data" / declaration["resource"]["path"]
            assert path.read_bytes() == payload


def test_author_input_tokens_errors_order_and_scale_groups_are_preserved():
    original = list(
        csv.DictReader(
            io.StringIO((AUDIT / "pyrope/author-2026-11-14-inputs.csv").read_text())
        )
    )
    material = get_material("pyrope")
    for kind, identifier, value, error in [
        (
            "V",
            PYROPE_IDS[0],
            "molar_volume_cm3_mol",
            "molar_volume_assigned_error_cm3_mol",
        ),
        (
            "KS",
            PYROPE_IDS[1],
            "bulk_modulus_adiabatic_gpa",
            "bulk_modulus_assigned_error_gpa",
        ),
    ]:
        dataset = material.get_dataset(identifier)
        source = [row for row in original if row["kind"] == kind]
        assert len(source) == len(dataset) == (115 if kind == "V" else 79)
        for key, name in [
            ("T", "temperature_k"),
            ("sT", "temperature_assigned_error_k"),
            ("P", "pressure_gpa"),
            ("sP", "pressure_assigned_error_gpa"),
            ("Y", value),
            ("sY", error),
            ("group", "group"),
            ("fit_row", "fit_row"),
        ]:
            np.testing.assert_array_equal(
                dataset[name], [float(row[key]) for row in source]
            )
        rows = list(
            csv.DictReader(
                io.StringIO(
                    (ROOT / "peritheos/data" / dataset.resource.path).read_text()
                )
            )
        )
        for row, old in zip(rows, source):
            for key, name in [
                ("T", "temperature_k"),
                ("sT", "temperature_assigned_error_k"),
                ("P", "pressure_gpa"),
                ("sP", "pressure_assigned_error_gpa"),
                ("Y", value),
                ("sY", error),
            ]:
                assert row[name] == old[key + "_token"]
        assert dataset.get_column(error).role == "uncertainty"
    modulus = material.get_dataset(PYROPE_IDS[1])
    selected = modulus["group"] == 3
    assert selected.sum() == 15
    assert set(modulus["source_dataset"][selected]) == {"sinogeikin_2000_ks"}
    ambient = (modulus["source_dataset"] == "sinogeikin_2000_ks") & (
        modulus["pressure_gpa"] < 0.001
    )
    assert ambient.sum() == 2 and np.all(modulus["group"][ambient] == 2)
    pv = material.get_dataset(PYROPE_IDS[0]).as_pressure_volume()
    assert pv.volume_unit == "cm^3/mol" and pv.pressure_uncertainty is not None
    assert (
        pv.pressure_sigma is None
    )  # Assigned weights are not relabelled as measured SDs.


def test_grossular_missing_errors_specimen_and_digitization_remain_explicit():
    material = get_material("grossular")
    isaak = material.get_dataset("grossular_isaak_1992_table3b_ggd_moduli")
    assert np.isnan(isaak["ks_gpa_reported_error"]).sum() == 19
    assert set(
        isaak["ks_gpa_reported_error"][np.isfinite(isaak["ks_gpa_reported_error"])]
    ) == {0.71, 1.02, 1.20}
    strain = material.get_dataset("grossular_isaak_1992_table2_gr2b_strain")
    assert all("76.1" in value for value in strain["observed_specimen"])
    assert np.isnan(strain["source_error_volume_strain"]).all()
    assert strain["observed_volume_strain"][0] == pytest.approx(0.013 / 100)
    assert strain.used_by_eos_records == ()
    conrad = material.get_dataset("grossular_conrad_1999_table6_stiffness")
    assert (
        len(conrad) == 5 and (conrad["selection_angel_2022"] == "selected").sum() == 3
    )
    assert "ks_gpa" not in conrad.data
    gwanmesia = material.get_dataset("grossular_gwanmesia_2014_table1")
    assert len(gwanmesia) == 52
    assert gwanmesia.values("temperature_c", unit="K")[0] == pytest.approx(303.15)
    assert not any(
        column.role in {"uncertainty", "standard_deviation"}
        for column in gwanmesia.columns
    )
    assert (
        gwanmesia.metadata["uncertainty"]["reported_relative_error_estimates_percent"][
            "ks"
        ]
        == 0.8
    )
    figure = material.get_dataset("grossular_thieblot_1998_angel_2022_figure4")
    assert (figure["selection"] == "selected").sum() == 7
    assert (figure["selection"] == "excluded_visible").sum() == 6
    assert not {"x_pdf_point", "y_pdf_point"} & set(figure.data)
    assert all(column.role != "uncertainty" for column in figure.columns)


@pytest.mark.parametrize("name", ["pyrope", "grossular"])
def test_published_records_and_frozen_predecessors_are_preserved(name):
    receipt = json.loads((ROOT / RECEIPT).read_bytes())
    path = ROOT / f"peritheos/data/materials/{name}.eosmat"
    previous = prior_dataset_bytes(path.read_bytes(), name, ROOT)
    assert (
        hashlib.sha256(previous).hexdigest()
        == receipt["materials"][name]["before_sha256"]
    )
    old = json.loads(previous)
    current = json.loads(path.read_bytes())
    for a, b in zip(old["eos_records"], current["eos_records"]):
        for key in (
            "model",
            "parameters",
            "thermal",
            "parameter_covariance",
            "default",
            "reference",
            "validity",
            "temperature_ref",
        ):
            assert a.get(key) == b.get(key)
        old_validation = copy.deepcopy(a["evidence"]["scientific_validation"])
        new_validation = copy.deepcopy(b["evidence"]["scientific_validation"])
        if name == "grossular":
            old_check = old_validation["primary_data_check"]
            new_check = new_validation["primary_data_check"]
            assert new_check["scope"] == "original_author_processed_joint_fit_inputs"
            assert old_check["status"] == new_check["status"] == "not_available"
            for field in ("finding", "audit_date", "scope"):
                old_check.pop(field, None)
                new_check.pop(field, None)
        assert old_validation == new_validation
    wrong = copy.deepcopy(current)
    wrong["eos_records"][0]["parameters"]["K0"]["value"] += 1
    with pytest.raises(ValueError, match="source changed"):
        prior_dataset_bytes(json.dumps(wrong).encode(), name, ROOT)
    assert_snapshot_lineage(scientific_snapshot(), receipt["scientific_sha256_before"])


@pytest.mark.parametrize("version", [3, 4, 5])
def test_dataset_metadata_and_missing_errors_survive_supported_exports(version):
    for name, entries in tables().items():
        material = get_material(name)
        exported = material.to_eosmat(format_version=version)
        Draft202012Validator(eosmat_schema(version)).validate(exported)
        reloaded = Material.from_eosmat(exported)
        for ds, _, _ in entries:
            a = material.get_dataset(ds["identifier"])
            b = reloaded.get_dataset(ds["identifier"])
            assert [vars(c) for c in a.columns] == [vars(c) for c in b.columns]
            for key in a.data:
                if a.get_column(key).unit in {"text", "identifier"}:
                    assert list(a[key]) == list(b[key])
                else:
                    np.testing.assert_array_equal(a[key], b[key])


def test_new_resources_reject_byte_tampering(tmp_path):
    ds, payload, _ = tables()["pyrope"][0]
    path = tmp_path / ds["resource"]["path"]
    path.parent.mkdir(parents=True)
    path.write_bytes(payload.replace(b"111.061569", b"111.061568", 1))
    with pytest.raises(DatasetError, match="SHA-256 mismatch"):
        Dataset.from_mapping(ds, resource_root=tmp_path)


def test_exported_installed_layout_loads_new_tables_without_repository_files(tmp_path):
    core = stage_core(tmp_path / "core")
    optional = stage_datasets(tmp_path / "optional")
    from peritheos import _rust

    (core / "peritheos" / Path(_rust.__file__).name).symlink_to(
        Path(_rust.__file__).resolve()
    )
    env = {**os.environ, "PYTHONPATH": os.pathsep.join([str(core), str(optional)])}
    program = """
from peritheos import get_material
p = get_material("pyrope")
g = get_material("grossular")
assert len(p.get_dataset("pyrope_angel_2022_author_pvt")) == 115
assert len(p.get_dataset("pyrope_angel_2022_author_pkt")) == 79
assert len(g.get_dataset("grossular_gwanmesia_2014_table1")) == 52
assert len(g.datasets) == 12
print("Installed-layout garnet data access passed")
"""
    result = subprocess.run(
        [sys.executable, "-c", program],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=True,
    )
    assert "Installed-layout garnet data access passed" in result.stdout
