"""Distribution boundaries preserve models without redistributing observations."""

import copy
import json
from types import SimpleNamespace

import numpy as np
import pytest

from peritheos import Material, get_material_document, validate_eosmat_document
from peritheos.errors import DatasetError
from peritheos_build import (
    ROOT,
    digest,
    runtime_document,
    stage_core,
    stage_datasets,
    version,
)


def test_every_runtime_export_preserves_models_structures_and_scientific_status():
    for path in (ROOT / "peritheos/data/materials").glob("*.eosmat"):
        original = json.loads(path.read_bytes())
        exported = runtime_document(original)
        validate_eosmat_document(exported)
        assert all("rows" not in d for d in exported.get("datasets", []))
        for key in ("formula", "lattice", "atom_sites", "space_group_number", "peaks"):
            assert exported.get(key) == original.get(key), (path.name, key)
        for before, after in zip(original["eos_records"], exported["eos_records"]):
            for key in (
                "identifier",
                "model",
                "parameters",
                "thermal",
                "parameter_covariance",
                "reference",
                "validity",
                "pressure_calibration",
                "record_kind",
                "default",
                "default_for",
                "catalog_access",
                "experimental_pressure_range_gpa",
                "experimental_temperature_range_k",
            ):
                assert after.get(key) == before.get(key), (path.name, key)
            for key in ("status", "note", "usage_recommendation"):
                assert after["evidence"]["scientific_validation"].get(key) == before[
                    "evidence"
                ]["scientific_validation"].get(key)
        assert exported["data_availability"]["source_document_sha256"] == digest(
            original
        )


def test_core_export_contains_no_observation_resources_or_research(tmp_path):
    stage = stage_core(tmp_path / "core")
    assert not (stage / "docs/data").exists()
    assert not (stage / "peritheos/data/research").exists()
    assert not list((stage / "peritheos/data").rglob("*.csv"))
    assert len(list((stage / "peritheos/data/datasets").glob("*.cif"))) == 3
    assert not (stage / "peritheos/data/primary-source-audit.json").exists()
    for path in (stage / "peritheos/data/materials").glob("*.eosmat"):
        document = json.loads(path.read_bytes())
        assert all("rows" not in d for d in document.get("datasets", []))
    assert 'build-backend = "maturin"' in (stage / "pyproject.toml").read_text()


def test_optional_export_is_exactly_allowlisted_and_source_preserving(tmp_path):
    stage = stage_datasets(tmp_path / "optional")
    package = stage / "peritheos_datasets"
    manifest = json.loads((package / "manifest.json").read_bytes())
    assert manifest["dataset_count"] == 44
    expected_files = set()
    for entry in manifest["policy"]["datasets"]:
        assert entry["source_status"] in {
            "open_license_declared",
            "open_source_declared",
            "author_permission_declared",
        }
        material = json.loads(
            (
                ROOT / "peritheos/data/materials" / f"{entry['material']}.eosmat"
            ).read_bytes()
        )
        expected = next(
            d for d in material["datasets"] if d["identifier"] == entry["dataset"]
        )
        exported = json.loads(
            (package / "declarations" / f"{entry['material']}.json").read_bytes()
        )
        assert expected == next(
            d for d in exported if d["identifier"] == entry["dataset"]
        )
        if "resource" in expected:
            name = expected["resource"]["path"]
            expected_files.add(name)
            assert (package / "data" / name).read_bytes() == (
                ROOT / "peritheos/data" / name
            ).read_bytes()
    actual_files = {
        p.relative_to(package / "data").as_posix()
        for p in (package / "data").rglob("*")
        if p.is_file()
    }
    assert actual_files == expected_files


def test_optional_export_refuses_declaration_drift(tmp_path, monkeypatch):
    import peritheos_build

    real = peritheos_build.digest
    monkeypatch.setattr(
        peritheos_build,
        "digest",
        lambda value: "bad" if "columns" in value else real(value),
    )
    with pytest.raises(ValueError, match="declaration changed"):
        stage_datasets(tmp_path / "optional")


def test_runtime_eos_works_without_data_and_missing_data_has_actionable_error(
    tmp_path, monkeypatch
):
    import peritheos.datasets

    original_files = peritheos.datasets.resources.files
    monkeypatch.setattr(
        peritheos.datasets.resources,
        "files",
        lambda package: (
            tmp_path if package == "peritheos.data" else original_files(package)
        ),
    )
    full = get_material_document("brucite", format_version=4)
    exported = runtime_document(full)
    material = Material.from_eosmat(exported)
    before = Material.from_eosmat(full)
    volumes = np.array([40.0, 38.0, 35.0])
    np.testing.assert_array_equal(
        material.default_record().pressure(volumes, 300),
        before.default_record().pressure(volumes, 300),
    )
    with pytest.raises(DatasetError, match=r"peritheos\[datasets\]") as error:
        material.get_dataset("brucite_fukui_2003_table4")
    assert error.value.code == "dataset.package_not_installed"


def test_user_authored_standalone_datasets_are_not_changed_by_package_split():
    original = get_material_document("brucite", format_version=4)
    saved = copy.deepcopy(original)
    material = Material.from_eosmat(original)
    assert len(material.get_dataset("brucite_fukui_2003_table4")) == 37
    assert original == saved


def test_optional_install_hydrates_inline_and_csv_data_without_changing_models(
    tmp_path, monkeypatch
):
    import peritheos._optional_data as optional
    import peritheos.datasets

    core = stage_core(tmp_path / "core")
    data = stage_datasets(tmp_path / "optional") / "peritheos_datasets"
    original_files = peritheos.datasets.resources.files

    def files(package):
        if package == "peritheos.data":
            return core / "peritheos/data"
        if package == "peritheos.data.materials":
            return core / "peritheos/data/materials"
        return original_files(package)

    monkeypatch.setattr(peritheos.datasets.resources, "files", files)
    monkeypatch.setattr(optional, "data_package", lambda: data)
    for name, identifier in [
        ("b4c", "b4c_somayazulu_2023_table4_pvt"),
        ("mgo", "zhu_2025_mgo_pvt"),
    ]:
        document = get_material_document(name)
        declaration = next(
            d for d in document["datasets"] if d["identifier"] == identifier
        )
        assert "distribution" not in declaration
        assert len(peritheos.datasets.Dataset.from_mapping(declaration)) > 0
    # Version-compatible installation does not imply every research dataset is admitted.
    document = get_material_document("brucite")
    with pytest.raises(DatasetError) as error:
        Material.from_eosmat(document).get_dataset("brucite_fukui_2003_table4")
    assert error.value.code == "dataset.not_distributed"
    # Corrupt optional bytes fail the same checksum verification as core resources.
    declaration = next(
        d
        for d in get_material_document("mgo")["datasets"]
        if d["identifier"] == "zhu_2025_mgo_pvt"
    )
    (data / "data" / declaration["resource"]["path"]).write_text("corrupt")
    with pytest.raises(DatasetError) as error:
        peritheos.datasets.Dataset.from_mapping(declaration)
    assert error.value.code == "dataset.checksum_mismatch"


def test_mismatched_optional_version_does_not_block_eos_document_loading(monkeypatch):
    import peritheos._optional_data as optional

    monkeypatch.setattr(optional.importlib.util, "find_spec", lambda name: object())
    monkeypatch.setattr(
        optional, "import_module", lambda name: SimpleNamespace(__version__="0.0.0")
    )
    assert optional.dataset_declarations("mgo") == []
    with pytest.raises(DatasetError) as error:
        optional.data_package()
    assert error.value.code == "dataset.package_version_mismatch"


def test_optional_dependency_and_data_project_versions_are_release_consistent():
    release = version()
    assert (
        f'datasets = ["peritheos-datasets=={release}"]'
        in (ROOT / "pyproject.toml").read_text()
    )
    project = (ROOT / "packages/peritheos-datasets/pyproject.toml").read_text()
    assert f'version = "{release}"' in project
    assert f'dependencies = ["peritheos=={release}"]' in project


def test_optional_export_refuses_a_restricted_source_scope(tmp_path):
    root = tmp_path / "source"
    (root / "peritheos").mkdir(parents=True)
    (root / "peritheos/__init__.py").write_text('__version__ = "0.17.4"\n')
    policy = root / "packages/peritheos-datasets/datasets.json"
    policy.parent.mkdir(parents=True)
    policy.write_text(
        json.dumps(
            {"datasets": [{"source_status": "restricted_source", "dataset": "blocked"}]}
        )
    )
    with pytest.raises(ValueError, match="unapproved source scope"):
        stage_datasets(tmp_path / "optional", root=root)
