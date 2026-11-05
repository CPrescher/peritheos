"""Distribution boundaries preserve models without redistributing observations."""

import copy
import hashlib
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
    assert manifest["dataset_count"] == 451
    review_path = package / "redistribution-review.json"
    review = json.loads(review_path.read_bytes())
    assert (
        hashlib.sha256(review_path.read_bytes()).hexdigest()
        == manifest["policy"]["redistribution_review"]["sha256"]
    )
    reviewed = {e["review_key"]: e for e in review["entries"]}
    originals = {}
    for path in (ROOT / "peritheos/data/materials").glob("*.eosmat"):
        document = json.loads(path.read_bytes())
        originals.update(
            {f"{path.stem}/{d['identifier']}": d for d in document.get("datasets", [])}
        )
    assert set(reviewed) == set(originals)
    for key, entry in reviewed.items():
        assert entry["declaration_sha256"] == digest(originals[key])
        assert entry["resource"] == originals[key].get("resource")
    selected = {
        f"{e['material']}/{e['dataset']}" for e in manifest["policy"]["datasets"]
    }
    assert selected == {
        k for k, e in reviewed.items() if e["decision"] in {"admit", "retain"}
    }
    held_files = {
        e["resource"]["path"]
        for e in review["entries"]
        if e["decision"] == "hold" and e["resource"]
    }
    assert held_files.isdisjoint(manifest["resources"])
    assert (package / "notices/numerical-transcriptions.LICENSE.md").is_file()
    expected_files = set()
    for entry in manifest["policy"]["datasets"]:
        if entry.get("review_key"):
            assessment = reviewed[entry["review_key"]]
            assert assessment["decision"] == "admit"
            assert assessment["basis"] == entry["redistribution_basis"]
            assert assessment["source_status"] == entry["source_status"]
        else:
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
        ("brucite", "brucite_fukui_2003_table4"),
        ("coesite", "coesite_i_iii_bykova_2018_author_fit_export"),
        ("fes_vi", "fes_morard_2026_author_literature_cold_pv"),
        ("nickel", "campbell_2009_independent_buffer_polynomial_grid"),
        ("gold", "gold_sakai_2025_table_s4_au_w_volumes"),
        ("gold", "gold_sakai_2025_independent_calibration_lineage"),
        ("nickel", "campbell_2009_independent_nno_pressure_integral"),
        ("iridium", "iridium_anzellini_2025_tables_s1_s3_pvt"),
        (
            "fe088sio3_bridgmanite",
            "fe088sio3_bridgmanite_ismailova_2016_table_s2_selected_crystallography",
        ),
    ]:
        document = get_material_document(name)
        declaration = next(
            d for d in document["datasets"] if d["identifier"] == identifier
        )
        assert "distribution" not in declaration
        assert len(peritheos.datasets.Dataset.from_mapping(declaration)) > 0
    # Version-compatible installation does not imply every research dataset is admitted.
    document = get_material_document("gold")
    with pytest.raises(DatasetError) as error:
        Material.from_eosmat(document).get_dataset("gold_sakai_2025_table_s8_lineage")
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


@pytest.mark.parametrize(
    "status",
    [
        "restricted_source",
        "source_terms_not_recorded",
        "limited_or_unspecified_source",
        "mixed_author_permission",
    ],
)
def test_optional_export_refuses_an_unreviewed_source_scope(tmp_path, status):
    root = tmp_path / "source"
    (root / "peritheos").mkdir(parents=True)
    (root / "peritheos/__init__.py").write_text('__version__ = "0.17.4"\n')
    policy = root / "packages/peritheos-datasets/datasets.json"
    policy.parent.mkdir(parents=True)
    policy.write_text(
        json.dumps({"datasets": [{"source_status": status, "dataset": "blocked"}]})
    )
    with pytest.raises(ValueError, match="unapproved source scope"):
        stage_datasets(tmp_path / "optional", root=root)


@pytest.mark.parametrize(
    "change",
    ["review_bytes", "missing_assessment", "held", "wrong_resource", "wrong_basis"],
)
def test_numerical_admission_requires_exact_review_binding(tmp_path, change):
    root = tmp_path / "source"
    (root / "peritheos").mkdir(parents=True)
    (root / "peritheos/__init__.py").write_text('__version__ = "0.17.4"\n')
    policy = json.loads(
        (ROOT / "packages/peritheos-datasets/datasets.json").read_bytes()
    )
    entry = next(
        e
        for e in policy["datasets"]
        if e.get("redistribution_basis") == "numerical_facts_assessment"
    )
    policy["datasets"] = [entry]
    receipt = policy["redistribution_review"]
    review = json.loads((ROOT / receipt["path"]).read_bytes())
    assessment = next(
        e for e in review["entries"] if e["review_key"] == entry["review_key"]
    )
    if change == "missing_assessment":
        review["entries"] = []
    elif change == "held":
        assessment["decision"] = "hold"
    elif change == "wrong_resource":
        assessment["resource"]["sha256"] = "bad"
        material_path = (
            root / "peritheos/data/materials" / f"{entry['material']}.eosmat"
        )
        material_path.parent.mkdir(parents=True)
        material_path.write_bytes(
            (ROOT / "peritheos/data/materials" / material_path.name).read_bytes()
        )
    elif change == "wrong_basis":
        entry["redistribution_basis"] = "unreviewed"
    payload = json.dumps(review, ensure_ascii=False, indent=2).encode() + b"\n"
    review_path = root / receipt["path"]
    review_path.parent.mkdir(parents=True)
    review_path.write_bytes(payload)
    if change != "review_bytes":
        receipt["sha256"] = hashlib.sha256(payload).hexdigest()
    else:
        review_path.write_bytes(payload + b" ")
    policy_path = root / "packages/peritheos-datasets/datasets.json"
    policy_path.parent.mkdir(parents=True)
    policy_path.write_text(json.dumps(policy))
    message = {
        "review_bytes": "review changed",
        "wrong_resource": "resource binding changed",
    }.get(change, "unapproved source scope")
    with pytest.raises(ValueError, match=message):
        stage_datasets(tmp_path / "optional", root=root)
