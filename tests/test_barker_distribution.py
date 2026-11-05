"""Supporting-study admission preserves source meaning and package boundaries."""

import copy
import hashlib
import json

import numpy as np
import pytest

from peritheos import (
    DatasetError,
    DatasetLookupError,
    get_study_dataset,
    get_study_document,
    list_studies,
)
from peritheos_build import ROOT, digest, stage_core, stage_datasets
from scripts.scientific_source_snapshot import scientific_sha256

STUDY = "argon_barker_1987"
LIQUID = "argon_barker_1987_table1_liquid"
FIGURE = "argon_barker_1987_figure5_selected"
PROOF_PATH = ROOT / "docs/data/barker-1987-redistribution.json"
PROOF = json.loads(PROOF_PATH.read_bytes())
SOURCE = ROOT / "peritheos/data/studies/argon-barker-1987.json"


def test_original_study_and_all_numerical_evidence_survive_additive_metadata():
    current = json.loads(SOURCE.read_bytes())
    assert (
        hashlib.sha256(SOURCE.read_bytes()).hexdigest() == PROOF["current_study_sha256"]
    )
    restored = copy.deepcopy(current)
    for key in PROOF["added_study_fields"]:
        restored.pop(key)
    for dataset in restored["datasets"]:
        original_keys = PROOF["original_dataset_fields"][dataset["identifier"]]
        for key in list(dataset):
            if key not in original_keys:
                dataset.pop(key)
    assert set(restored) == set(PROOF["original_study_fields"])
    assert (
        hashlib.sha256((json.dumps(restored, indent=2) + "\n").encode()).hexdigest()
        == PROOF["original_study_sha256"]
    )
    assert scientific_sha256() == PROOF["catalog_scientific_sha256"]
    assert current["phase"] is None
    assert current["reproduction"]["status"] == "not_reproduced"
    assert current["reproduction"]["eos_record_identifiers"] == []
    for dataset in current["datasets"]:
        key = f"study/{STUDY}/{dataset['identifier']}"
        assert digest(dataset) == PROOF["declaration_bindings"][key]
        relative = dataset["resource"]["path"]
        assert (
            hashlib.sha256(
                (ROOT / "peritheos/data" / relative).read_bytes()
            ).hexdigest()
            == PROOF["resource_bindings"][relative]
        )


def test_public_study_api_preserves_phase_theory_and_graphical_bounds():
    assert STUDY in list_studies()
    liquid = get_study_dataset(STUDY, LIQUID)
    figure = get_study_dataset(STUDY, FIGURE)
    assert len(liquid) == 7
    assert liquid.metadata["phase"] == "liquid"
    assert (
        liquid["determination_method"].tolist()
        == ["experimental"] + ["theoretical"] * 6
    )
    np.testing.assert_array_equal(
        liquid["pv_over_nkt"], [2.15, 0.90, 2.09, 0.77, 1.95, 0.71, 1.89]
    )
    assert len(figure) == 32
    assert figure.metadata["phase"] == "solid_unspecified"
    assert sum(figure["determination_method"] == "experimental") == 11
    assert sum(figure["determination_method"] == "theoretical") == 21
    np.testing.assert_allclose(
        figure.values("pressure_kbar", unit="GPa"),
        figure["pressure_gpa"],
        # Both stored pressure columns were independently rounded to six decimals.
        atol=5.5e-7,
        rtol=1e-14,
    )
    assert figure.get_column("digitization_bound_pressure_gpa").role == "bound"
    assert figure.get_column("molar_volume_cm3_mol").unit == "cm^3/mol"
    assert figure.get_column("volume_atomic_a3").unit == "angstrom^3/atom"
    assert liquid.used_by_eos_records == figure.used_by_eos_records == ()
    assert liquid.checksum_verified and figure.checksum_verified
    with pytest.raises(ValueError):
        liquid["pv_over_nkt"][0] = 0
    document = get_study_document(STUDY)
    document["datasets"][0]["resource"]["sha256"] = "0" * 64
    assert get_study_document(STUDY)["datasets"][0]["resource"]["sha256"] != "0" * 64


def test_study_lookup_does_not_interpret_unknown_identifiers_as_paths():
    with pytest.raises(DatasetLookupError) as error:
        get_study_document("../materials/gold")
    assert error.value.code == "dataset.study_not_found"
    with pytest.raises(DatasetLookupError) as error:
        get_study_dataset(STUDY, "missing")
    assert error.value.code == "dataset.not_found"


@pytest.fixture
def exported_study(tmp_path, monkeypatch):
    from peritheos import _optional_data as optional
    from peritheos import datasets, studies

    core = stage_core(tmp_path / "core") / "peritheos/data"
    package = stage_datasets(tmp_path / "optional") / "peritheos_datasets"
    files = studies.resources.files

    def resources(package):
        return core if package == "peritheos.data" else files(package)

    monkeypatch.setattr(datasets.resources, "files", resources)
    monkeypatch.setattr(optional, "data_package", lambda: package_path[0])
    package_path = [None]
    return core, package, package_path


def test_core_metadata_and_optional_hydration_are_separate(exported_study):
    core, package, package_path = exported_study
    assert list_studies() == (STUDY,)
    document = get_study_document(STUDY)
    assert document["phase"] is None
    assert document["reproduction"]["status"] == "not_reproduced"
    assert "source_audit" not in document
    assert all(d["distribution"]["storage"] == "optional" for d in document["datasets"])
    assert not list(core.rglob("*.csv"))
    assert "/Users/" not in json.dumps(document)
    with pytest.raises(DatasetError) as error:
        get_study_dataset(STUDY, LIQUID)
    assert error.value.code == "dataset.package_not_installed"
    package_path[0] = package
    assert len(get_study_dataset(STUDY, LIQUID)) == 7
    assert len(get_study_dataset(STUDY, FIGURE)) == 32
    assert "distribution" not in get_study_document(STUDY)["datasets"][0]
    manifest = json.loads((package / "manifest.json").read_bytes())
    review = json.loads((package / "redistribution-review.json").read_bytes())
    assert manifest["dataset_count"] == 451
    assert manifest["study_dataset_count"] == 2
    assert manifest["total_dataset_count"] == 453
    assert len(manifest["resources"]) == 392
    assert review["barker_review"] == PROOF
    assert all(e["decision"] == "admit" for e in review["supporting_studies"])
    assert len([e for e in review["entries"] if e["decision"] == "hold"]) == 5
    for entry in manifest["policy"]["study_datasets"]:
        assert entry["study"] == STUDY and "material" not in entry
        assessment = next(
            e
            for e in review["supporting_studies"]
            if e["review_key"] == entry["review_key"]
        )
        assert (
            assessment["followup_evidence"]["sha256"]
            == hashlib.sha256(PROOF_PATH.read_bytes()).hexdigest()
        )
    assert (package / "notices/argon-barker-1987.LICENSE.md").read_bytes() == (
        ROOT / "peritheos/data/datasets/argon-barker-1987.LICENSE.md"
    ).read_bytes()
    assert not any(
        p.suffix.lower() in {".pdf", ".png", ".xlsx", ".xls", ".zip"}
        for p in package.rglob("*")
    )


def test_study_corruption_and_version_mismatch_fail_without_hiding_metadata(
    exported_study, monkeypatch
):
    from peritheos import _optional_data as optional

    _, package, package_path = exported_study
    package_path[0] = package
    resource = package / "data/datasets/argon-barker-1987-table1-liquid.csv"
    resource.write_bytes(resource.read_bytes() + b"\n")
    with pytest.raises(DatasetError) as error:
        get_study_dataset(STUDY, LIQUID)
    assert error.value.code == "dataset.checksum_mismatch"

    def mismatch():
        raise DatasetError(
            "Wrong optional version", code="dataset.package_version_mismatch"
        )

    monkeypatch.setattr(optional, "data_package", mismatch)
    assert get_study_document(STUDY)["reproduction"]["status"] == "not_reproduced"
    with pytest.raises(DatasetError) as error:
        get_study_dataset(STUDY, FIGURE)
    assert error.value.code == "dataset.package_version_mismatch"


def test_study_export_rejects_borrowed_material_approval(tmp_path):
    root = tmp_path / "root"
    root.mkdir()
    (root / "peritheos").symlink_to(ROOT / "peritheos", target_is_directory=True)
    (root / "docs/data").mkdir(parents=True)
    (root / "docs/data/numerical-redistribution-review.json").symlink_to(
        ROOT / "docs/data/numerical-redistribution-review.json"
    )
    policy = json.loads(
        (ROOT / "packages/peritheos-datasets/datasets.json").read_bytes()
    )
    policy["study_datasets"][0]["material"] = "argon_fcc"
    destination = root / "packages/peritheos-datasets/datasets.json"
    destination.parent.mkdir(parents=True)
    destination.write_text(json.dumps(policy))
    with pytest.raises(ValueError, match="Unapproved supporting-study scope"):
        stage_datasets(tmp_path / "rejected", root=root)
