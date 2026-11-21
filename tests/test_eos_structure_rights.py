"""Component assessments must preserve science and never become blanket grants."""

import csv
import hashlib
import json
import shutil
from collections import Counter

import pytest

from peritheos_build import ROOT, digest, source_document, stage_core, stage_datasets
from scripts.inventory_distribution_licenses import components
from scripts.review_eos_structure_rights import (
    POLICY,
    REVIEW,
    build_review,
)
from scripts.scientific_source_snapshot import scientific_sha256
from scripts.split_eosmat_audits import pre_split_view
from tests.catalog_state import (
    added_records,
    expected,
    material_count,
    record_count,
)


def test_complete_review_is_source_bound_and_scientific_content_is_unchanged():
    outputs = build_review()
    for path, content in outputs.items():
        assert (ROOT / path).read_bytes() == content, path
    review = json.loads(outputs[REVIEW])
    assert scientific_sha256() == review["catalog_scientific_sha256"]
    assert review["counts"]["materials"] == material_count()
    assert review["counts"]["eos_records"] == record_count()
    # The review of the catalog curated before the catalog revision ledger;
    # each added record carries its own rights statement in the ledger.
    added = {
        e["key"]
        for e in review["entries"]
        if e["kind"] == "eos" and e["identifier"] in added_records()
    }
    entries = [e for e in review["entries"] if e["key"] not in added]
    groups = {
        e["reference"].get("doi") or digest(e["reference"])
        for e in entries
        if e["kind"] == "eos"
    }
    assert {
        "materials": len({e["material"] for e in entries}),
        "eos_records": sum(e["kind"] == "eos" for e in entries),
        "structural_components": sum(
            e["kind"] == "structure_and_diffraction" for e in entries
        ),
        "eos_source_groups": len(groups),
        "migration_records": sum(bool(e.get("migration_source")) for e in entries),
        "assessments": dict(sorted(Counter(e["assessment"] for e in entries).items())),
    } == {
        "materials": 235,
        "eos_records": 687,
        "structural_components": 235,
        "eos_source_groups": 273,
        "migration_records": 144,
        "assessments": {
            "cod_source_and_numerical_changes": 5,
            "independently_calculated_pattern": 10,
            "no_structural_values": 36,
            "numerical_parameter_assessment": 680,
            "numerical_structure_assessment": 184,
            "thermodynamic_collection_review_required": 7,
        },
    }
    # With revisions: base plus the source groups the added records open.
    new_groups = {
        e["reference"].get("doi") or digest(e["reference"])
        for e in review["entries"]
        if e["key"] in added
    } - groups
    want = expected()
    assert review["counts"]["eos_source_groups"] == (
        want["eos_source_groups"] + len(new_groups)
    )
    assert (
        review["counts"]["assessments"]["numerical_parameter_assessment"]
        == (want["numerical_parameter_assessment"])
    )
    by_key = {e["key"]: e for e in review["entries"]}
    # An unsuccessful ICSD search is not a redistributed ICSD structure.
    assert by_key["structure/coesite_v"]["assessment"] == "no_structural_values"
    for material in (
        "hydrous_stishovite_nisr_2017",
        "hydrous_silica_cacl2_nisr_2017",
        "phase_e_fe0072",
    ):
        assert (
            by_key[f"structure/{material}"]["assessment"]
            == "numerical_structure_assessment"
        )
    assert not any(e.get("legacy_source_comparison") for e in review["entries"])
    assert (
        sum(
            e["peak_count"]
            for e in review["entries"]
            if e["assessment"] == "independently_calculated_pattern"
        )
        == 1154
    )


def test_scientific_or_provenance_changes_require_review_not_automatic_approval(
    tmp_path,
):
    (tmp_path / "docs/data").mkdir(parents=True)
    shutil.copy2(ROOT / POLICY, tmp_path / POLICY)
    shutil.copytree(
        ROOT / "peritheos/data/materials", tmp_path / "peritheos/data/materials"
    )
    path = tmp_path / "peritheos/data/materials/mgsio3_post_perovskite.eosmat"
    original = path.read_bytes()
    for mutate in ("coefficient", "provenance"):
        document = json.loads(original)
        if mutate == "coefficient":
            document["eos_records"][0]["parameters"]["K0"]["value"] += 1
        else:
            document["notes"] += " Unverified claim of open permission."
        path.write_text(json.dumps(document))
        with pytest.raises(ValueError, match="Reviewed material changed"):
            build_review(tmp_path)


def test_current_assessments_do_not_relicense_sources_or_backdate_release_rows():
    path = ROOT / "docs/data/license-inventory/catalog-components.csv"
    with path.open(encoding="utf-8", newline="") as stream:
        rows = list(csv.DictReader(stream))
    for row in rows:
        if row["state"] == "release" or row["kind"] == "dataset":
            assert not row["redistribution_assessment"]
            assert not row["review_evidence"]
        elif row["kind"] == "eos":
            assert row["status"] == "source_terms_not_recorded"
            assert row["redistribution_assessment"]
            assert row["review_evidence"]
    # No matching DOI is used as a grant to the EOS or its structural model.
    sakai = next(
        r
        for r in rows
        if r["state"] == "local"
        and r["identifier"] == "gold_sakai_2025_rydberg_stacey_1"
    )
    assert sakai["status"] == "source_terms_not_recorded"
    assert sakai["redistribution_assessment"] == "numerical_parameter_assessment"
    assert not sakai["remaining_issue"]
    replacement = next(
        r
        for r in rows
        if r["state"] == "local"
        and r["kind"] == "structure_and_diffraction"
        and r["material"] == "mgsio3_post_perovskite"
    )
    assert not replacement["remaining_issue"]
    assert replacement["status"] == "open_source_declared"
    assert (
        replacement["redistribution_assessment"] == "cod_source_and_numerical_changes"
    )
    # Dedicated structural sources are not silently omitted by the inventory.
    super_b = next(
        r
        for r in rows
        if r["state"] == "local"
        and r["kind"] == "structure_and_diffraction"
        and r["material"] == "superhydrous_phase_b_lt"
    )
    assert "structure_provenance" in json.loads(super_b["source"])


def test_stale_component_receipt_cannot_be_applied_to_changed_sources():
    path = "peritheos/data/materials/gold.eosmat"
    content = (ROOT / path).read_bytes() + b"\n"
    review = json.loads((ROOT / REVIEW).read_bytes())
    with pytest.raises(ValueError, match="EOS/structure review source changed"):
        components({path: content}, [], True, review)


def test_both_exports_preserve_upstream_notices_review_and_runtime_science(tmp_path):
    review = json.loads((ROOT / REVIEW).read_bytes())
    core = stage_core(tmp_path / "core")
    optional = stage_datasets(tmp_path / "optional")
    core_data = core / "peritheos/data"
    optional_data = optional / "peritheos_datasets"
    for base in (core_data, optional_data):
        assert (base / "eos-structure-rights-review.json").read_bytes() == (
            ROOT / REVIEW
        ).read_bytes()
    for entry in review["entries"]:
        # Staged files are split like the repository; compare the reviewed
        # pre-split runtime content.
        document = pre_split_view(
            source_document(
                core_data / "materials" / f"{entry['material']}.eosmat", core
            )
        )
        if entry["kind"] == "eos":
            record = next(
                r
                for r in document["eos_records"]
                if r["identifier"] == entry["identifier"]
            )
            assert digest(record) == entry["runtime_sha256"]
    for name in (
        "dioptas-materials.LICENSE.md",
        "eos-structure-compilations.LICENSE.md",
    ):
        original = (ROOT / "peritheos/data/datasets" / name).read_bytes()
        assert (core_data / "datasets" / name).read_bytes() == original
        assert (optional_data / "notices" / name).read_bytes() == original
    notice = (core_data / "datasets/dioptas-materials.LICENSE.md").read_bytes()
    license_bytes = notice[notice.index(b"MIT License\n") :]
    assert (
        hashlib.sha256(license_bytes).hexdigest() == review["dioptas"]["license_sha256"]
    )
    # Review metadata must not become a route for source artifacts or data values.
    assert (
        b"/Users/" not in (core_data / "eos-structure-rights-review.json").read_bytes()
    )
    assert not list(core_data.rglob("*.csv"))
    manifest = json.loads((optional_data / "manifest.json").read_bytes())
    assert manifest["total_dataset_count"] == 467
    assert len(manifest["resources"]) == 406
    assert not any(
        "campbell-2009-tables-s4-s5" in path
        or "sakai-2025-table-s8-calibration-lineage" in path
        for path in manifest["resources"]
    )
