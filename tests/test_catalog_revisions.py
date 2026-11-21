"""Catalog revisions: reversible, receipt-backed changes of bundled content.

Every test that changes data works on a temporary copy of the repository
parts the tool reads and writes, and picks its targets at run time (a record
lacking the field, a free identifier), so it passes in any ledger state. The
records written here are mechanism fixtures, not curated scientific claims.
"""

from __future__ import annotations

import copy
import hashlib
import json
import os
import shutil
from pathlib import Path

import pytest

from peritheos._eosmat_audit import merge_audit
from scripts import apply_catalog_revision as tool
from scripts import catalog_revisions
from scripts.apply_check_text import apply_plan as apply_check_text
from scripts.apply_check_text import plan_changes as plan_check_text
from scripts.apply_check_text import verify as verify_check_text
from scripts.apply_record_text import DEFAULT_REVISION, _git_file
from scripts.apply_record_text import apply_plan as apply_record_text
from scripts.apply_record_text import plan_changes as plan_record_text
from scripts.apply_record_text import verify as verify_record_text
from scripts.catalog_revisions import (
    LEDGER,
    apply_operations,
    check_bindings,
    prior_revision_bytes,
    prior_revision_snapshot,
    read_ledger,
    verify_ledger,
)
from scripts.migrate_eosmat_v5 import prior_material_bytes
from scripts.review_eos_structure_rights import build_review
from scripts.scientific_source_snapshot import (
    assert_snapshot_lineage,
    scientific_snapshot,
)
from scripts.split_eosmat_audits import (
    LAYOUTS,
    RECEIPT,
    pre_split_bytes,
    split_added_materials,
    verify_catalog,
)

ROOT = Path(__file__).resolve().parents[1]
RECEIPT_PATH = RECEIPT.relative_to(ROOT).as_posix()
SIDERITE = "siderite_fe095mn005"
ONE_SIGMA = 0.6826894921370859
GARNET = json.loads(
    (
        ROOT / "docs/data/garnet-angel-2022/dataset-integration-2026-11-15.json"
    ).read_bytes()
)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


def _copy(target: Path) -> Path:
    tool.copy_repository(ROOT, target)
    return target


@pytest.fixture
def repo(tmp_path, monkeypatch):
    """A writable copy of everything a catalog revision reads or writes."""
    root = _copy(tmp_path / "repo")
    # assert_snapshot_lineage reads the revision ledger of this copy.
    monkeypatch.setattr(catalog_revisions, "ROOT", root)
    return root


@pytest.fixture(scope="module")
def pristine(tmp_path_factory):
    """A copy for tests that only plan (and are refused); never written."""
    return _copy(tmp_path_factory.mktemp("pristine") / "repo")


def _tree(root: Path) -> dict[str, bytes]:
    return {k: v.read_bytes() for k, v in tool._tree(root).items()}


def _records(root: Path, material: str) -> list[dict]:
    path = root / f"peritheos/data/materials/{material}.eosmat"
    return json.loads(path.read_bytes())["eos_records"]


def _identifiers(root: Path) -> set[str]:
    return {
        r["identifier"]
        for path in (root / "peritheos/data/materials").glob("*.eosmat")
        for r in json.loads(path.read_bytes())["eos_records"]
    }


def _common(identifier: str, material: str, kind: str) -> dict:
    return {
        "id": identifier,
        "material": material,
        "kind": kind,
        "date": "2026-11-20",
        "decided_by": "Test",
        "reason": "Mechanism test on a temporary copy.",
        "source": {"citation": "Test fixture", "location": "none"},
        "rights": {
            "reviewed_by": "Test",
            "date": "2026-11-20",
            "statement": "Test fixture; no third-party content is added.",
            **({"same_basis": True} if kind == "add_record" else {}),
        },
    }


def add_request(root: Path, revision: str = "t-add") -> dict:
    """A thermal BM3 record for siderite, under a free identifier and label."""
    taken = _identifiers(root)
    number = next(n for n in range(1, 1000) if f"{SIDERITE}_test_{n}" not in taken)
    identifier = f"{SIDERITE}_test_{number}"
    sidecar = json.loads(
        (root / f"peritheos/data/audits/{SIDERITE}.audit.json").read_bytes()
    )
    first = _records(root, SIDERITE)[0]
    record = {
        "identifier": identifier,
        "label": f"Test fixture {number}, thermal BM3",
        "model": "BM3",
        "temperature_ref": 300.0,
        "parameters": {
            "V0": {"value": 293.4, "error": 0.1, "fixed": True},
            "K0": {"value": 120.0, "error": 2.0},
            "K0_prime": {"value": 3.57, "error": 0.09},
        },
        "thermal": {
            "model": "AlphaKT",
            "parameters": {
                "Tr": {"value": 300.0, "error": None, "fixed": True},
                "alpha0": {"value": 3.77e-05, "error": 0.09e-05},
                "alpha1": {"value": 0.06e-08, "error": 0.14e-08},
                "dK_dT": {"value": -0.015, "error": 0.001},
            },
            "thermal_expansion_law": "linear_temperature",
            "reference_volume_law": "integrated_expansivity",
        },
        "parameter_error_confidence": None,
        "experimental_pressure_range_gpa": [0.0, 33.01],
        "experimental_temperature_range_k": [300.0, 1673.0],
        "determination_method": "experimental",
        "record_kind": "published",
        "notes": "Mechanism test fixture.",
        "reference": copy.deepcopy(first["reference"]),
        "pressure_calibration": copy.deepcopy(
            sidecar["records"][first["identifier"]]["pressure_calibration"]
        ),
        "evidence": {
            "scientific_validation": {
                "status": "primary_source_validated",
                "audit_date": "2026-11-20",
                "primary_source_check": {"finding": "Test fixture: Table 1."},
                # The validator refits the bundled observations.
                "primary_data_check": {
                    "status": "bundled",
                    "dataset_identifiers": [
                        "siderite_fe095mn005_litasov_2013_table2_pvt"
                    ],
                    "finding": "Test fixture: Table 2 observations.",
                },
            },
            "parameter_provenance": {"parameters": "Test fixture: Table 1."},
        },
    }
    return {**_common(revision, SIDERITE, "add_record"), "record": record}


def set_request(revision, material, changes) -> dict:
    return {**_common(revision, material, "set_fields"), "changes": changes}


def change(record, path, expected, value) -> dict:
    result = {"record": record, "path": path, "value": value}
    if expected is catalog_revisions.ABSENT:
        result["expected_absent"] = True
    else:
        result["expected"] = expected
    return result


def lacking(root: Path, material: str, field: str) -> list[str]:
    """Records of a material without a top-level field."""
    return [r["identifier"] for r in _records(root, material) if field not in r]


def confidence_request(root: Path, revision: str = "t-confidence") -> dict:
    """Set the 1-sigma confidence on a record that lacks it.

    Prefers qandilite and larnite, whose refit-ledger rows are computed by the
    validator's generic path, then any material.
    """
    materials = sorted(
        p.stem for p in (root / "peritheos/data/materials").glob("*.eosmat")
    )
    for material in ["qandilite", "larnite", *materials]:
        records = lacking(root, material, "parameter_error_confidence")
        if records:
            return set_request(
                revision,
                material,
                [
                    change(
                        records[0],
                        ["parameter_error_confidence"],
                        catalog_revisions.ABSENT,
                        ONE_SIGMA,
                    )
                ],
            )
    raise AssertionError("every record has parameter_error_confidence")


def fixed_request(root: Path, revision: str = "t-fixed") -> dict:
    """Mark a larnite parameter fixed (mechanism only)."""
    for record in _records(root, "larnite"):
        for name, item in record["parameters"].items():
            if "fixed" not in item:
                return set_request(
                    revision,
                    "larnite",
                    [
                        change(
                            record["identifier"],
                            ["parameters", name, "fixed"],
                            catalog_revisions.ABSENT,
                            True,
                        )
                    ],
                )
    raise AssertionError("larnite has no parameter without a fixed flag")


def apply(root, request, derived=False):
    plan = tool.plan_revision(request, root, derived=derived)
    tool.apply_plan(plan, root)
    return plan


def undo(root, identifier, derived=False):
    plan = tool.plan_undo(identifier, root, derived=derived)
    tool.apply_plan(plan, root)
    return plan


def _read_at_revision(relative):
    return _git_file(DEFAULT_REVISION, relative, ROOT)


def assert_history(root, materials):
    """Every historical check a content change must keep."""
    assert verify_ledger(root) == []
    assert check_bindings(root) == []
    assert verify_catalog(root) == len(
        json.loads((root / RECEIPT_PATH).read_bytes())["materials"]
    )
    receipt = json.loads((root / RECEIPT_PATH).read_bytes())
    for material in materials:
        reference = pre_split_bytes(material, root=root)
        binding = receipt["materials"][material]
        assert _sha(reference) == binding["catalog_revision_sha256"]
        original = prior_revision_bytes(reference, material, root)
        assert _sha(original) == binding["before_sha256"]
        assert prior_material_bytes(material, root=root) == prior_material_bytes(
            material, root=ROOT
        )
    for name, content in build_review(root).items():
        assert (root / name).read_bytes() == content
    snapshot = scientific_snapshot(root)
    assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])
    assert verify_check_text(root)[1] == []
    assert verify_record_text(root=root, read=_read_at_revision)[1] == []


# Ledger and zero revisions -------------------------------------------------


def test_the_bundled_ledger_verifies():
    assert verify_ledger() == []
    assert check_bindings() == []
    assert tool.main(["--verify"]) == 0


def test_zero_revisions_change_no_check(repo):
    (repo / LEDGER).write_text(
        catalog_revisions.ledger_text(catalog_revisions.empty_ledger())
    )
    receipt = json.loads((repo / RECEIPT_PATH).read_bytes())
    for binding in receipt["materials"].values():
        binding.pop("catalog_revision_sha256", None)
        binding.pop(catalog_revisions.HEAD_KEY, None)
    receipt.pop(catalog_revisions.SPLIT_RECEIPT_STATEMENT[0], None)
    (repo / RECEIPT_PATH).write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n"
    )
    policy = json.loads((repo / tool.POLICY).read_bytes())
    policy.pop("catalog_revisions", None)
    (repo / tool.POLICY).write_text(
        json.dumps(policy, indent=2, ensure_ascii=False) + "\n"
    )
    data = b"any reference"
    assert prior_revision_bytes(data, SIDERITE, repo) is data
    assert prior_revision_snapshot({"documents": []}, repo) is None
    assert check_bindings(repo) == []
    assert verify_ledger(repo) == []


# Adding records and setting fields -----------------------------------------


def test_add_record_derives_its_content_and_undoes_exactly(repo):
    before = _tree(repo)
    request = add_request(repo)
    identifier = request["record"]["identifier"]
    plan = apply(repo, request, derived=True)
    assert_history(repo, [SIDERITE])
    # Derived, not written by hand: the assessment, check, refit-ledger row,
    # source-audit entry, manifest counts and generated documentation.
    sidecar = json.loads(
        (repo / f"peritheos/data/audits/{SIDERITE}.audit.json").read_bytes()
    )
    assessment = sidecar["records"][identifier]["evidence"]["reproduction"][
        "assessment"
    ]
    rows = {
        r["record_identifier"]: r
        for r in json.loads((repo / tool.REFITS).read_bytes())["records"]
    }
    from scripts.standardize_reconstruction_metadata import assessment_for_record

    merged = merge_audit(
        json.loads((repo / f"peritheos/data/materials/{SIDERITE}.eosmat").read_bytes()),
        sidecar,
    )
    record = next(r for r in merged["eos_records"] if r["identifier"] == identifier)
    assert record["evidence"]["reproduction"]["assessment"] == assessment_for_record(
        record, rows[identifier]
    )
    assert (
        assessment["status"]
        == record["evidence"]["reproduction"]["assessment"]["status"]
    )
    manifest = json.loads((repo / tool.MANIFEST).read_bytes())
    assert manifest["eos_records"] == len(_identifiers(repo))
    audit = json.loads((repo / tool.SOURCE_AUDIT).read_bytes())
    assert identifier in {e["record"] for e in audit["records"]}
    assert f"docs/reconstruction-reports/{SIDERITE}.md" in plan.files
    assert f"docs/reconstruction-reports/figures/{identifier}.svg" in plan.files
    assert f"docs/material-summaries/{SIDERITE}.md" in plan.files
    revision = read_ledger(repo)["revisions"][-1]
    assert set(revision["written_files"]) >= {
        tool.MANIFEST,
        tool.SOURCE_AUDIT,
        tool.REFITS,
    }
    assert list(revision["refit_rows"]) == [identifier]
    # The documented --add route is closed for an existing material.
    from peritheos._eosmat_format import format_eosmat

    path = repo / f"peritheos/data/materials/{SIDERITE}.eosmat"
    split = path.read_bytes()
    path.write_text(format_eosmat(merged), encoding="utf-8")
    with pytest.raises(ValueError, match="only for new materials"):
        split_added_materials([SIDERITE], root=repo, check_snapshot=False)
    path.write_bytes(split)
    undo(repo, request["id"], derived=True)
    assert _tree(repo) == before


def test_set_a_snapshot_field_and_a_fixed_flag_on_two_materials(repo):
    before = _tree(repo)
    confidence = confidence_request(repo)
    fixed = fixed_request(repo)
    apply(repo, confidence)
    apply(repo, fixed)
    assert_history(repo, sorted({confidence["material"], "larnite"}))
    target = confidence["changes"][0]["record"]
    material = confidence["material"]
    reference = json.loads(pre_split_bytes(material, root=repo))
    record = next(r for r in reference["eos_records"] if r["identifier"] == target)
    assert record["parameter_error_confidence"] == ONE_SIGMA
    original = json.loads(pre_split_bytes(material, root=ROOT))
    assert "parameter_error_confidence" not in next(
        r for r in original["eos_records"] if r["identifier"] == target
    )
    apply(
        repo,
        set_request(
            "t-later",
            material,
            [change(target, ["parameter_error_confidence"], ONE_SIGMA, 0.95)],
        ),
    )
    with pytest.raises(tool.Refused, match="out_of_order"):
        tool.plan_undo(confidence["id"], repo)
    undo(repo, "t-later")
    undo(repo, fixed["id"])
    undo(repo, confidence["id"])
    assert _tree(repo) == before


def test_one_request_changes_several_records(pristine):
    records = _records(pristine, "qandilite")
    changes = [
        change(
            record["identifier"],
            ["parameters", "V0", "error"],
            record["parameters"]["V0"]["error"],
            record["parameters"]["V0"]["error"] + 0.1,
        )
        for record in records
    ]
    plan = tool.plan_revision(
        set_request("t-many", "qandilite", changes), pristine, derived=False
    )
    for record in records:
        assert f"Record {record['identifier']} (changed):" in plan.shown


def test_a_headline_change_needs_an_exact_acknowledgement(pristine):
    """Adding the missing errors moves a published fit inside the combined
    error criterion; the headline change must be acknowledged exactly."""
    material, identifier = "coo", "coo_clendenen_1966_murnaghan_1"
    record = next(
        r for r in _records(pristine, material) if r["identifier"] == identifier
    )
    changes = [
        change(
            identifier,
            ["parameters", name, "error"],
            record["parameters"][name].get("error"),
            value,
        )
        for name, value in (("K0", 30.0), ("K0_prime", 3.0))
    ]
    request = set_request("t-ack", material, changes)
    with pytest.raises(tool.Refused, match="check_change") as refusal:
        tool.plan_revision(request, pristine, derived=False)
    expected = {
        identifier: {
            "reproduction": {"from": record["check"]["reproduction"], "to": "close"}
        }
    }
    assert json.dumps(expected) in str(refusal.value)
    request[tool.ACKNOWLEDGE] = {
        identifier: {"reproduction": {"from": "x", "to": "close"}}
    }
    with pytest.raises(tool.Refused, match="check_change"):
        tool.plan_revision(request, pristine, derived=False)
    request[tool.ACKNOWLEDGE] = expected
    plan = tool.plan_revision(request, pristine, derived=False)
    assert plan.revision[tool.ACKNOWLEDGE] == expected
    assert any(
        line.strip().startswith("CHECK CHANGE reproduction") for line in plan.shown
    )
    # The acknowledgement is part of the digest the binding pins.
    altered = {**plan.revision, tool.ACKNOWLEDGE: {}}
    assert catalog_revisions.content_sha256(altered) != plan.revision["content_sha256"]


def test_a_row_that_does_not_take_the_stored_values_is_refused(pristine):
    """The validator takes this record's published values from a fixed table,
    so its row and ``why`` would keep the old error; refuse instead."""
    material, identifier = (
        "akimotoite",
        "akimotoite_dorogokupets_2015_298k_rydberg_stacey",
    )
    record = next(
        r for r in _records(pristine, material) if r["identifier"] == identifier
    )
    current = record["parameters"]["K0"].get("error", catalog_revisions.ABSENT)
    request = set_request(
        "t-row",
        material,
        [change(identifier, ["parameters", "K0", "error"], current, 7.0)],
    )
    with pytest.raises(tool.Refused, match="row_disagrees") as refusal:
        tool.plan_revision(request, pristine, derived=False)
    assert "pinned in scripts/" in str(refusal.value)


# Text revisions and undo ----------------------------------------------------


def _text_change(root, material, identifier, field, value):
    record = next(r for r in _records(root, material) if r["identifier"] == identifier)
    current = record.get(field) if field == "notes" else None
    batches = {
        material: {identifier: {"before": {field: current}, "after": {field: value}}}
    }
    plan = plan_record_text(batches, {}, root=root)
    assert not plan.refusals, plan.refusals
    apply_record_text(plan, "2026-11-20", root)


def _check_text_change(root, material, identifier, note):
    plan = plan_check_text({material: {identifier: {"source_note": note}}}, root=root)
    assert not plan.refusals, plan.refusals
    apply_check_text(plan, root)


def test_text_then_content_then_text_then_undo_all(repo):
    original = _tree(repo)
    sidecar = json.loads(
        (repo / f"peritheos/data/audits/{SIDERITE}.audit.json").read_bytes()
    )
    # A record whose notes already have a stored original, so returning to
    # the earlier text changes neither the stored original nor its date.
    first = next(
        r["identifier"]
        for r in _records(repo, SIDERITE)
        if "notes"
        in sidecar["records"].get(r["identifier"], {}).get("superseded_record_text", {})
        and r.get("notes")
    )
    earlier = next(r for r in _records(repo, SIDERITE) if r["identifier"] == first)[
        "notes"
    ]
    _text_change(repo, SIDERITE, first, "notes", "Rewritten before a content revision.")
    request = add_request(repo, "t-c1")
    added = request["record"]["identifier"]
    apply(repo, request)
    _check_text_change(repo, SIDERITE, added, "Checked against Table 1 of the paper.")
    _text_change(repo, SIDERITE, added, "notes", "Rewritten note of the added record.")
    _text_change(
        repo, SIDERITE, first, "notes", "Rewritten after the content revision."
    )
    assert_history(repo, [SIDERITE])
    revision = read_ledger(repo)["revisions"][-1]
    # The pre-split reference at each stage.
    binding = json.loads((repo / RECEIPT_PATH).read_bytes())["materials"][SIDERITE]
    from peritheos._eosmat_format import format_eosmat

    reference = pre_split_bytes(SIDERITE, root=repo)
    assert _sha(reference) == revision["reference_sha256"]["after"]
    document = json.loads(reference)
    apply_operations(document["eos_records"], revision["reverse"], "reference")
    reference = format_eosmat(document, layout=LAYOUTS[binding.get("layout")]).encode()
    assert _sha(reference) == revision["reference_sha256"]["before"]
    assert reference == pre_split_bytes(SIDERITE, root=ROOT)
    # Undo goes through the later text revisions; they stay on other records.
    plan = undo(repo, "t-c1")
    assert any("later text revisions" in line for line in plan.shown)
    assert added not in _identifiers(repo)
    assert_history(repo, [])
    _text_change(repo, SIDERITE, first, "notes", earlier)
    assert _tree(repo) == original


def test_undo_keeps_a_later_notes_rewrite_of_a_record_whose_row_changed(repo):
    """A changed comparison rewrites the refit row; a later notes rewrite only
    touches the row's text copy, so undo goes through and keeps it."""
    material, identifier = "qandilite", "qandilite_lv_2016_bm3"
    record = next(r for r in _records(repo, material) if r["identifier"] == identifier)
    error = record["parameters"]["K0"]["error"]
    request = set_request(
        "t-row",
        material,
        [change(identifier, ["parameters", "K0", "error"], error, None)],
    )
    try:
        tool.plan_revision(request, repo, derived=False)
    except tool.Refused as refusal:
        assert refusal.reason == "check_change"
        request[tool.ACKNOWLEDGE] = json.loads(
            refusal.detail.split("changes ", 1)[1].split("; the request", 1)[0]
        )
    apply(repo, request)
    assert identifier in read_ledger(repo)["revisions"][-1]["refit_rows"]
    _text_change(repo, material, identifier, "notes", "Rewritten after the revision.")
    plan = undo(repo, "t-row")
    assert any("later text revisions" in line for line in plan.shown)
    record = next(r for r in _records(repo, material) if r["identifier"] == identifier)
    assert record["notes"] == "Rewritten after the revision."
    assert record["parameters"]["K0"]["error"] == error
    rows = {
        r["record_identifier"]: r
        for r in json.loads((repo / tool.REFITS).read_bytes())["records"]
    }
    original = {
        r["record_identifier"]: r
        for r in json.loads((ROOT / tool.REFITS).read_bytes())["records"]
    }
    assert rows[identifier]["source_notes"] == "Rewritten after the revision."
    assert tool._without_text(rows[identifier]) == tool._without_text(
        original[identifier]
    )
    assert_history(repo, [])


# Tampering ------------------------------------------------------------------


@pytest.fixture(scope="module")
def revised(tmp_path_factory):
    """A copy with two revisions of siderite: an added record, then a field of it."""
    root = _copy(tmp_path_factory.mktemp("revised") / "repo")
    saved = catalog_revisions.ROOT
    catalog_revisions.ROOT = root
    try:
        request = add_request(root, "t-c1")
        apply(root, request)
        added = request["record"]["identifier"]
        apply(
            root,
            set_request(
                "t-c2",
                SIDERITE,
                [change(added, ["parameter_error_confidence"], None, ONE_SIGMA)],
            ),
        )
    finally:
        catalog_revisions.ROOT = saved
    return root


@pytest.mark.parametrize(
    "edit",
    [
        "reverse_value",
        "consistent_value",
        "hash",
        "metadata_of_latest",
        "reorder",
        "delete_first",
        "delete_last",
        # Every stored hash of the latest revision is covered by its digest.
        ("files", 0, "before"),
        ("files", 1, "after"),
        ("document_sha256", "before"),
        ("scientific", "catalog_before"),
        ("scientific", "catalog_after"),
        ("policy", "material_sha256", "before"),
        ("policy", "catalog_scientific_sha256", "after"),
        ("written_files",),
        ("split_binding", "before", "after_record_sha256"),
        ("refit_rows",),
    ],
)
def test_a_tampered_ledger_is_detected_everywhere(revised, tmp_path, monkeypatch, edit):
    root = tmp_path / "repo"
    shutil.copytree(revised, root)
    monkeypatch.setattr(catalog_revisions, "ROOT", root)
    snapshot = scientific_snapshot(root)
    ledger = json.loads((root / LEDGER).read_bytes())
    revisions = ledger["revisions"]
    # The fixture's two revisions are the last two, whatever came before.
    first, second = len(revisions) - 2, len(revisions) - 1
    if edit == "reverse_value":
        revisions[second]["reverse"][0]["after"] = 0.95
    elif edit == "consistent_value":
        revisions[second]["forward"][0]["before"] = 0.95
        revisions[second]["reverse"][0]["after"] = 0.95
        revisions[second]["content_sha256"] = catalog_revisions.content_sha256(
            revisions[second]
        )
    elif edit == "hash":
        revisions[first]["reference_sha256"]["after"] = "0" * 64
    elif edit == "metadata_of_latest":
        revisions[second]["decided_by"] = "Someone else"
    elif edit == "reorder":
        revisions[first], revisions[second] = revisions[second], revisions[first]
    elif edit == "delete_first":
        del revisions[first]
    elif edit == "delete_last":
        del revisions[second]
    else:
        node = revisions[second]
        for key in edit[:-1]:
            node = node[list(node)[key]] if isinstance(key, int) else node[key]
        last = edit[-1]
        if last in ("written_files", "refit_rows"):
            node[last] = (
                []
                if last == "written_files"
                else {"tampered": {"before": None, "after": None}}
            )
        else:
            node[last] = "0" * 64
    (root / LEDGER).write_text(catalog_revisions.ledger_text(ledger))
    assert verify_ledger(root)
    assert tool.main(["--verify"], root=root) == 1
    with pytest.raises(ValueError):
        verify_catalog(root)
    with pytest.raises(ValueError):
        prior_material_bytes(SIDERITE, root=root)
    from scripts.review_eos_structure_rights import build_review as review

    if edit not in ("hash", "delete_first"):  # the latest revision changed
        with pytest.raises(ValueError):
            review(root)
    if isinstance(edit, str) and edit not in ("hash", "metadata_of_latest"):
        with pytest.raises((ValueError, AssertionError)):
            assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])


def test_a_binding_without_its_ledger_fails_every_verify(
    revised, tmp_path, monkeypatch
):
    root = tmp_path / "repo"
    shutil.copytree(revised, root)
    monkeypatch.setattr(catalog_revisions, "ROOT", root)
    for content in (
        catalog_revisions.ledger_text(catalog_revisions.empty_ledger()),
        None,
    ):
        if content is None:
            (root / LEDGER).unlink()
        else:
            (root / LEDGER).write_text(content)
        message = "names a catalog revision"
        with pytest.raises(ValueError, match=message):
            verify_catalog(root)
        assert any(message in f for f in verify_check_text(root)[1])
        assert any(
            message in f
            for f in verify_record_text(root=root, read=_read_at_revision)[1]
        )
        with pytest.raises(ValueError, match=message):
            build_review(root)


def test_a_change_no_revision_accounts_for_still_fails(revised, tmp_path, monkeypatch):
    root = tmp_path / "repo"
    shutil.copytree(revised, root)
    monkeypatch.setattr(catalog_revisions, "ROOT", root)
    snapshot = scientific_snapshot(root)
    for document in snapshot["documents"]:
        if document["identifier"] == SIDERITE:
            document["eos_records"][0]["parameters"]["K0"]["value"] += 1
    with pytest.raises(ValueError, match="snapshot changed"):
        assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])
    snapshot = scientific_snapshot(root)
    for document in snapshot["documents"]:
        if document["identifier"] == "mgo":
            document["eos_records"][0]["parameters"]["K0"]["value"] += 1
    with pytest.raises(AssertionError):
        assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])
    path = root / f"peritheos/data/materials/{SIDERITE}.eosmat"
    path.write_bytes(path.read_bytes().replace(b'"value": 120.0', b'"value": 121.0', 1))
    with pytest.raises(ValueError, match="record file changed"):
        prior_material_bytes(SIDERITE, root=root)


def test_a_hand_edited_added_record_text_is_caught(revised, tmp_path, monkeypatch):
    root = tmp_path / "repo"
    shutil.copytree(revised, root)
    monkeypatch.setattr(catalog_revisions, "ROOT", root)
    added = _records(root, SIDERITE)[-1]["identifier"]
    _text_change(root, SIDERITE, added, "notes", "Rewritten note of the added record.")
    assert verify_record_text(root=root, read=_read_at_revision)[1] == []
    path = root / f"peritheos/data/materials/{SIDERITE}.eosmat"
    path.write_bytes(
        path.read_bytes().replace(
            b"Rewritten note of the added record.", b"Silently edited note."
        )
    )
    receipt = json.loads((root / RECEIPT_PATH).read_bytes())
    receipt["materials"][SIDERITE]["after_record_sha256"] = _sha(path.read_bytes())
    (root / RECEIPT_PATH).write_text(
        json.dumps(receipt, indent=2, ensure_ascii=False) + "\n"
    )
    failures = verify_record_text(root=root, read=_read_at_revision)[1]
    assert any("differ in docs/data/primary-eos-refits.json" in f for f in failures)


# Refusals ---------------------------------------------------------------------


@pytest.mark.parametrize(
    ("path", "value", "reason"),
    [
        (["model"], "BM3", "protected_field"),
        (["validity"], {"pressure_range_gpa": [0, 1]}, "protected_field"),
        (["record_kind"], "refit", "protected_field"),
        (["frobnicate"], 1, "protected_field"),
        (["aliases"], ["mgo_b1_tange_2009_vinet"], "protected_field"),
        (["catalog_access"], "explicit_selection", "protected_field"),
        (["fit_datasets"], [], "protected_field"),
        (["reference"], "x", "protected_field"),
        (["determination_method"], "theoretical", "protected_field"),
        (["label"], "x", "protected_field"),
        (["notes"], "x", "protected_field"),
        (["validity", "notes"], ["x"], "protected_field"),
        (["thermal", "model"], "MGD", "protected_field"),
        (["parameters", "Q", "error"], 1.0, "protected_field"),
        # A published coefficient is never replaced by a revision.
        (["parameters", "K0", "value"], 176.797, "published_value"),
        (["parameters", "K0_prime", "value"], 5.148, "published_value"),
        (["thermal", "parameters", "Tr", "value"], 298.0, "published_value"),
        (["parameters", "K0", "error"], -1, "value"),
        (["parameter_error_confidence"], 1.5, "value"),
        (["experimental_pressure_range_gpa"], [5, 1], "range"),
    ],
)
def test_set_fields_refuses_paths_outside_the_allow_list(pristine, path, value, reason):
    record = _records(pristine, "qandilite")[0]
    try:
        expected = catalog_revisions.read_path(record, path)
    except ValueError:
        expected = catalog_revisions.ABSENT
    request = set_request(
        "t", "qandilite", [change(record["identifier"], path, expected, value)]
    )
    with pytest.raises(tool.Refused, match=reason) as refusal:
        tool.plan_revision(request, pristine, derived=False)
    if reason == "protected_field" and path[0] not in ("label", "notes", "validity"):
        assert "settable:" in str(refusal.value)


def _covariance_record(root: Path) -> tuple[str, dict]:
    """A record whose covariance names free parameters (picked at run time)."""
    for path in sorted((root / "peritheos/data/materials").glob("*.eosmat")):
        for record in json.loads(path.read_bytes())["eos_records"]:
            covariance = record.get("parameter_covariance")
            if (
                isinstance(covariance, dict)
                and len(covariance.get("parameter_order", [])) == 2
                and all(
                    n in record["parameters"] for n in covariance["parameter_order"]
                )
                and not record.get("thermal")
            ):
                return path.stem, record
    raise AssertionError("no record with a two-parameter covariance")


@pytest.mark.parametrize(
    ("edit", "reason"),
    [
        (lambda c: c.update(parameter_order=["K0", "Q"]), "not a parameter"),
        (lambda c: c.update(matrix=[[1.0, 0.0]]), "must be 2x2"),
        (lambda c: c.update(matrix=[[1.0, 0.5], [0.4, 1.0]]), "not symmetric"),
        (lambda c: c.update(matrix=[[1.0, 2.0], [2.0, 1.0]]), "positive semidefinite"),
        (lambda c: c.update(matrix=[[1.0, "x"], ["x", 1.0]]), "finite numbers"),
    ],
)
def test_a_covariance_is_checked_before_staging(pristine, edit, reason):
    material, record = _covariance_record(pristine)
    covariance = copy.deepcopy(record["parameter_covariance"])
    edit(covariance)
    request = set_request(
        "t",
        material,
        [
            change(
                record["identifier"],
                ["parameter_covariance"],
                record["parameter_covariance"],
                covariance,
            )
        ],
    )
    with pytest.raises(tool.Refused, match=reason):
        tool.plan_revision(request, pristine, derived=False)


def test_fixing_a_covariance_parameter_is_refused(pristine):
    material, record = _covariance_record(pristine)
    name = record["parameter_covariance"]["parameter_order"][0]
    current = record["parameters"][name].get("fixed", catalog_revisions.ABSENT)
    request = set_request(
        "t",
        material,
        [change(record["identifier"], ["parameters", name, "fixed"], current, True)],
    )
    with pytest.raises(tool.Refused, match="is fixed, not a free parameter"):
        tool.plan_revision(request, pristine, derived=False)


@pytest.mark.parametrize(
    ("path", "value", "reason"),
    [
        (["experimental_pressure_range_gpa"], [0, 3000], "outside"),
        (["experimental_pressure_range_gpa"], [0, float("inf")], "finite numbers"),
        (["experimental_temperature_range_k"], [0, 300], "outside"),
        (["pressure_range_status"], "approximately", "is not one of"),
        (["parameters", "K0", "error"], float("nan"), "value"),
        (["temperature_ref"], 310.0, "outside the stated temperature range"),
    ],
)
def test_ranges_and_values_are_validated(pristine, path, value, reason):
    record = _records(pristine, "qandilite")[0]
    try:
        expected = catalog_revisions.read_path(record, path)
    except ValueError:
        expected = catalog_revisions.ABSENT
    request = set_request(
        "t", "qandilite", [change(record["identifier"], path, expected, value)]
    )
    with pytest.raises(tool.Refused, match=reason):
        tool.plan_revision(request, pristine, derived=False)


def test_an_added_record_is_validated_like_a_change(pristine):
    for edit, reason in (
        (lambda r: r.update(temperature_ref=298.0), "differs from thermal Tr"),
        (lambda r: r.update(experimental_pressure_range_gpa=[0.0, 5000.0]), "outside"),
        (lambda r: r.update(pressure_range_status="roughly"), "is not one of"),
        (lambda r: r["parameters"]["K0"].update(error=-2.0), "error must be"),
    ):
        request = add_request(pristine)
        edit(request["record"])
        with pytest.raises(tool.Refused, match=reason):
            tool.plan_revision(request, pristine, derived=False)


def test_removing_or_emptying_validity_is_refused(pristine):
    for material in sorted(
        p.stem for p in (pristine / "peritheos/data/materials").glob("*.eosmat")
    ):
        record = next(
            (
                r
                for r in _records(pristine, material)
                if isinstance(r.get("validity"), dict)
            ),
            None,
        )
        if record is not None:
            break
    validity = record["validity"]
    for change_ in (
        {
            "record": record["identifier"],
            "path": ["validity"],
            "expected": validity,
            "remove": True,
        },
        {
            "record": record["identifier"],
            "path": ["validity"],
            "expected": validity,
            "value": {},
        },
    ):
        with pytest.raises(tool.Refused, match="protected_field"):
            tool.plan_revision(
                set_request("t", material, [change_]), pristine, derived=False
            )


def test_values_compare_numerically_and_no_ops_are_refused(pristine):
    record = next(
        r
        for r in _records(pristine, "qandilite")
        if isinstance(r["parameters"]["K0"].get("error"), float)
        and r["parameters"]["K0"]["error"].is_integer()
    )
    error = record["parameters"]["K0"]["error"]
    path = ["parameters", "K0", "error"]
    # An integer expectation matches the stored float...
    request = set_request(
        "t", "qandilite", [change(record["identifier"], path, int(error), error)]
    )
    # ...and setting the same number again changes nothing.
    with pytest.raises(tool.Refused, match="no_change"):
        tool.plan_revision(request, pristine, derived=False)
    request["changes"][0]["value"] = int(error)
    with pytest.raises(tool.Refused, match="no_change"):
        tool.plan_revision(request, pristine, derived=False)
    for expected in (error + 1, catalog_revisions.ABSENT):
        request = set_request(
            "t",
            "qandilite",
            [change(record["identifier"], path, expected, error + 2)],
        )
        with pytest.raises(tool.Refused, match="expected_mismatch"):
            tool.plan_revision(request, pristine, derived=False)


@pytest.mark.parametrize(
    ("mutation", "reason"),
    [
        (lambda r: r.pop("pressure_calibration"), "record"),
        (lambda r: r["evidence"].pop("parameter_provenance"), "record"),
        (
            lambda r: r["evidence"]["scientific_validation"].pop(
                "primary_source_check"
            ),
            "record",
        ),
        (
            lambda r: r["evidence"].update(
                reproduction={"assessment": {"status": "parity"}}
            ),
            "derived_field",
        ),
        (lambda r: r.update(frobnicate=1), "record"),
        (lambda r: r.update(check={}), "record"),
        (lambda r: r.update(default=True), "undeclared_change"),
        (lambda r: r.update(identifier="qandilite_lv_2016_bm2"), "identifier_taken"),
        (lambda r: r.update(aliases=["mgo_b1_tange_2009_vinet"]), "identifier_taken"),
        (lambda r: r.update(model="NoSuchModel"), "invalid|not_loadable"),
    ],
)
def test_add_record_refuses_what_the_catalog_rejects(pristine, mutation, reason):
    request = add_request(pristine)
    mutation(request["record"])
    with pytest.raises(tool.Refused, match=reason):
        tool.plan_revision(request, pristine, derived=False)


@pytest.mark.parametrize(
    "edit",
    [
        lambda r: r["parameters"]["V0"].update(value=-5.0),
        lambda r: r["parameters"]["K0"].update(value=1e-10),
        lambda r: r["parameters"]["K0_prime"].update(value=-50.0),
    ],
)
def test_an_unusable_record_is_refused_without_a_traceback(
    pristine, tmp_path, capsys, edit
):
    request = add_request(pristine)
    edit(request["record"])
    with pytest.raises(tool.Refused):
        tool.plan_revision(request, pristine, derived=False)
    path = tmp_path / "request.json"
    path.write_text(json.dumps(request))
    assert tool.main([str(path), "--check"], root=pristine) == 1
    assert capsys.readouterr().out.startswith("REFUSED ")


def test_add_on_a_bundled_material_is_a_clean_refusal():
    import subprocess
    import sys

    result = subprocess.run(
        [sys.executable, "scripts/split_eosmat_audits.py", "--add", "qandilite"],
        cwd=ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 1
    assert result.stdout.startswith("REFUSED ")
    assert "Traceback" not in result.stderr


def test_requests_need_a_complete_rights_statement(pristine):
    request = add_request(pristine)
    del request["rights"]
    with pytest.raises(tool.Refused, match=r"^rights: missing \['rights'\]"):
        tool.check_request(request)
    request = add_request(pristine)
    del request["rights"]["same_basis"]
    with pytest.raises(tool.Refused, match=r"^rights: add_record needs"):
        tool.check_request(request)
    request = add_request(pristine)
    request["rights"]["same_basis"] = False
    with pytest.raises(tool.Refused, match="differences"):
        tool.check_request(request)
    request["rights"]["differences"] = "A different table of the same article."
    tool.check_request(request)
    request = confidence_request(pristine)
    del request["rights"]["statement"]
    with pytest.raises(tool.Refused, match=r"^rights: set_fields needs"):
        tool.check_request(request)
    request = add_request(pristine)
    request["default_records"] = {}
    with pytest.raises(tool.Refused, match="undeclared_change"):
        tool.check_request(request)


# Command line ------------------------------------------------------------------


def test_command_line_check_apply_verify_and_undo(repo, tmp_path, capsys):
    before = _tree(repo)
    path = tmp_path / "request.json"
    request = confidence_request(repo)
    path.write_text(json.dumps(request))
    assert tool.main([str(path), "--check"], root=repo) == 0
    output = capsys.readouterr().out
    assert "parameter_error_confidence: absent -> 0.6826894921370859" in output
    assert "check block (record file):" in output
    assert output.count('"reproduction"') == 1
    assert _tree(repo) == before
    assert tool.main([str(path)], root=repo) == 0
    assert tool.main(["--verify"], root=repo) == 0
    assert tool.main(["--undo", request["id"], "--check"], root=repo) == 0
    assert tool.main(["--undo", request["id"]], root=repo) == 0
    assert "returns to absent" in capsys.readouterr().out
    assert _tree(repo) == before


def test_bad_input_is_refused_without_a_traceback(repo, tmp_path, capsys):
    assert tool.main([str(tmp_path / "missing.json")], root=repo) == 1
    assert "REFUSED request: cannot read" in capsys.readouterr().out
    (tmp_path / "bad.json").write_text("{")
    assert tool.main([str(tmp_path / "bad.json")], root=repo) == 1
    assert "not valid JSON" in capsys.readouterr().out
    request = confidence_request(repo)
    (tmp_path / "ok.json").write_text(json.dumps(request))
    ledger = (repo / LEDGER).read_bytes()
    (repo / LEDGER).write_text("[]")
    assert tool.main([str(tmp_path / "ok.json")], root=repo) == 1
    assert "REFUSED" in capsys.readouterr().out
    (repo / LEDGER).write_bytes(ledger)
    before = _tree(repo)
    material = repo / f"peritheos/data/materials/{request['material']}.eosmat"
    os.chmod(material, 0o444)
    try:
        assert tool.main([str(tmp_path / "ok.json")], root=repo) == 1
        assert "REFUSED not_writable" in capsys.readouterr().out
    finally:
        os.chmod(material, 0o644)
    assert _tree(repo) == before
