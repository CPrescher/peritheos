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
    """Set the 1-sigma confidence on a qandilite record that lacks it."""
    records = lacking(root, "qandilite", "parameter_error_confidence")
    assert records, "qandilite has no record without parameter_error_confidence"
    return set_request(
        revision,
        "qandilite",
        [
            change(
                records[0],
                ["parameter_error_confidence"],
                catalog_revisions.ABSENT,
                ONE_SIGMA,
            )
        ],
    )


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
    assert f"docs/material-summaries/{SIDERITE}.md" in plan.files
    revision = read_ledger(repo)["revisions"][-1]
    assert set(revision["derived_files"]) >= {
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
    assert_history(repo, ["qandilite", "larnite"])
    target = confidence["changes"][0]["record"]
    reference = json.loads(pre_split_bytes("qandilite", root=repo))
    record = next(r for r in reference["eos_records"] if r["identifier"] == target)
    assert record["parameter_error_confidence"] == ONE_SIGMA
    original = json.loads(pre_split_bytes("qandilite", root=ROOT))
    assert "parameter_error_confidence" not in next(
        r for r in original["eos_records"] if r["identifier"] == target
    )
    apply(
        repo,
        set_request(
            "t-later",
            "qandilite",
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
    else:
        del revisions[second]
    (root / LEDGER).write_text(catalog_revisions.ledger_text(ledger))
    assert verify_ledger(root)
    assert tool.main(["--verify"], root=root) == 1
    with pytest.raises(ValueError):
        verify_catalog(root)
    with pytest.raises(ValueError):
        prior_material_bytes(SIDERITE, root=root)
    if edit not in ("hash", "metadata_of_latest"):  # not part of the snapshot
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
        (["parameters", "Q", "value"], 1.0, "protected_field"),
        (["parameters", "K0", "value"], "x", "value"),
        (["parameters", "K0", "error"], -1, "value"),
        (["parameter_error_confidence"], 1.5, "value"),
        (["experimental_pressure_range_gpa"], [5, 1], "value"),
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
    record = _records(pristine, "qandilite")[0]
    k0 = record["parameters"]["K0"]["value"]
    assert isinstance(k0, float) and k0.is_integer()
    # An integer expectation matches the stored float...
    request = set_request(
        "t",
        "qandilite",
        [change(record["identifier"], ["parameters", "K0", "value"], int(k0), k0)],
    )
    # ...and setting the same number again changes nothing.
    with pytest.raises(tool.Refused, match="no_change"):
        tool.plan_revision(request, pristine, derived=False)
    request["changes"][0]["value"] = int(k0)
    with pytest.raises(tool.Refused, match="no_change"):
        tool.plan_revision(request, pristine, derived=False)
    for expected in (k0 + 1, catalog_revisions.ABSENT):
        request = set_request(
            "t",
            "qandilite",
            [
                change(
                    record["identifier"],
                    ["parameters", "K0", "value"],
                    expected,
                    k0 + 2,
                )
            ],
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
    (tmp_path / "ok.json").write_text(json.dumps(confidence_request(repo)))
    ledger = (repo / LEDGER).read_bytes()
    (repo / LEDGER).write_text("[]")
    assert tool.main([str(tmp_path / "ok.json")], root=repo) == 1
    assert "REFUSED" in capsys.readouterr().out
    (repo / LEDGER).write_bytes(ledger)
    before = _tree(repo)
    material = repo / "peritheos/data/materials/qandilite.eosmat"
    os.chmod(material, 0o444)
    try:
        assert tool.main([str(tmp_path / "ok.json")], root=repo) == 1
        assert "REFUSED not_writable" in capsys.readouterr().out
    finally:
        os.chmod(material, 0o644)
    assert _tree(repo) == before
