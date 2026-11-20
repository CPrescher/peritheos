"""Catalog revisions: reversible, receipt-backed changes of bundled content.

Every test that changes data works on a temporary copy of the catalog files,
the receipts and the rights review. The records written here are mechanism
fixtures, not curated scientific claims.
"""

import copy
import hashlib
import json
import shutil
from pathlib import Path

import pytest

from peritheos._eosmat_audit import merge_audit
from peritheos.materials import Material
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
    prior_revision_bytes,
    prior_revision_snapshot,
    read_ledger,
    verify_ledger,
)
from scripts.migrate_eosmat_v5 import prior_material_bytes
from scripts.review_eos_structure_rights import (
    GROUPS,
    INDEX,
    POLICY,
    REVIEW,
    build_review,
)
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
ADDED = "siderite_fe095mn005_litasov_2013_bm3_thermal"
ONE_SIGMA = 0.6826894921370859
GARNET = json.loads(
    (
        ROOT / "docs/data/garnet-angel-2022/dataset-integration-2026-11-15.json"
    ).read_bytes()
)


def _sha(data: bytes) -> str:
    return hashlib.sha256(data).hexdigest()


@pytest.fixture
def catalog(tmp_path, monkeypatch):
    """A writable copy of everything a catalog revision reads or writes."""
    for relative in ("peritheos/data/materials", "peritheos/data/audits"):
        shutil.copytree(ROOT / relative, tmp_path / relative)
    (tmp_path / "peritheos/data/datasets").symlink_to(
        ROOT / "peritheos/data/datasets", target_is_directory=True
    )
    for relative in (RECEIPT_PATH, LEDGER, POLICY, REVIEW, INDEX, GROUPS):
        (tmp_path / relative).parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(ROOT / relative, tmp_path / relative)
    # assert_snapshot_lineage reads the revision ledger of this copy.
    monkeypatch.setattr(catalog_revisions, "ROOT", tmp_path)
    return tmp_path


def _tree(root):
    """Bytes of every file a revision may write."""
    paths = [
        *sorted((root / "peritheos/data/materials").glob("*.eosmat")),
        *sorted((root / "peritheos/data/audits").glob("*.audit.json")),
        *(root / p for p in (RECEIPT_PATH, LEDGER, POLICY, REVIEW, INDEX, GROUPS)),
    ]
    return {p.relative_to(root).as_posix(): p.read_bytes() for p in paths}


def _common(identifier, material, kind):
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


def add_request(root, identifier="t-add"):
    """Add the Litasov (2013) Table 1 thermal BM3 record to siderite."""
    record_file = json.loads(
        (root / f"peritheos/data/materials/{SIDERITE}.eosmat").read_bytes()
    )
    sidecar = json.loads(
        (root / f"peritheos/data/audits/{SIDERITE}.audit.json").read_bytes()
    )
    first = record_file["eos_records"][0]
    record = {
        "identifier": ADDED,
        "label": "Litasov et al. (2013), natural Fe0.95Mn0.05CO3, thermal BM3",
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
            "reproduction": {
                "assessment": {
                    "status": "not_refittable",
                    "scope": "parameter_agreement",
                    "source_status": "not_refittable",
                    "basis": "Test fixture: not compared.",
                }
            },
        },
    }
    return {**_common(identifier, SIDERITE, "add_record"), "record": record}


def set_request(identifier, material, record, path, expected, value):
    change = {"record": record, "path": path, "value": value}
    if expected is catalog_revisions.ABSENT:
        change["expected_absent"] = True
    else:
        change["expected"] = expected
    return {**_common(identifier, material, "set_fields"), "changes": [change]}


def confidence_request(identifier="t-confidence"):
    return set_request(
        identifier,
        "qandilite",
        "qandilite_lv_2016_bm2",
        ["parameter_error_confidence"],
        catalog_revisions.ABSENT,
        ONE_SIGMA,
    )


def fixed_request(identifier="t-fixed"):
    # A fixed flag on another material (mechanism only).
    return set_request(
        identifier,
        "larnite",
        "larnite_sagatova_2021_gga_300k_vinet",
        ["parameters", "K0_prime", "fixed"],
        catalog_revisions.ABSENT,
        True,
    )


def apply(root, request):
    plan = tool.plan_revision(request, root)
    tool.apply_plan(plan, root)
    return plan


def undo(root, identifier):
    plan = tool.plan_undo(identifier, root)
    tool.apply_plan(plan, root)


def _read_at_revision(relative):
    return _git_file(DEFAULT_REVISION, relative, ROOT)


def assert_history(root, materials):
    """Every historical check that a content change must keep."""
    assert verify_ledger(root) == []
    assert verify_catalog(root) == len(
        json.loads((root / RECEIPT_PATH).read_bytes())["materials"]
    )
    receipt = json.loads((root / RECEIPT_PATH).read_bytes())
    for material in materials:
        reference = pre_split_bytes(material, root=root)
        binding = receipt["materials"][material]
        assert _sha(reference) == binding["catalog_revision_sha256"]
        # The original pre-split file and the v5-migration predecessor.
        assert prior_revision_bytes(reference, material, root) == pre_split_bytes(
            material
        )
        assert prior_material_bytes(material, root=root) == prior_material_bytes(
            material
        )
    for name, content in build_review(root).items():
        assert (root / name).read_bytes() == content
    review = json.loads((root / REVIEW).read_bytes())
    assert review["catalog_scientific_sha256"] == _sha(
        json.dumps(
            scientific_snapshot(root), sort_keys=True, separators=(",", ":")
        ).encode()
    )
    snapshot = scientific_snapshot(root)
    assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])
    reversed_ = prior_revision_snapshot(snapshot, root)
    assert catalog_revisions.fingerprint(
        reversed_ or snapshot
    ) == catalog_revisions.fingerprint(scientific_snapshot())
    assert verify_check_text(root)[1] == []
    assert verify_record_text(root=root, read=_read_at_revision)[1] == []


def test_ledger_verifies_and_zero_revisions_change_nothing():
    assert read_ledger()["revisions"] == []
    assert verify_ledger() == []
    assert tool.main(["--verify"]) == 0
    receipt = json.loads((ROOT / RECEIPT_PATH).read_bytes())
    assert not any(
        "catalog_revision_sha256" in b for b in receipt["materials"].values()
    )
    assert "catalog_revision" not in receipt
    assert "catalog_revisions" not in json.loads((ROOT / POLICY).read_bytes())
    data = pre_split_bytes(SIDERITE)
    assert prior_revision_bytes(data, SIDERITE) is data
    assert prior_revision_snapshot(scientific_snapshot()) is None


def test_add_record_keeps_every_historical_check_and_undoes_exactly(catalog):
    before = _tree(catalog)
    plan = apply(catalog, add_request(catalog))
    assert set(plan.files) == {
        f"peritheos/data/materials/{SIDERITE}.eosmat",
        f"peritheos/data/audits/{SIDERITE}.audit.json",
        RECEIPT_PATH,
        LEDGER,
        POLICY,
        REVIEW,
        INDEX,
        GROUPS,
    }
    assert_history(catalog, [SIDERITE])
    revision = read_ledger(catalog)["revisions"][0]
    assert revision["rights"]["same_basis"] is True
    assert list(revision["scientific"]["added_records"]) == [ADDED]
    # The record is executable, its check generated, defaults unchanged.
    record_file = json.loads(
        (catalog / f"peritheos/data/materials/{SIDERITE}.eosmat").read_bytes()
    )
    sidecar = json.loads(
        (catalog / f"peritheos/data/audits/{SIDERITE}.audit.json").read_bytes()
    )
    merged = merge_audit(record_file, sidecar)
    assert record_file["eos_records"][-1]["check"]["reproduction"] == "untested"
    assert (
        merged["default_records"]
        == json.loads(
            (ROOT / f"peritheos/data/materials/{SIDERITE}.eosmat").read_bytes()
        )["default_records"]
    )
    material = Material.from_eosmat(merged, require_primary_validation=False)
    added = next(r for r in material.eos_records if r.identifier == ADDED)
    assert added.pressure(280.0, 1000.0) > 0
    binding = json.loads((catalog / RECEIPT_PATH).read_bytes())["materials"][SIDERITE]
    assert binding["before_sha256"] == _sha(pre_split_bytes(SIDERITE))
    assert len(record_file["eos_records"]) - 1 in binding["added_check"]
    review = json.loads((catalog / REVIEW).read_bytes())
    assert f"eos/{SIDERITE}/{ADDED}" in {e["key"] for e in review["entries"]}
    # The documented --add route would supersede the binding that roots the
    # chain: a complete inline file over the record file is refused.
    from peritheos._eosmat_format import format_eosmat

    path = catalog / f"peritheos/data/materials/{SIDERITE}.eosmat"
    split = path.read_bytes()
    path.write_text(format_eosmat(merged), encoding="utf-8")
    with pytest.raises(ValueError, match="catalog revisions"):
        split_added_materials([SIDERITE], root=catalog, check_snapshot=False)
    with pytest.raises(ValueError, match="catalog revisions"):
        split_added_materials(root=catalog, check_snapshot=False)
    path.write_bytes(split)
    undo(catalog, "t-add")
    assert _tree(catalog) == before


def test_set_snapshot_field_and_fixed_flag_on_two_materials(catalog):
    before = _tree(catalog)
    apply(catalog, confidence_request())
    apply(catalog, fixed_request())
    assert_history(catalog, ["qandilite", "larnite"])
    # The revised reference carries the change; history does not.
    reference = json.loads(pre_split_bytes("qandilite", root=catalog))
    record = next(
        r
        for r in reference["eos_records"]
        if r["identifier"] == "qandilite_lv_2016_bm2"
    )
    assert record["parameter_error_confidence"] == ONE_SIGMA
    original = json.loads(pre_split_bytes("qandilite"))
    assert "parameter_error_confidence" not in original["eos_records"][0]
    larnite = json.loads(pre_split_bytes("larnite", root=catalog))
    assert larnite["eos_records"][0]["parameters"]["K0_prime"]["fixed"] is True
    # The recorded key order of the qandilite record names the new field.
    binding = json.loads((catalog / RECEIPT_PATH).read_bytes())["materials"][
        "qandilite"
    ]
    assert "parameter_error_confidence" in binding["object_orders"][0]["keys"]
    undo(catalog, "t-fixed")
    undo(catalog, "t-confidence")
    assert _tree(catalog) == before


def _record_text_change(root, material, identifier, notes):
    record = next(
        r
        for r in json.loads(
            (root / f"peritheos/data/materials/{material}.eosmat").read_bytes()
        )["eos_records"]
        if r["identifier"] == identifier
    )
    batches = {
        material: {
            identifier: {
                "before": {"notes": record.get("notes")},
                "after": {"notes": notes},
            }
        }
    }
    plan = plan_record_text(batches, {}, root=root)
    assert not plan.refusals, plan.refusals
    apply_record_text(plan, "2026-11-20", root)


def _check_text_change(root, material, identifier, note):
    plan = plan_check_text({material: {identifier: {"source_note": note}}}, root=root)
    assert not plan.refusals, plan.refusals
    apply_check_text(plan, root)


def test_mixed_sequence_with_text_revisions_recovers_every_stage(catalog):
    original = _tree(catalog)
    receipt = json.loads((catalog / RECEIPT_PATH).read_bytes())
    binding = receipt["materials"][SIDERITE]
    # The material already has bound text and check-text revisions.
    assert binding["superseded_record_text"] and binding["check_text"]
    first = f"{SIDERITE}_litasov_2013_bm3"
    apply(catalog, add_request(catalog, "t-c1"))
    _record_text_change(catalog, SIDERITE, first, "Rewritten after a content revision.")
    _check_text_change(
        catalog, SIDERITE, ADDED, "Checked against Table 1 of the paper."
    )
    _record_text_change(catalog, SIDERITE, ADDED, "Rewritten note of the added record.")
    apply(
        catalog,
        set_request(
            "t-c2", SIDERITE, first, ["parameter_error_confidence"], None, ONE_SIGMA
        ),
    )
    _check_text_change(
        catalog, SIDERITE, first, "Checked again after the second revision."
    )
    assert_history(catalog, [SIDERITE])
    revisions = read_ledger(catalog)["revisions"]
    assert [r["id"] for r in revisions] == ["t-c1", "t-c2"]
    # The pre-split reference at each stage.
    receipt = json.loads((catalog / RECEIPT_PATH).read_bytes())
    layout = LAYOUTS[receipt["materials"][SIDERITE].get("layout")]
    from peritheos._eosmat_format import format_eosmat

    reference = pre_split_bytes(SIDERITE, root=catalog)
    assert _sha(reference) == revisions[1]["reference_sha256"]["after"]
    for revision in reversed(revisions):
        document = json.loads(reference)
        apply_operations(document["eos_records"], revision["reverse"], "reference")
        reference = format_eosmat(document, layout=layout).encode()
        assert _sha(reference) == revision["reference_sha256"]["before"]
    assert reference == pre_split_bytes(SIDERITE)
    # Undo order: only the latest, and only while its files are unchanged.
    with pytest.raises(tool.Refused, match="out_of_order"):
        tool.plan_undo("t-c1", catalog)
    with pytest.raises(tool.Refused, match="changed_since"):
        tool.plan_undo("t-c2", catalog)
    assert _tree(catalog) != original


def _tamper(root, change):
    path = root / LEDGER
    ledger = json.loads(path.read_bytes())
    change(ledger["revisions"])
    path.write_text(catalog_revisions.ledger_text(ledger), encoding="utf-8")


@pytest.mark.parametrize(
    "change",
    [
        "reverse_value",
        "consistent_value",
        "hash",
        "reorder",
        "delete_first",
        "delete_last",
    ],
)
def test_tampered_ledger_is_detected_everywhere(catalog, change):
    apply(catalog, add_request(catalog, "t-c1"))
    first = f"{SIDERITE}_litasov_2013_bm3"
    apply(
        catalog,
        set_request(
            "t-c2", SIDERITE, first, ["parameter_error_confidence"], None, ONE_SIGMA
        ),
    )
    snapshot = scientific_snapshot(catalog)

    def edit(revisions):
        if change == "reverse_value":
            revisions[1]["reverse"][0]["after"] = 0.95
        elif change == "consistent_value":
            revisions[1]["forward"][0]["before"] = 0.95
            revisions[1]["reverse"][0]["after"] = 0.95
        elif change == "hash":
            revisions[0]["reference_sha256"]["after"] = "0" * 64
        elif change == "reorder":
            revisions.reverse()
        elif change == "delete_first":
            del revisions[0]
        else:
            del revisions[1]

    _tamper(catalog, edit)
    assert verify_ledger(catalog)
    with pytest.raises(ValueError):
        prior_material_bytes(SIDERITE, root=catalog)
    with pytest.raises(ValueError):
        verify_catalog(catalog)
    if change != "hash":  # reference hashes are not part of the snapshot
        with pytest.raises((ValueError, AssertionError)):
            assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])


def test_a_change_no_revision_accounts_for_still_fails(catalog):
    apply(catalog, confidence_request())
    snapshot = scientific_snapshot(catalog)
    for document in snapshot["documents"]:
        if document["identifier"] == "qandilite":
            document["eos_records"][0]["parameters"]["K0"]["value"] += 1
    with pytest.raises(ValueError, match="snapshot changed"):
        assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])
    snapshot = scientific_snapshot(catalog)
    for document in snapshot["documents"]:
        if document["identifier"] == "mgo":
            document["eos_records"][0]["parameters"]["K0"]["value"] += 1
    with pytest.raises(AssertionError):
        assert_snapshot_lineage(snapshot, GARNET["scientific_sha256_before"])
    # A hand edit of the revised files breaks the binding, too.
    path = catalog / "peritheos/data/materials/qandilite.eosmat"
    path.write_bytes(path.read_bytes().replace(b"172.0", b"173.0", 1))
    with pytest.raises(ValueError, match="record file changed"):
        prior_material_bytes("qandilite", root=catalog)


def test_undeclared_changes_are_refused(catalog):
    before = _tree(catalog)
    request = add_request(catalog)
    request["record"]["default"] = True
    with pytest.raises(tool.Refused, match="undeclared_change"):
        tool.plan_revision(request, catalog)
    request = add_request(catalog)
    request["default_records"] = {"equilibrium": ADDED}
    with pytest.raises(tool.Refused, match="undeclared_change"):
        tool.plan_revision(request, catalog)
    for path in (["label"], ["notes"], ["check"], ["evidence"], ["validity", "notes"]):
        request = set_request(
            "t", "qandilite", "qandilite_lv_2016_bm2", path, None, "x"
        )
        with pytest.raises(tool.Refused, match="protected_field"):
            tool.plan_revision(request, catalog)
    request = add_request(catalog)
    request["record"]["identifier"] = "qandilite_lv_2016_bm2"
    with pytest.raises(tool.Refused, match="identifier_taken"):
        tool.plan_revision(request, catalog)
    request = add_request(catalog)
    del request["record"]["evidence"]["reproduction"]
    with pytest.raises(tool.Refused, match="audit"):
        tool.plan_revision(request, catalog)
    assert _tree(catalog) == before


def test_a_mismatched_expected_value_is_refused(catalog):
    before = _tree(catalog)
    for expected in (0.95, None):
        request = set_request(
            "t",
            "qandilite",
            "qandilite_lv_2016_bm2",
            ["parameter_error_confidence"],
            expected,
            ONE_SIGMA,
        )
        with pytest.raises(tool.Refused, match="expected_mismatch"):
            tool.plan_revision(request, catalog)
    request = set_request(
        "t",
        SIDERITE,
        f"{SIDERITE}_litasov_2013_bm3",
        ["parameter_error_confidence"],
        catalog_revisions.ABSENT,
        ONE_SIGMA,
    )
    with pytest.raises(tool.Refused, match="expected_mismatch"):
        tool.plan_revision(request, catalog)
    # null and 0 are different values; so are 120 and 120.0.
    request = set_request(
        "t",
        SIDERITE,
        f"{SIDERITE}_litasov_2013_bm3",
        ["parameters", "K0", "value"],
        120,
        121.0,
    )
    with pytest.raises(tool.Refused, match="expected_mismatch"):
        tool.plan_revision(request, catalog)
    assert _tree(catalog) == before


def test_an_added_record_needs_a_rights_statement(catalog):
    request = add_request(catalog)
    del request["rights"]
    with pytest.raises(tool.Refused, match="request|rights"):
        tool.check_request(request)
    request = add_request(catalog)
    del request["rights"]["same_basis"]
    with pytest.raises(tool.Refused, match="rights"):
        tool.check_request(request)
    request = add_request(catalog)
    request["rights"]["same_basis"] = False
    with pytest.raises(tool.Refused, match="differences"):
        tool.check_request(request)
    request["rights"]["differences"] = "A different table of the same article."
    tool.check_request(request)
    request = confidence_request()
    del request["rights"]["statement"]
    with pytest.raises(tool.Refused, match="rights"):
        tool.check_request(request)


def test_dry_run_writes_nothing(catalog, capsys):
    before = _tree(catalog)
    path = catalog / "request.json"
    path.write_text(json.dumps(confidence_request()))
    plan = tool.plan_revision(tool.read_request(path), catalog)
    assert "qandilite_lv_2016_bm2/parameter_error_confidence" in "\n".join(plan.shown)
    assert _tree(catalog) == before
