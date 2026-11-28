"""Material notes revisions keep the replaced original in the audit sidecar.

Covers the material-level ``superseded_material_text`` audit block and the
``"(material)"`` input of ``scripts/apply_record_text.py``: split and merge,
the audit-split receipt binding and byte-exact pre-split recovery, a second
rewrite, a revert, the guards, ``--verify`` of every maintenance tool, and a
catalog revision applied after and before a notes rewrite. Every change is
made on a temporary copy.
"""

import copy
import json
import shutil
from pathlib import Path

import jsonschema
import pytest

from peritheos import convert_eosmat
from peritheos._eosmat_audit import (
    MATERIAL_AUDIT_KEY,
    merge_audit,
    split_document,
)
from peritheos._eosmat_format import format_audit, format_eosmat
from peritheos.errors import EosmatError
from scripts import apply_catalog_revision as revision_tool
from scripts import apply_record_text as tool
from scripts.apply_check_text import verify as verify_check_text
from scripts.split_eosmat_audits import (
    RECEIPT,
    pre_split_bytes,
    pre_split_view,
    verify_catalog,
)
from tests.test_catalog_revisions import (
    HAS_TEXT_BASE,
    SIDERITE,
    _read_at_revision,
    _tree,
    add_request,
    apply,
    assert_history,
    record_text_failures,
    repo,  # noqa: F401 - fixture
    undo,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
AUDITS = ROOT / "peritheos/data/audits"
KEY = tool.MATERIAL_KEY
BLOCK = {"notes": "Original material wording.", "superseded_on": "2026-11-28"}


def standalone(**extra):
    return {
        "format_version": 5,
        "identifier": "example",
        "formula": "MgO",
        "notes": "Rewritten material wording.",
        **extra,
        "eos_records": [
            {
                "identifier": "example_bm3",
                "label": "Example (2020), BM3",
                "model": "BM3",
                "temperature_ref": 300,
                "parameters": {"V0": 74.7, "K0": 160, "K0_prime": 4},
                "reference": "Fictional example; no experimental validation.",
                "evidence": {"scientific_validation": {"status": "deferred"}},
            }
        ],
    }


def test_material_block_is_split_merged_formatted_and_validated():
    document = standalone(superseded_material_text=copy.deepcopy(BLOCK))
    records, sidecar = split_document(document, add_check=False)
    assert MATERIAL_AUDIT_KEY not in records
    assert sidecar[MATERIAL_AUDIT_KEY] == BLOCK
    # The block follows "material" in a sidecar and is last in a full file.
    assert list(json.loads(format_audit(sidecar))) == [
        "format",
        "format_version",
        "material",
        MATERIAL_AUDIT_KEY,
        "records",
    ]
    merged = merge_audit(records, sidecar)
    assert merged[MATERIAL_AUDIT_KEY] == BLOCK
    assert list(json.loads(format_eosmat(merged)))[-1] == MATERIAL_AUDIT_KEY
    for version in (3, 4):
        assert convert_eosmat(merged, version=version)[MATERIAL_AUDIT_KEY] == BLOCK
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v5.schema.json").read_text(encoding="utf-8")
    )
    jsonschema.Draft202012Validator(schema).validate(merged)
    with pytest.raises(EosmatError, match="both inline"):
        merge_audit({**records, MATERIAL_AUDIT_KEY: BLOCK}, sidecar)


@pytest.mark.parametrize(
    ("block", "message"),
    [
        ({"notes": "x"}, "superseded_on must be an ISO date"),
        ({"notes": "x", "label": "y", "superseded_on": "2026-11-28"}, "unknown"),
        ({"superseded_on": "2026-11-28"}, "names no superseded field"),
        ({"notes": 3, "superseded_on": "2026-11-28"}, "notes must be a string"),
    ],
)
def test_invalid_material_block_is_rejected_inline_and_in_a_sidecar(block, message):
    with pytest.raises(EosmatError, match=message):
        convert_eosmat(standalone(superseded_material_text=block), version=3)
    records, sidecar = split_document(standalone(), add_check=False)
    sidecar[MATERIAL_AUDIT_KEY] = block
    with pytest.raises(EosmatError, match=message):
        merge_audit(records, sidecar)


def _catalog_copy(root, materials):
    """A minimal repository root with some materials and their bindings."""
    receipt = json.loads(RECEIPT.read_bytes())
    receipt["materials"] = {m: receipt["materials"][m] for m in materials}
    receipt.pop("superseded", None)
    # The statement is kept while some binding holds a material block.
    if not any(b.get(MATERIAL_AUDIT_KEY) for b in receipt["materials"].values()):
        receipt.pop("material_text_revision", None)
    for directory in ("materials", "audits"):
        (root / "peritheos/data" / directory).mkdir(parents=True)
    for material in materials:
        shutil.copy2(
            MATERIALS / f"{material}.eosmat", root / "peritheos/data/materials"
        )
        shutil.copy2(AUDITS / f"{material}.audit.json", root / "peritheos/data/audits")
    target = root / RECEIPT.relative_to(ROOT)
    target.parent.mkdir(parents=True)
    target.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
    # Notes already revised in the bundled catalog go back to their stored
    # original first, so the copy starts without a material block.
    plan = tool.Plan()
    for material in materials:
        block = _sidecar(root, material).get(MATERIAL_AUDIT_KEY)
        if block is not None:
            plan.changes[material] = {KEY: {"notes": block["notes"]}}
    tool.apply_plan(plan, "2026-11-28", root)
    return root


def _document(root, material):
    return json.loads(
        (root / f"peritheos/data/materials/{material}.eosmat").read_bytes()
    )


def _sidecar(root, material):
    return json.loads(
        (root / f"peritheos/data/audits/{material}.audit.json").read_bytes()
    )


def _change(root, material, notes):
    current = _document(root, material).get("notes")
    return {material: {KEY: {"before": {"notes": current}, "after": {"notes": notes}}}}


def _rewrite(root, material, notes, date, allow=()):
    plan = tool.plan_changes(_change(root, material, notes), {}, allow=allow, root=root)
    assert plan.refusals == [] and plan.changes == {material: {KEY: {"notes": notes}}}
    tool.apply_plan(plan, date, root)


def _files(root, materials):
    paths = [RECEIPT.relative_to(ROOT)]
    for material in materials:
        paths += [
            Path(f"peritheos/data/materials/{material}.eosmat"),
            Path(f"peritheos/data/audits/{material}.audit.json"),
        ]
    return {path: (root / path).read_bytes() for path in paths}


def test_rewrite_keeps_original_recovers_pre_split_and_reverts_exactly(tmp_path):
    root = _catalog_copy(tmp_path, ["akimotoite", "larnite"])
    before = _files(root, ["akimotoite", "larnite"])
    original = _document(root, "akimotoite")
    pre_split = pre_split_bytes("akimotoite", root=root)
    assert pre_split == pre_split_bytes("akimotoite")
    first = "Ilmenite-type MgSiO3; the hexagonal cell has Z=6."
    _rewrite(root, "akimotoite", first, "2026-11-28")
    document = _document(root, "akimotoite")
    assert document["notes"] == first
    # Only the material notes changed in the record file; the stored original
    # is verbatim and the sidecar entries are untouched.
    assert {k: v for k, v in document.items() if k != "notes"} == {
        k: v for k, v in original.items() if k != "notes"
    }
    sidecar = _sidecar(root, "akimotoite")
    assert sidecar[MATERIAL_AUDIT_KEY] == {
        "notes": original["notes"],
        "superseded_on": "2026-11-28",
    }
    assert (
        sidecar["records"]
        == json.loads(before[Path("peritheos/data/audits/akimotoite.audit.json")])[
            "records"
        ]
    )
    for path in (
        "peritheos/data/materials/akimotoite.eosmat",
        "peritheos/data/audits/akimotoite.audit.json",
    ):
        text = (root / path).read_text()
        formatter = format_audit if "audits" in path else format_eosmat
        assert formatter(json.loads(text)) == text
    # The receipt binds the block; the pre-split file recovers byte for byte.
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    assert receipt["materials"]["akimotoite"][MATERIAL_AUDIT_KEY] is True
    assert "material_text_revision" in receipt
    assert pre_split_bytes("akimotoite", root=root) == pre_split
    assert verify_catalog(root) == 2
    merged = merge_audit(document, sidecar)
    view = pre_split_view(merged, root)
    assert view["notes"] == original["notes"] and MATERIAL_AUDIT_KEY not in view

    if HAS_TEXT_BASE:
        checked, failures = tool.verify(root=root, read=_read_at_revision)
        assert failures == [] and checked >= 1
    # A second rewrite keeps the first original and its date.
    _rewrite(root, "akimotoite", "Ilmenite-type MgSiO3 with Z=6.", "2026-11-30")
    assert _sidecar(root, "akimotoite")[MATERIAL_AUDIT_KEY] == {
        "notes": original["notes"],
        "superseded_on": "2026-11-28",
    }
    assert pre_split_bytes("akimotoite", root=root) == pre_split
    # Returning to the original removes the block and the binding flag.
    _rewrite(root, "akimotoite", original["notes"], "2026-12-01")
    assert _files(root, ["akimotoite", "larnite"]) == before
    # An edit without a stored original is reported.
    document = _document(root, "larnite")
    document["notes"] = "Edited by hand."
    (root / "peritheos/data/materials/larnite.eosmat").write_text(
        format_eosmat(document)
    )
    if HAS_TEXT_BASE:
        _, failures = tool.verify(root=root, read=_read_at_revision)
        assert failures == [
            "larnite: unrecorded change of material notes differ from c26e8d81"
        ]


def test_material_notes_may_be_added_and_removed(tmp_path):
    root = _catalog_copy(tmp_path, ["larnite"])
    before = _files(root, ["larnite"])
    pre_split = pre_split_bytes("larnite", root=root)
    original = _document(root, "larnite")["notes"]
    _rewrite(root, "larnite", None, "2026-11-28")
    assert "notes" not in _document(root, "larnite")
    assert _sidecar(root, "larnite")[MATERIAL_AUDIT_KEY]["notes"] == original
    assert pre_split_bytes("larnite", root=root) == pre_split
    _rewrite(root, "larnite", original, "2026-11-29")
    assert _files(root, ["larnite"]) == before


def _refusals(root, batches, allow=()):
    plan = tool.plan_changes(batches, {}, allow=allow, root=root)
    return {(r.record, r.reason) for r in plan.refusals}, plan


@pytest.mark.parametrize(
    ("notes", "guard"),
    [
        ("Ilmenite-type MgSiO3 with Z=7 and 123.45 GPa.", "numbers"),
        ("Ilmenite-type MgSiO3 with Z=6 from Horiuchi (1983).", "numbers"),
        ("x" * 301, "material_notes_length"),
        ("See docs/eos-labels.md.", "forbidden_text"),
        ("LitCurate import.", "forbidden_text"),
    ],
)
def test_each_material_guard_refuses_and_allow_overrides_it(tmp_path, notes, guard):
    root = _catalog_copy(tmp_path, ["akimotoite"])
    batches = _change(root, "akimotoite", notes)
    refusals, plan = _refusals(root, batches)
    assert (KEY, guard) in refusals and not plan.changes
    refusals, plan = _refusals(root, batches, allow=[f"akimotoite/{KEY}:{guard}"])
    assert (KEY, guard) not in refusals
    assert any(guard in line for line in plan.allowed)


def test_material_numbers_come_from_the_record_file_and_stored_notes(tmp_path):
    root = _catalog_copy(tmp_path, ["akimotoite"])
    # 1982 and 788 are in the current notes, V0 values are in the records.
    record = _document(root, "akimotoite")["eos_records"][0]
    volume = record["parameters"]["V0"]
    volume = volume["value"] if isinstance(volume, dict) else volume
    notes = f"Structure from Horiuchi et al. (1982), p.788; V0 = {volume} Å³."
    refusals, plan = _refusals(root, _change(root, "akimotoite", notes))
    assert refusals == set() and plan.changes


def test_material_input_is_checked_against_the_file(tmp_path):
    root = _catalog_copy(tmp_path, ["akimotoite"])
    batches = _change(root, "akimotoite", "Ilmenite-type MgSiO3.")
    batches["akimotoite"][KEY]["before"]["notes"] = "Different text."
    refusals, plan = _refusals(root, batches)
    assert refusals == {(KEY, "before_mismatch")} and not plan.changes
    path = tmp_path / "batch.json"
    path.write_text(
        json.dumps(
            {"akimotoite": {KEY: {"before": {"label": "x"}, "after": {"label": "y"}}}}
        )
    )
    with pytest.raises(ValueError, match="notes only"):
        tool.read_batches([path], tool.Plan())
    # A sidecar block the receipt does not bind is refused.
    sidecar = _sidecar(root, "akimotoite")
    sidecar[MATERIAL_AUDIT_KEY] = dict(BLOCK)
    sidecar = {k: sidecar[k] for k in ("format", "format_version", "material")} | {
        MATERIAL_AUDIT_KEY: BLOCK,
        "records": sidecar["records"],
    }
    (root / "peritheos/data/audits/akimotoite.audit.json").write_text(
        format_audit(sidecar)
    )
    refusals, _ = _refusals(root, _change(root, "akimotoite", "Ilmenite-type."))
    assert refusals == {(KEY, "predates_split")}


def _all_verifies_pass(root):
    assert verify_catalog(root) == len(
        json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())["materials"]
    )
    assert record_text_failures(root) == []
    assert verify_check_text(root)[1] == []
    assert revision_tool.verify(root) == []


SIDERITE_NOTES = (
    "Diffraction model combining the EOS sample lattice with occupied sites of "
    "matching composition; lattice volume 293.4853 A^3."
)


def test_notes_rewrite_then_catalog_revision_then_undo(repo):  # noqa: F811
    original = _tree(repo)
    _rewrite(repo, SIDERITE, SIDERITE_NOTES, "2026-11-28")
    _all_verifies_pass(repo)
    request = add_request(repo, "t-notes-first")
    apply(repo, request)
    assert_history(repo, [SIDERITE])
    _all_verifies_pass(repo)
    assert _document(repo, SIDERITE)["notes"] == SIDERITE_NOTES
    reference = json.loads(pre_split_bytes(SIDERITE, root=repo))
    assert reference["notes"] == json.loads(pre_split_bytes(SIDERITE))["notes"]
    undo(repo, request["id"])
    assert_history(repo, [])
    _all_verifies_pass(repo)
    previous = json.loads(original[f"peritheos/data/materials/{SIDERITE}.eosmat"])
    # The original exceeds the length limit; a revert names the override.
    _rewrite(
        repo,
        SIDERITE,
        previous["notes"],
        "2026-11-29",
        allow=[f"{SIDERITE}/{KEY}:material_notes_length"],
    )
    assert _tree(repo) == original


def test_catalog_revision_then_notes_rewrite_then_undo(repo):  # noqa: F811
    original = _tree(repo)
    request = add_request(repo, "t-notes-after")
    apply(repo, request)
    _rewrite(repo, SIDERITE, SIDERITE_NOTES, "2026-11-28")
    assert_history(repo, [SIDERITE])
    _all_verifies_pass(repo)
    # Undo goes through the later notes rewrite and keeps it.
    plan = undo(repo, request["id"])
    assert any("later text revisions" in line for line in plan.shown)
    assert _document(repo, SIDERITE)["notes"] == SIDERITE_NOTES
    assert_history(repo, [])
    _all_verifies_pass(repo)
    previous = json.loads(original[f"peritheos/data/materials/{SIDERITE}.eosmat"])
    # The original exceeds the length limit; a revert names the override.
    _rewrite(
        repo,
        SIDERITE,
        previous["notes"],
        "2026-11-29",
        allow=[f"{SIDERITE}/{KEY}:material_notes_length"],
    )
    assert _tree(repo) == original


def test_a_lazily_loaded_material_exports_its_superseded_material_text():
    """The record file alone and the complete document export the same content."""
    import json

    from peritheos import Material, get_material_document
    from peritheos._eosmat_audit import MATERIAL_AUDIT_KEY
    from peritheos.materials import _material_from_eosmat

    audits = Path(__file__).resolve().parents[1] / "peritheos/data/audits"
    rewritten = [
        path.name.removesuffix(".audit.json")
        for path in sorted(audits.glob("*.audit.json"))
        if MATERIAL_AUDIT_KEY in json.loads(path.read_bytes())
    ]
    assert rewritten
    for identifier in rewritten[:5]:
        record_file = get_material_document(
            identifier, format_version=None, include_audit=False
        )
        assert MATERIAL_AUDIT_KEY not in record_file
        lazy = _material_from_eosmat(
            record_file, require_primary_validation=False, bundled_audit=True
        )
        eager = Material.from_eosmat(
            get_material_document(identifier, format_version=5),
            require_primary_validation=False,
        )
        exported = lazy.to_eosmat()
        assert exported == eager.to_eosmat()
        assert MATERIAL_AUDIT_KEY in lazy.to_eosmat(format_version=5)
