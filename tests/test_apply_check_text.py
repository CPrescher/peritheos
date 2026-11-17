"""``scripts/apply_check_text.py`` on temporary copies of bundled materials.

Covers every guard, the dry run, a round trip on two materials (apply, the
lossless audit-split reverse, removal restoring the verbatim ``finding``),
``split_eosmat_audits.py --add`` with the fields present and the
``check_record_text.py`` cases. Bundled files are never written.
"""

import json
import shutil
from pathlib import Path

import pytest

from peritheos import save_eosmat
from peritheos._eosmat_audit import merge_audit
from scripts import apply_check_text as tool
from scripts.apply_record_text import Plan
from scripts.check_record_text import current_cases, record_cases
from scripts.split_eosmat_audits import (
    RECEIPT,
    pre_split_bytes,
    split_added_materials,
    verify_catalog,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
AUDITS = ROOT / "peritheos/data/audits"

# Two small materials, each with one restricted and one unrestricted record;
# both already carry a bound superseded_record_text revision.
FUNAMORI = "mgal2o4_cati2o4_funamori_1998_bm2_1"  # Not recommended ...
ONO = "mgal2o4_cati2o4_ono_2006_bm2_2"  # no warning
BYKOVA_REF = "coesite_iv_v_bykova_2018_experimental_bm2_reference"  # warning
BYKOVA_FIT = "coesite_iv_bykova_2018_am05_static_bm3_refit"  # no warning
NOTE_F = "Table 1 prints V0 = 240.3(4) Å³ and K0 = 206(3) GPa with K0' = 4 fixed."
NOTE_O = "The 300 K BM2 values V0 = 238.9(9) Å³ and K0 = 219(6) GPa match the paper."
NOTE_B = "The drawn red curve uses an explicit 32.7 GPa reference pressure."
REASON_F = "K0 rests on one usable compressed cell."
REASON_B = "The author fit details are unresolved."
BATCH = {
    "mgal2o4_cati2o4": {
        FUNAMORI: {"source_note": NOTE_F, "warning_reason": REASON_F, "x": 1},
        ONO: {"source_note": NOTE_O},
    },
    "coesite_iv": {BYKOVA_REF: {"source_note": NOTE_B, "warning_reason": REASON_B}},
}


@pytest.fixture
def root(tmp_path):
    """A minimal repository root: two split materials and their receipt."""
    materials = ["coesite_iv", "mgal2o4_cati2o4"]
    receipt = json.loads(RECEIPT.read_bytes())
    receipt["materials"] = {m: receipt["materials"][m] for m in materials}
    receipt.pop("superseded", None)
    (tmp_path / "peritheos/data/materials").mkdir(parents=True)
    (tmp_path / "peritheos/data/audits").mkdir(parents=True)
    for material in materials:
        shutil.copy2(
            MATERIALS / f"{material}.eosmat", tmp_path / "peritheos/data/materials"
        )
        shutil.copy2(
            AUDITS / f"{material}.audit.json", tmp_path / "peritheos/data/audits"
        )
    target = tmp_path / RECEIPT.relative_to(ROOT)
    target.parent.mkdir(parents=True)
    target.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
    # Start from the uncurated state: remove the bundled curated texts, so
    # each check shows the verbatim finding and no check_text is bound.
    removal = {
        material: {
            record: {"source_note": None, "warning_reason": None}
            for record in json.loads(
                (tmp_path / f"peritheos/data/audits/{material}.audit.json").read_bytes()
            )["records"]
        }
        for material in materials
    }
    result = plan(removal, tmp_path)
    assert not result.refusals, [r.line() for r in result.refusals]
    tool.apply_plan(result, tmp_path)
    assert tool.verify(tmp_path) == (0, [])
    for material in materials:
        binding = json.loads(target.read_bytes())["materials"][material]
        assert "check_text" not in binding
    return tmp_path


def files(root, material):
    return (
        (root / f"peritheos/data/materials/{material}.eosmat").read_bytes(),
        (root / f"peritheos/data/audits/{material}.audit.json").read_bytes(),
    )


def checks(root, material):
    document = json.loads(files(root, material)[0])
    return {r["identifier"]: r["check"] for r in document["eos_records"]}


def evidence(root, material, record):
    document, sidecar = (json.loads(b) for b in files(root, material))
    return sidecar["records"][record]["evidence"]["scientific_validation"]


def plan(batch, root, allow=()):
    result = Plan()
    tool.plan_changes(batch, allow=allow, root=root, plan=result)
    return result


def refusals(result):
    return {(r.record, r.reason) for r in result.refusals}


def test_round_trip_reverse_and_removal_on_two_materials(root):
    before = {m: files(root, m) for m in ("coesite_iv", "mgal2o4_cati2o4")}
    findings = {
        record: evidence(root, material, record)["primary_source_check"]["finding"]
        for material, records in BATCH.items()
        for record in records
    }
    old_checks = {m: checks(root, m) for m in before}
    result = plan(BATCH, root)
    assert not result.refusals, [r.line() for r in result.refusals]
    tool.apply_plan(result, root)

    # Sidecar fields written; finding unchanged; check regenerated.
    for material, records in BATCH.items():
        new_checks = checks(root, material)
        for record, entry in records.items():
            validation = evidence(root, material, record)
            source = validation["primary_source_check"]
            assert source["finding"] == findings[record]
            assert source["record_note"] == entry["source_note"]
            assert (
                list(source).index("record_note") == list(source).index("finding") + 1
            )
            assert new_checks[record]["source_note"] == entry["source_note"]
            if "warning_reason" in entry:
                assert validation["warning_reason"] == entry["warning_reason"]
                assert new_checks[record]["warning_reason"] == entry["warning_reason"]
                assert (
                    new_checks[record]["warning"]
                    == old_checks[material][record]["warning"]
                )
                keys = list(new_checks[record])
                assert keys[keys.index("warning") + 1] == "warning_reason"
            # Nothing else in check changed.
            assert {
                k: v
                for k, v in new_checks[record].items()
                if k not in ("source_note", "warning_reason")
            } == {
                k: v
                for k, v in old_checks[material][record].items()
                if k not in ("source_note", "warning_reason")
            }
        # Untouched records keep their check.
        for record, check in old_checks[material].items():
            if record not in records:
                assert new_checks[record] == check

    # The reverse recovers every pre-split file byte for byte.
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    binding = receipt["materials"]["mgal2o4_cati2o4"]
    assert binding["check_text"] == [
        [0, {"record_note": None, "warning_reason": None}],
        [1, {"record_note": None}],
    ]
    assert "check_text_revision" in receipt
    for material in before:
        assert pre_split_bytes(material, root=root) == pre_split_bytes(material)
    assert verify_catalog(root) == 2
    assert tool.verify(root) == (3, [])

    # A second change keeps the first original; returning to it drops the entry.
    second = {"coesite_iv": {BYKOVA_REF: {"source_note": NOTE_B + " Checked."}}}
    tool.apply_plan(plan(second, root), root)
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    assert receipt["materials"]["coesite_iv"]["check_text"] == [
        [1, {"record_note": None, "warning_reason": None}]
    ]
    assert pre_split_bytes("coesite_iv", root=root) == pre_split_bytes("coesite_iv")

    # Removal (null) restores the verbatim finding and the original files.
    removal = {
        material: {
            record: {"source_note": None, "warning_reason": None} for record in records
        }
        for material, records in BATCH.items()
    }
    result = plan(removal, root)
    assert not result.refusals
    tool.apply_plan(result, root)
    for material in before:
        assert files(root, material) == before[material]
        assert checks(root, material) == old_checks[material]
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    assert "check_text" not in receipt["materials"]["coesite_iv"]
    assert tool.verify(root) == (0, [])


def test_dry_run_writes_nothing(root, tmp_path, capsys):
    snapshot = {m: files(root, m) for m in ("coesite_iv", "mgal2o4_cati2o4")}
    result = plan(BATCH, root)
    assert set(result.changes) == {"coesite_iv", "mgal2o4_cati2o4"}
    assert {m: files(root, m) for m in snapshot} == snapshot
    # The command line dry run on the bundled catalog changes no file.
    batch = tmp_path / "batch.json"
    batch.write_text(json.dumps(BATCH))
    bundled = {m: files(ROOT, m) for m in snapshot}
    receipt = RECEIPT.read_bytes()
    assert tool.main([str(batch), "--check"]) == 0
    assert "Would change 3 records (5 fields) in 2 materials" in capsys.readouterr().out
    assert {m: files(ROOT, m) for m in snapshot} == bundled
    assert RECEIPT.read_bytes() == receipt


@pytest.mark.parametrize(
    ("record", "entry", "guard"),
    [
        (FUNAMORI, {"source_note": "x" * 301}, "source_note_length"),
        (FUNAMORI, {"warning_reason": "x" * 141}, "warning_reason_length"),
        (FUNAMORI, {"source_note": "V0 = 240.5(4) Å³."}, "numbers"),
        (FUNAMORI, {"warning_reason": "Only 17 cells."}, "numbers"),
        (FUNAMORI, {"source_note": "See docs/validation.md."}, "forbidden_text"),
        (FUNAMORI, {"source_note": "LitCurate checked Table 1."}, "forbidden_text"),
        (FUNAMORI, {"warning_reason": "Audit 2026 found a gap."}, "forbidden_text"),
    ],
)
def test_each_guard_refuses_and_can_be_overridden(root, record, entry, guard):
    batch = {"mgal2o4_cati2o4": {record: entry}}
    result = plan(batch, root)
    assert (record, guard) in refusals(result)
    assert not result.changes
    lines = [r.line() for r in result.refusals]
    assert any(f"guard {guard}" in line for line in lines), lines
    other_guards = {r.reason for r in result.refusals} - {guard}
    allow = [f"{record}:{name}" for name in {guard, *other_guards}]
    allowed = plan(batch, root, allow=allow)
    assert not allowed.refusals
    assert allowed.changes == {"mgal2o4_cati2o4": {record: entry}}
    assert any(guard in line for line in allowed.allowed)


def test_numbers_must_occur_in_the_record(root):
    # 240.3, 206 and 58 do not all occur: 240.3(4) and 206(3) are the stored
    # values with their printed errors, single digits are not checked.
    ok = {"source_note": "V0 = 240.3(4) Å³, K0 = 206(3) GPa, K0' = 4."}
    assert not plan({"mgal2o4_cati2o4": {FUNAMORI: ok}}, root).refusals
    # Numbers from the finding, other evidence and stored originals count.
    finding = evidence(root, "coesite_iv", BYKOVA_FIT)["primary_source_check"]
    assert "168.52" in finding["finding"]
    note = {"source_note": "Supplementary Table 11 prints K0 = 168.52 GPa."}
    assert not plan({"coesite_iv": {BYKOVA_FIT: note}}, root).refusals
    bad = {"source_note": "Supplementary Table 11 prints K0 = 168.53 GPa."}
    result = plan({"coesite_iv": {BYKOVA_FIT: bad}}, root)
    assert [r.detail for r in result.refusals] == ["not in the record: 168.53"]


def test_shared_notes_are_a_conscious_exception(root):
    same = "The source table values match."
    batch = {
        "mgal2o4_cati2o4": {
            FUNAMORI: {"source_note": same},
            ONO: {"source_note": same},
        }
    }
    result = plan(batch, root)
    assert refusals(result) == {(FUNAMORI, "shared_note"), (ONO, "shared_note")}
    # Each record needs its own override.
    result = plan(batch, root, allow=[f"{FUNAMORI}:shared_note"])
    assert refusals(result) == {(ONO, "shared_note")}
    assert result.changes == {"mgal2o4_cati2o4": {FUNAMORI: {"source_note": same}}}
    both = [f"{FUNAMORI}:shared_note", f"{ONO}:shared_note"]
    assert not plan(batch, root, allow=both).refusals
    # The existing note (here the verbatim finding) of another record counts.
    ono_finding = evidence(root, "mgal2o4_cati2o4", ONO)["primary_source_check"][
        "finding"
    ]
    copied = {"mgal2o4_cati2o4": {FUNAMORI: {"source_note": ono_finding}}}
    result = plan(copied, root, allow=[f"{FUNAMORI}:numbers"])
    assert refusals(result) == {(FUNAMORI, "shared_note")}


def test_reason_without_warning_and_bad_input_are_refused(root):
    result = plan({"mgal2o4_cati2o4": {ONO: {"warning_reason": "No reason."}}}, root)
    assert refusals(result) == {(ONO, "no_warning")}
    # Not a guard: the library rejects such a record, so --allow cannot help.
    assert not result.refusals[0].guard
    batch = {
        "mgal2o4_cati2o4": {
            FUNAMORI: {"source_note": 3},
            ONO: {"source_note": " padded "},
            "missing": {"source_note": "x"},
        },
        "no_such_material": {"r": {"source_note": "x"}},
    }
    assert refusals(plan(batch, root)) == {
        (FUNAMORI, "not_text"),
        (ONO, "not_text"),
        ("missing", "unknown_record"),
        ("r", "unknown_material"),
    }
    # Removing what is absent, or writing the current text, changes nothing.
    noop = {"mgal2o4_cati2o4": {ONO: {"source_note": None, "warning_reason": None}}}
    result = plan(noop, root)
    assert not result.refusals and not result.changes


def test_added_inline_file_with_the_fields_splits_and_recovers(root, monkeypatch):
    tool.apply_plan(plan(BATCH, root), root)
    applied = files(root, "mgal2o4_cati2o4")
    # A curator writes the complete inline document (record_note and
    # warning_reason in its evidence) and adds it again with --add.
    record_file, sidecar = (json.loads(b) for b in applied)
    complete = merge_audit(record_file, sidecar)
    path = root / "peritheos/data/materials/mgal2o4_cati2o4.eosmat"
    save_eosmat(path, complete)
    inline = path.read_bytes()
    assert split_added_materials(
        ["mgal2o4_cati2o4"], root=root, check_snapshot=False
    ) == ["mgal2o4_cati2o4"]
    assert files(root, "mgal2o4_cati2o4") == applied
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    binding = receipt["materials"]["mgal2o4_cati2o4"]
    # The fields are part of the new pre-split source, not a later revision.
    assert "check_text" not in binding
    assert pre_split_bytes("mgal2o4_cati2o4", root=root) == inline
    assert tool.verify(root) == (3, [])


def test_check_record_text_counts_notes_and_reasons(root):
    # A record over the source-note limit is resolved once a record_note
    # within the limit replaces its finding in check.
    assert "source_note_length" in current_cases(root)[f"mgal2o4_cati2o4/{FUNAMORI}"]
    tool.apply_plan(plan(BATCH, root), root)
    assert f"mgal2o4_cati2o4/{FUNAMORI}" not in current_cases(root)
    check = {"warning": "Do not use.", "warning_reason": "x" * 141}
    assert record_cases({"check": check}) == ["warning_reason_length"]
    check["warning_reason"] = "x" * 140
    assert record_cases({"check": check}) == []
