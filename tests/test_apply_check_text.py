"""``scripts/apply_check_text.py`` on temporary copies of bundled materials.

Covers every guard, the dry run, a round trip on two materials (apply, the
lossless audit-split reverse, removal restoring the verbatim ``finding``), a
coded usage recommendation set, changed and removed on two materials,
``split_eosmat_audits.py --add`` with the fields present and the
``check_record_text.py`` cases. Bundled files are never written.
"""

import hashlib
import json
import shutil
from pathlib import Path

import pytest

from peritheos import save_eosmat
from peritheos._eosmat_audit import merge_audit
from peritheos._eosmat_format import format_audit
from scripts import apply_check_text as tool
from scripts.apply_record_text import Plan
from scripts.check_record_text import current_cases, record_cases
from scripts.split_eosmat_audits import (
    RECEIPT,
    check_text,
    pre_split_bytes,
    pre_split_view,
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
    receipt_path = root / RECEIPT.relative_to(ROOT)
    with pytest.raises(ValueError, match="only for new materials"):
        split_added_materials(["mgal2o4_cati2o4"], root=root, check_snapshot=False)
    # As a new material (no binding yet) the same file splits and recovers.
    receipt = json.loads(receipt_path.read_bytes())
    del receipt["materials"]["mgal2o4_cati2o4"]
    receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
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


def _tree(root):
    return {
        path.relative_to(root): path.read_bytes()
        for path in sorted(root.rglob("*"))
        if path.is_file()
    }


def test_a_failure_on_a_later_material_writes_nothing(root):
    # coesite_iv is applied first; a broken binding of mgal2o4_cati2o4 makes
    # its pre-split file unrecoverable. Nothing may be written: neither
    # material's files nor the receipt.
    receipt_path = root / RECEIPT.relative_to(ROOT)
    receipt = json.loads(receipt_path.read_bytes())
    receipt["materials"]["mgal2o4_cati2o4"]["before_sha256"] = "0" * 64
    receipt_path.write_text(json.dumps(receipt, indent=2, ensure_ascii=False) + "\n")
    before = _tree(root)
    result = plan(BATCH, root)
    assert list(result.changes) == ["coesite_iv", "mgal2o4_cati2o4"]
    with pytest.raises(ValueError, match="predecessor changed"):
        tool.apply_plan(result, root)
    assert _tree(root) == before


def test_dry_run_refuses_a_file_outside_the_canonical_layout(root, capsys):
    audit = root / "peritheos/data/audits/coesite_iv.audit.json"
    audit.write_text(json.dumps(json.loads(audit.read_bytes()), indent=4) + "\n")
    result = plan(BATCH, root)
    assert refusals(result) == {(BYKOVA_REF, "not_canonical")}
    assert list(result.changes) == ["mgal2o4_cati2o4"]


@pytest.mark.parametrize(
    ("note", "unknown"),
    [
        # A correct rounding at the precision written is accepted.
        ("Table 1 prints V0 = 240 Å³ and K0 = 206 GPa.", None),
        ("Table 1 prints V0 = 240.30 Å³.", None),
        # A wrong rounding, a changed sign or an exponent are numbers too.
        ("Table 1 prints V0 = 241 Å³.", "241"),
        ("Table 1 prints V0 = -240.3 Å³.", "-240.3"),
        ("Table 1 prints V0 = −240.3 Å³.", "-240.3"),
        ("Table 1 prints a = 2.403e2 Å³.", None),
        ("Table 1 prints a = 2.41e2 Å³.", "2.41e2"),
        ("Table 1 prints 1e-5 per K.", "1e-5"),
        # A locator's full stop does not hide the number after it.
        ("See p.750 of the paper.", "750"),
        ("See Eq.12 and Fig.3.", "12"),
    ],
)
def test_number_guard_reads_signs_exponents_locators_and_roundings(root, note, unknown):
    result = plan({"mgal2o4_cati2o4": {FUNAMORI: {"source_note": note}}}, root)
    details = [r.detail for r in result.refusals]
    assert details == ([] if unknown is None else [f"not in the record: {unknown}"])


def test_bundled_check_texts_pass_every_guard_without_overrides():
    """Every bundled record_note and warning_reason passes the guards as is."""
    failed = []
    for path in sorted(MATERIALS.glob("*.eosmat")):
        merged = merge_audit(
            json.loads(path.read_bytes()),
            json.loads((AUDITS / f"{path.stem}.audit.json").read_bytes()),
        )
        for record in merged["eos_records"]:
            new = {
                key: check_text(record, field)
                for key, field in tool.INPUT_FIELDS.items()
                if check_text(record, field) is not None
            }
            if new:
                failed += [
                    (record["identifier"], guard) for guard in tool._guards(record, new)
                ]
    assert failed == []


@pytest.mark.parametrize(
    ("material", "identifier"),
    [
        # Both rebuild primary_source_check from the script's own tables.
        ("neon_fcc", "neon_fcc_hemley_1989_bm3_refit"),
        ("aragonite", "aragonite_martinez_1996_bm2_2"),
    ],
)
def test_primary_source_audit_keeps_the_record_note(material, identifier):
    from peritheos import get_material_document
    from scripts.apply_primary_source_audit import audit_record

    document = get_material_document(material)
    record = next(r for r in document["eos_records"] if r["identifier"] == identifier)
    note = record["scientific_validation"]["primary_source_check"]["record_note"]
    rebuilt = audit_record(json.loads(json.dumps(record)), f"{material}.eosmat")
    check = rebuilt["scientific_validation"]["primary_source_check"]
    assert check["record_note"] == note
    assert list(check).index("record_note") == list(check).index("finding") + 1


# --- Coded usage recommendations ------------------------------------------

ARGON = "argon_hcp_wittlinger_1997_bm2_1"  # free-text usage_recommendation
MAO = "ca_perovskite_mao_1989_bm3_2"  # no recommendation, no qualifications
FU = "ca_perovskite_fu_2023_bm3_mgd_refit"  # coded, with a warning_reason
FU_CANDIDATE = "ca_perovskite_fu_2023_candidate_data_unweighted_bm3_mgd_refit"
REASON_A = "The published fit is not reproduced from the digitized Figure 3 data."
REASON_M = "Superseded by later measurements of the same phase."
USAGE_BATCH = {
    "argon_hcp": {
        ARGON: {"usage_recommendation": "qualified_use", "warning_reason": REASON_A}
    },
    "ca_perovskite": {
        MAO: {
            "usage_recommendation": "not_recommended_for_quantitative_use",
            "warning_reason": REASON_M,
        },
        FU: {"usage_recommendation": "qualified_use"},
    },
}
LEDGER = "docs/data/primary-eos-refits.json"
LEDGER_MARKDOWN = "docs/primary-eos-refits.md"
REVISIONS = "docs/data/catalog-revisions.json"


@pytest.fixture
def usage_root(tmp_path):
    """Two split materials as bundled, their receipt and the refit ledger."""
    materials = ["argon_hcp", "ca_perovskite"]
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
    shutil.copy2(ROOT / LEDGER, tmp_path / LEDGER)
    shutil.copy2(ROOT / LEDGER_MARKDOWN, tmp_path / LEDGER_MARKDOWN)
    # The catalog revisions of these materials, which their bindings pin.
    revisions = json.loads((ROOT / REVISIONS).read_bytes())
    revisions["revisions"] = [
        r for r in revisions["revisions"] if r["material"] in materials
    ]
    (tmp_path / REVISIONS).write_text(
        json.dumps(revisions, indent=2, ensure_ascii=False) + "\n"
    )
    # Start argon from its original free-text recommendation, whatever code the
    # bundled record carries now: the tool itself returns the copy to the
    # stored original, and no usage binding remains.
    removal = {
        "argon_hcp": {ARGON: {"usage_recommendation": None, "warning_reason": None}}
    }
    if "check_text" in receipt["materials"]["argon_hcp"]:
        result = plan(removal, tmp_path)
        assert not result.refusals, [r.line() for r in result.refusals]
        tool.apply_plan(result, tmp_path)
    assert tool.verify(tmp_path)[1] == []
    binding = json.loads(target.read_bytes())["materials"]["argon_hcp"]
    assert "check_text" not in binding
    argon = evidence(tmp_path, "argon_hcp", ARGON)
    prose = argon["usage_recommendation"]
    assert not prose.islower() and "8.5 GPa" in prose
    assert "recommendation_reason" not in argon and "warning_reason" not in argon
    assert checks(tmp_path, "argon_hcp")[ARGON]["warning"] is None
    return tmp_path


def _usage_tree(root):
    names = [
        f"peritheos/data/{kind}/{m}{suffix}"
        for m in ("argon_hcp", "ca_perovskite")
        for kind, suffix in (("materials", ".eosmat"), ("audits", ".audit.json"))
    ]
    names += [str(RECEIPT.relative_to(ROOT)), LEDGER, LEDGER_MARKDOWN]
    return {name: (root / name).read_bytes() for name in names}


def _merged(root, material):
    record_file, sidecar = (json.loads(b) for b in files(root, material))
    return {
        r["identifier"]: r for r in merge_audit(record_file, sidecar)["eos_records"]
    }


def _historical_checks(root, material):
    """The receipts that failed when a record was re-split with --add."""
    from scripts.migrate_eosmat_v5 import prior_material_bytes

    recovered = pre_split_bytes(material, root=root)
    assert recovered == pre_split_bytes(material)
    assert prior_material_bytes(material, root=root) == prior_material_bytes(material)
    policy = json.loads(
        (ROOT / "docs/data/eos-structure-rights-policy.json").read_bytes()
    )
    relative = f"peritheos/data/materials/{material}.eosmat"
    assert hashlib.sha256(recovered).hexdigest() == policy["material_sha256"][relative]


def _ledger_consistent(root):
    """Ledger copies, the rendered page and every assessment agree with the files."""
    from scripts.standardize_reconstruction_metadata import assessment_for_record
    from scripts.validate_primary_eos_refits import render_markdown

    ledger = json.loads((root / LEDGER).read_bytes())
    rows = {row["record_identifier"]: row for row in ledger["records"]}
    for material in ("argon_hcp", "ca_perovskite"):
        for identifier, record in _merged(root, material).items():
            validation = record["evidence"]["scientific_validation"]
            for key in ("usage_recommendation", "recommendation_reason"):
                assert rows[identifier].get(key) == validation.get(key)
            assessment = record["evidence"]["reproduction"]["assessment"]
            assert assessment == assessment_for_record(record, rows[identifier])
    assert (root / LEDGER_MARKDOWN).read_text(encoding="utf-8") == render_markdown(
        ledger
    )


def test_usage_recommendation_set_change_and_removal(usage_root):
    root = usage_root
    before = _usage_tree(root)
    old = {m: _merged(root, m) for m in ("argon_hcp", "ca_perovskite")}
    prose = old["argon_hcp"][ARGON]["evidence"]["scientific_validation"][
        "usage_recommendation"
    ]
    assert not prose.islower() and "8.5 GPa" in prose
    result = plan(USAGE_BATCH, root)
    assert not result.refusals, [r.line() for r in result.refusals]
    warnings = tool.apply_plan(result, root)

    # The warning is generated from the code; nothing else in check changes.
    assert warnings == {
        "argon_hcp": {ARGON: "Qualified use only; see the audit."},
        "ca_perovskite": {
            MAO: "Not recommended for quantitative use.",
            FU: "Qualified use only; see the audit. Never chosen automatically; "
            "select it explicitly.",
        },
    }
    for material, records in USAGE_BATCH.items():
        new_checks = checks(root, material)
        old_checks = {k: v["check"] for k, v in old[material].items()}
        for record in records:
            assert new_checks[record]["warning"] == warnings[material][record]
        for record, check in old_checks.items():
            kept = ("warning", "warning_reason")
            assert {k: v for k, v in new_checks[record].items() if k not in kept} == {
                k: v for k, v in check.items() if k not in kept
            }
    assert checks(root, "argon_hcp")[ARGON]["warning_reason"] == REASON_A

    # Free text replaced by a code is kept verbatim as recommendation_reason and
    # in the qualifications, as in records that already carry a code.
    argon = evidence(root, "argon_hcp", ARGON)
    assert argon["usage_recommendation"] == "qualified_use"
    assert argon["recommendation_reason"] == prose
    assert list(argon).index("recommendation_reason") == (
        list(argon).index("warning_reason") + 1
    )
    sidecar = json.loads(files(root, "argon_hcp")[1])
    assessment = sidecar["records"][ARGON]["evidence"]["reproduction"]["assessment"]
    assert assessment["qualifications"] == [
        prose,
        "Usage recommendation: qualified_use",
    ]
    mao = _merged(root, "ca_perovskite")[MAO]["evidence"]
    assert mao["scientific_validation"]["usage_recommendation"] == (
        "not_recommended_for_quantitative_use"
    )
    assert "recommendation_reason" not in mao["scientific_validation"]
    assert mao["reproduction"]["assessment"]["qualifications"] == [
        "Usage recommendation: not_recommended_for_quantitative_use"
    ]
    fu = _merged(root, "ca_perovskite")[FU]["evidence"]
    fu_old = old["ca_perovskite"][FU]["evidence"]
    assert (
        fu["scientific_validation"]["recommendation_reason"]
        == (fu_old["scientific_validation"]["recommendation_reason"])
    )
    assert fu["reproduction"]["assessment"]["qualifications"] == [
        fu_old["scientific_validation"]["recommendation_reason"],
        "Usage recommendation: qualified_use",
    ]

    # The receipt binds the originals; every historical receipt still verifies.
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    assert receipt["materials"]["argon_hcp"]["check_text"] == [
        [
            0,
            {
                "usage_recommendation": prose,
                "recommendation_reason": None,
                "assessment_qualifications": [f"Usage recommendation: {prose}"],
                "warning_reason": None,
            },
        ]
    ]
    assert (
        receipt["materials"]["argon_hcp"]["before_sha256"]
        == (json.loads(RECEIPT.read_bytes())["materials"]["argon_hcp"]["before_sha256"])
    )
    for material in ("argon_hcp", "ca_perovskite"):
        _historical_checks(root, material)
        # The historical-fingerprint projection restores the originals too.
        current = merge_audit(*(json.loads(b) for b in files(root, material)))
        bundled = merge_audit(*(json.loads(b) for b in files(ROOT, material)))
        assert pre_split_view(current, root) == pre_split_view(bundled)
    assert verify_catalog(root) == 2
    assert tool.verify(root)[1] == []
    _ledger_consistent(root)
    row = next(
        r
        for r in json.loads((root / LEDGER).read_bytes())["records"]
        if r["record_identifier"] == ARGON
    )
    assert row["usage_recommendation"] == "qualified_use"
    assert row["recommendation_reason"] == row["qualification"] == prose

    # A second change keeps the first original and the moved text.
    second = {
        "argon_hcp": {
            ARGON: {"usage_recommendation": "not_recommended_for_quantitative_use"}
        }
    }
    result = plan(second, root)
    assert not result.refusals, [r.line() for r in result.refusals]
    tool.apply_plan(result, root)
    assert checks(root, "argon_hcp")[ARGON]["warning"] == (
        "Not recommended for quantitative use."
    )
    assert evidence(root, "argon_hcp", ARGON)["recommendation_reason"] == prose
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    assert receipt["materials"]["argon_hcp"]["check_text"][0][1][
        "usage_recommendation"
    ] == (prose)
    _historical_checks(root, "argon_hcp")
    _ledger_consistent(root)

    # null returns to the stored original; with the reasons removed, every file
    # is the original byte for byte.
    removal = {
        "argon_hcp": {ARGON: {"usage_recommendation": None, "warning_reason": None}},
        "ca_perovskite": {
            MAO: {"usage_recommendation": None, "warning_reason": None},
            FU: {"usage_recommendation": None},
        },
    }
    result = plan(removal, root)
    assert not result.refusals, [r.line() for r in result.refusals]
    tool.apply_plan(result, root)
    assert _usage_tree(root) == before
    assert tool.verify(root)[1] == []


def test_setting_the_original_code_again_restores_everything(usage_root):
    before = _usage_tree(usage_root)
    batch = {"ca_perovskite": {FU: {"usage_recommendation": "qualified_use"}}}
    tool.apply_plan(plan(batch, usage_root), usage_root)
    assert _usage_tree(usage_root) != before
    back = {
        "ca_perovskite": {
            FU: {"usage_recommendation": "not_recommended_for_quantitative_use"}
        }
    }
    result = plan(back, usage_root)
    assert not result.refusals
    tool.apply_plan(result, usage_root)
    assert _usage_tree(usage_root) == before


def test_usage_recommendation_refusals(usage_root):
    root = usage_root
    # Only the coded values the library turns into a warning are accepted.
    for value in ("Use with care below 8.5 GPa.", "use_with_care", 3, ""):
        batch = {"argon_hcp": {ARGON: {"usage_recommendation": value}}}
        assert refusals(plan(batch, root)) == {(ARGON, "not_coded")}
    # A first warning needs a warning_reason; --allow cannot override that.
    batch = {"ca_perovskite": {MAO: {"usage_recommendation": "qualified_use"}}}
    result = plan(batch, root, allow=[f"{MAO}:{g}" for g in tool.GUARDS])
    assert refusals(result) == {(MAO, "reason_required")}
    assert not result.refusals[0].guard and not result.changes
    # Returning to no warning while a warning_reason stays is refused.
    tool.apply_plan(
        plan(
            {
                "ca_perovskite": {
                    MAO: {
                        "usage_recommendation": "qualified_use",
                        "warning_reason": REASON_M,
                    }
                }
            },
            root,
        ),
        root,
    )
    batch = {"ca_perovskite": {MAO: {"usage_recommendation": None}}}
    assert refusals(plan(batch, root)) == {(MAO, "no_warning")}
    # The reason's text guards apply as before.
    batch = {
        "argon_hcp": {
            ARGON: {
                "usage_recommendation": "qualified_use",
                "warning_reason": "See docs/argon.md.",
            }
        }
    }
    assert refusals(plan(batch, root)) == {(ARGON, "forbidden_text")}
    # A change that would alter the assessment beyond its qualifications.
    batch = {"ca_perovskite": {FU_CANDIDATE: {"usage_recommendation": "qualified_use"}}}
    batch["ca_perovskite"][FU_CANDIDATE]["warning_reason"] = REASON_M
    result = plan(batch, root)
    assert refusals(result) == {(FU_CANDIDATE, "beyond_fields")}
    assert "basis" in result.refusals[0].detail
    # A different recommendation_reason would hide the free text: refused.
    # (Written consistently into sidecar, ledger row and assessment.)
    from scripts.standardize_reconstruction_metadata import assessment_for_record

    path = root / "peritheos/data/audits/argon_hcp.audit.json"
    sidecar = json.loads(path.read_bytes())
    ledger = json.loads((root / LEDGER).read_bytes())
    row = next(r for r in ledger["records"] if r["record_identifier"] == ARGON)
    row["recommendation_reason"] = row["qualification"] = "Another reason."
    (root / LEDGER).write_text(json.dumps(ledger, indent=2, ensure_ascii=False) + "\n")
    entry = sidecar["records"][ARGON]
    entry["evidence"]["scientific_validation"]["recommendation_reason"] = (
        "Another reason."
    )
    record = {**_merged(root, "argon_hcp")[ARGON], "evidence": entry["evidence"]}
    entry["evidence"]["reproduction"]["assessment"]["qualifications"] = (
        assessment_for_record(record, row)["qualifications"]
    )
    path.write_text(format_audit(sidecar), encoding="utf-8")
    batch = {"argon_hcp": {ARGON: USAGE_BATCH["argon_hcp"][ARGON]}}
    assert refusals(plan(batch, root)) == {(ARGON, "recommendation_reason_present")}


def test_dry_run_prints_the_resulting_warning(tmp_path, capsys):
    batch = tmp_path / "batch.json"
    batch.write_text(json.dumps({"argon_hcp": USAGE_BATCH["argon_hcp"]}))
    bundled = _usage_tree(ROOT)
    assert tool.main([str(batch), "--check"]) == 0
    out = capsys.readouterr().out
    assert (
        f"WOULD-CHANGE argon_hcp/{ARGON}: warning_reason {len(REASON_A)} characters, "
        'usage_recommendation qualified_use; warning: "Qualified use only; see the '
        'audit."'
    ) in out
    assert _usage_tree(ROOT) == bundled


def test_reference_lists_every_accepted_code_with_its_warning():
    from peritheos._eosmat_audit import USAGE_RECOMMENDATION_WARNINGS

    text = (ROOT / "docs/eosmat-reference.md").read_text(encoding="utf-8")
    section = text.split("### Giving a bundled record a warning", 1)[1]
    rows = {
        cells[0].strip("` "): cells[1].strip()
        for line in section.splitlines()
        if line.startswith("| `")
        for cells in [line.strip("|").split("|")]
    }
    rows.pop("usage_recommendation")  # the header
    assert rows == USAGE_RECOMMENDATION_WARNINGS
    assert tool.CODED_RECOMMENDATIONS == tuple(USAGE_RECOMMENDATION_WARNINGS)
