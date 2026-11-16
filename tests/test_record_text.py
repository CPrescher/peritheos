"""Record text revisions keep every replaced original in the audit sidecar.

Covers the ``superseded_record_text`` audit block, ``scripts/apply_record_text.py``
(round trip, guards, a second rewrite, ``--verify``), the reverse step of the
audit-split recovery and ``scripts/check_record_text.py``.
"""

import copy
import json
import shutil
from pathlib import Path

import pytest

from peritheos import (
    Material,
    convert_eosmat,
    get_eos_record,
    get_material_document,
)
from peritheos._eosmat_audit import (
    AUDIT_RECORD_KEYS,
    SUPERSEDED_TEXT_FIELDS,
    merge_audit,
    split_document,
)
from peritheos._eosmat_format import format_audit, format_eosmat
from peritheos.eosmat import list_material_documents
from peritheos.errors import EosmatError
from scripts import apply_record_text as tool
from scripts.check_record_text import compare, current_cases, load_baseline
from scripts.scientific_source_snapshot import scientific_sha256
from scripts.split_eosmat_audits import (
    RECEIPT,
    merged_document,
    pre_split_bytes,
    record_text,
    set_record_text,
    verify_catalog,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
AUDITS = ROOT / "peritheos/data/audits"
BLOCK = {
    "label": "Example (2020), BM3 (diagnostic only)",
    "notes": "Original wording, kept verbatim.",
    "validity_notes": ["Original range note."],
    "superseded_on": "2026-11-16",
}


def standalone(**extra):
    record = {
        "identifier": "example_bm3",
        "label": "Example (2020), BM3",
        "model": "BM3",
        "temperature_ref": 300,
        "parameters": {"V0": 74.7, "K0": 160, "K0_prime": 4},
        "reference": "Fictional example; no experimental validation.",
        **extra,
    }
    return {
        "format_version": 5,
        "identifier": "example",
        "formula": "MgO",
        "eos_records": [record],
    }


def test_superseded_text_is_an_audit_block_split_merged_and_exported():
    assert "superseded_record_text" in AUDIT_RECORD_KEYS
    document = standalone(superseded_record_text=copy.deepcopy(BLOCK))
    record = Material.from_eosmat(
        document, require_primary_validation=False
    ).eos_records[0]
    assert dict(record.audit)["superseded_record_text"] == BLOCK
    assert record.name == "Example (2020), BM3"
    records, sidecar = split_document(document)
    assert "superseded_record_text" not in records["eos_records"][0]
    assert sidecar["records"]["example_bm3"]["superseded_record_text"] == BLOCK
    merged = merge_audit(records, sidecar)
    assert merged["eos_records"][0]["superseded_record_text"] == BLOCK
    # The block is the last key of a saved record and of a sidecar entry.
    assert list(json.loads(format_eosmat(document))["eos_records"][0])[-1] == (
        "superseded_record_text"
    )
    for version in (3, 4):
        exported = convert_eosmat(merged, version=version)
        assert exported["eos_records"][0]["superseded_record_text"] == BLOCK
    clash = copy.deepcopy(records)
    clash["eos_records"][0]["superseded_record_text"] = BLOCK
    with pytest.raises(EosmatError, match="both inline"):
        merge_audit(clash, sidecar)


@pytest.mark.parametrize(
    ("block", "message"),
    [
        ({"notes": "x"}, "superseded_on must be an ISO date"),
        ({"notes": "x", "superseded_on": "2026-02-30"}, "superseded_on"),
        ({"note": "x", "superseded_on": "2026-11-16"}, "unknown fields"),
        ({"superseded_on": "2026-11-16"}, "names no superseded field"),
        ({"validity_notes": "x", "superseded_on": "2026-11-16"}, "validity_notes"),
        ({"label": 3, "superseded_on": "2026-11-16"}, "label must be a string"),
    ],
)
def test_invalid_superseded_text_is_rejected_inline_and_in_a_sidecar(block, message):
    with pytest.raises(EosmatError, match=message):
        Material.from_eosmat(
            standalone(superseded_record_text=block), require_primary_validation=False
        )
    records, sidecar = split_document(standalone())
    sidecar["records"]["example_bm3"] = {"superseded_record_text": block}
    with pytest.raises(EosmatError, match=message):
        merge_audit(records, sidecar)


def test_bundled_blocks_read_through_the_record_audit():
    # The Rust reader exposes the same blocks (eosmat_audit.rs).
    found = 0
    for material in list_material_documents():
        sidecar = json.loads((AUDITS / f"{material}.audit.json").read_bytes())
        full = get_material_document(material, format_version=5)
        for record in full["eos_records"]:
            stored = (
                sidecar["records"]
                .get(record["identifier"], {})
                .get("superseded_record_text")
            )
            assert record.get("superseded_record_text") == stored
            found += stored is not None
    assert found
    record = get_eos_record("palladium_fedotenko_2020_bm2_2")
    assert record.name == "Fedotenko et al. (2020), constrained BM2"
    assert dict(record.audit)["superseded_record_text"]["label"] == (
        "Fedotenko et al. (2020), constrained BM2 (poor fit; diagnostic only)"
    )


def test_record_text_fields_round_trip():
    record = {
        "identifier": "r",
        "validity": {"notes": ["a"]},
        "thermal": {"model": "MGD"},
    }
    assert record_text(record, "validity_notes") == ["a"]
    assert record_text(record, "thermal_notes") is None
    set_record_text(record, "validity_notes", [])
    assert "validity" not in record
    set_record_text(record, "validity_notes", ["b"])
    set_record_text(record, "thermal_notes", "t")
    set_record_text(record, "notes", "n")
    set_record_text(record, "notes", None)
    assert record == {
        "identifier": "r",
        "validity": {"notes": ["b"]},
        "thermal": {"model": "MGD", "notes": "t"},
    }
    kept = {"identifier": "r", "validity": {"pressure_gpa": [0, 1], "notes": ["a"]}}
    set_record_text(kept, "validity_notes", None)
    assert kept["validity"] == {"pressure_gpa": [0, 1]}


def _catalog_copy(tmp_path, materials):
    """A minimal repository root with some materials and their receipt bindings."""
    receipt = json.loads(RECEIPT.read_bytes())
    receipt["materials"] = {m: receipt["materials"][m] for m in materials}
    receipt.pop("superseded", None)
    for directory in ("materials", "audits"):
        (tmp_path / "peritheos/data" / directory).mkdir(parents=True)
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
    return tmp_path


def _record(root, material, identifier):
    document = json.loads(
        (root / f"peritheos/data/materials/{material}.eosmat").read_bytes()
    )
    return next(r for r in document["eos_records"] if r["identifier"] == identifier)


def _entry(root, material, identifier):
    sidecar = json.loads(
        (root / f"peritheos/data/audits/{material}.audit.json").read_bytes()
    )
    return sidecar["records"].get(identifier, {})


def _change(root, material, identifier, **after):
    record = _record(root, material, identifier)
    return {
        material: {
            identifier: {
                "before": {name: record_text(record, name) for name in after},
                "after": after,
            }
        }
    }


AKIMOTOITE = "akimotoite_karki_2002_static_bm4_1"
LARNITE = "larnite_sagatova_2021_gga_300k_vinet"


def _without_text(record):
    result = copy.deepcopy(record)
    for name in SUPERSEDED_TEXT_FIELDS:
        set_record_text(result, name, None)
    return result


def test_apply_round_trip_on_two_materials(tmp_path):
    root = _catalog_copy(tmp_path, ["akimotoite", "larnite"])
    before = {
        m: (
            (root / f"peritheos/data/materials/{m}.eosmat").read_bytes(),
            (root / f"peritheos/data/audits/{m}.audit.json").read_bytes(),
            pre_split_bytes(m, root=root),
        )
        for m in ("akimotoite", "larnite")
    }
    old_akimotoite = _record(root, "akimotoite", AKIMOTOITE)
    old_larnite = _record(root, "larnite", LARNITE)
    batches = {
        **_change(
            root,
            "akimotoite",
            AKIMOTOITE,
            notes="Static-lattice LDA BM4.",
            validity_notes=None,
        ),
        **_change(root, "larnite", LARNITE, validity_notes=["Range note."]),
    }
    labels = {
        "larnite": {
            LARNITE: (old_larnite["label"], "Sagatova et al. (2021), GGA 300 K Vinet")
        }
    }
    plan = tool.plan_changes(batches, labels, root=root)
    assert plan.refusals == []
    ledger_path = root / tool.LEDGER
    shutil.copy2(ROOT / tool.LEDGER, ledger_path)
    ledger_before = json.loads(ledger_path.read_bytes())
    tool.apply_plan(plan, "2026-11-16", root)
    # The refit ledger's copies of label and notes follow; nothing else changes.
    ledger = json.loads(ledger_path.read_bytes())
    for old, new in zip(ledger_before["records"], ledger["records"]):
        identifier = new["record_identifier"]
        if identifier in (AKIMOTOITE, LARNITE):
            material = "akimotoite" if identifier == AKIMOTOITE else "larnite"
            record = _record(root, material, identifier)
            assert (new["label"], new["source_notes"]) == (
                record["label"],
                record.get("notes"),
            )
            old = {**old, "label": new["label"], "source_notes": new["source_notes"]}
        assert new == old
    from scripts.validate_primary_eos_refits import render_markdown

    assert (root / tool.LEDGER_MARKDOWN).read_text() == render_markdown(ledger)
    new_akimotoite = _record(root, "akimotoite", AKIMOTOITE)
    assert new_akimotoite["notes"] == "Static-lattice LDA BM4."
    assert "notes" not in new_akimotoite["validity"]
    assert _entry(root, "akimotoite", AKIMOTOITE)["superseded_record_text"] == {
        "notes": old_akimotoite["notes"],
        "validity_notes": old_akimotoite["validity"]["notes"],
        "superseded_on": "2026-11-16",
    }
    assert _entry(root, "larnite", LARNITE)["superseded_record_text"] == {
        "label": old_larnite["label"],
        "validity_notes": old_larnite["validity"]["notes"],
        "superseded_on": "2026-11-16",
    }
    # Only the four text fields changed; the files stay canonical.
    for material, identifier, old in (
        ("akimotoite", AKIMOTOITE, old_akimotoite),
        ("larnite", LARNITE, old_larnite),
    ):
        assert _without_text(_record(root, material, identifier)) == _without_text(old)
        path = root / f"peritheos/data/materials/{material}.eosmat"
        audit = root / f"peritheos/data/audits/{material}.audit.json"
        assert format_eosmat(json.loads(path.read_bytes())) == path.read_text()
        assert format_audit(json.loads(audit.read_bytes())) == audit.read_text()
        # The pre-split file is still recovered byte for byte.
        assert pre_split_bytes(material, root=root) == before[material][2]
    assert verify_catalog(root) == 2
    binding = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())["materials"]
    assert binding["akimotoite"]["superseded_record_text"] == [0]

    def read(relative):
        return before[Path(relative).stem][0]

    assert tool.verify(root=root, read=read) == (2, [])
    # A text change without a stored original is reported.
    path = root / "peritheos/data/materials/akimotoite.eosmat"
    document = json.loads(path.read_bytes())
    document["eos_records"][1]["notes"] = "Edited by hand."
    path.write_text(format_eosmat(document))
    checked, failures = tool.verify(root=root, read=read)
    assert checked == 2 and failures and "unrecorded change of notes" in failures[0]


def test_second_rewrite_keeps_the_first_original_and_a_revert_drops_it(tmp_path):
    root = _catalog_copy(tmp_path, ["akimotoite"])
    original = _record(root, "akimotoite", AKIMOTOITE)
    for text, date in (
        ("First rewrite.", "2026-11-16"),
        ("Second rewrite.", "2026-11-20"),
    ):
        plan = tool.plan_changes(
            _change(root, "akimotoite", AKIMOTOITE, notes=text), {}, root=root
        )
        assert plan.refusals == []
        tool.apply_plan(plan, date, root)
    block = _entry(root, "akimotoite", AKIMOTOITE)["superseded_record_text"]
    assert block == {"notes": original["notes"], "superseded_on": "2026-11-16"}
    assert pre_split_bytes("akimotoite", root=root) == pre_split_bytes("akimotoite")
    plan = tool.plan_changes(
        _change(root, "akimotoite", AKIMOTOITE, notes=original["notes"]), {}, root=root
    )
    tool.apply_plan(plan, "2026-11-21", root)
    assert "superseded_record_text" not in _entry(root, "akimotoite", AKIMOTOITE)
    assert (root / "peritheos/data/materials/akimotoite.eosmat").read_bytes() == (
        MATERIALS / "akimotoite.eosmat"
    ).read_bytes()
    assert (root / "peritheos/data/audits/akimotoite.audit.json").read_bytes() == (
        AUDITS / "akimotoite.audit.json"
    ).read_bytes()
    receipt = json.loads((root / RECEIPT.relative_to(ROOT)).read_bytes())
    assert "superseded_record_text" not in receipt["materials"]["akimotoite"]


def _refusals(root, batches, labels=None, allow=()):
    plan = tool.plan_changes(batches, labels or {}, allow=allow, root=root)
    return {(r.record, r.reason) for r in plan.refusals}, plan


@pytest.mark.parametrize(
    ("after", "guard"),
    [
        ({"notes": "Static LDA BM4 to 123.45 GPa."}, "numbers"),
        ({"notes": "Static LDA BM4 with V0 = 58.868(2)."}, "numbers"),
        ({"notes": "x" * 401}, "notes_length"),
        ({"notes": "One. Two. Three. Four."}, "notes_length"),
        ({"validity_notes": ["a", "b", "c"]}, "validity_length"),
        ({"validity_notes": ["x" * 161]}, "validity_length"),
        ({"notes": "See docs/eos-labels.md."}, "forbidden_text"),
        ({"notes": "Moved to the .json sidecar."}, "forbidden_text"),
        ({"notes": "Checked in Audit 2026."}, "forbidden_text"),
        ({"notes": "LitCurate import."}, "forbidden_text"),
    ],
)
def test_each_text_guard_refuses_and_allow_overrides_it(tmp_path, after, guard):
    root = _catalog_copy(tmp_path, ["akimotoite"])
    batches = _change(root, "akimotoite", AKIMOTOITE, **after)
    refusals, plan = _refusals(root, batches)
    assert (AKIMOTOITE, guard) in refusals and not plan.changes
    refusals, plan = _refusals(root, batches, allow=[f"{AKIMOTOITE}:{guard}"])
    assert (AKIMOTOITE, guard) not in refusals
    assert any(guard in line for line in plan.allowed)


def test_number_guard_accepts_numbers_from_text_and_structured_fields(tmp_path):
    root = _catalog_copy(tmp_path, ["akimotoite"])
    record = _record(root, "akimotoite", AKIMOTOITE)
    volume = record["parameters"]["V0"]
    volume = volume["value"] if isinstance(volume, dict) else volume
    notes = f"Static LDA BM4 (2002); V0 = {volume} Å³; Z=6."
    refusals, plan = _refusals(
        root, _change(root, "akimotoite", AKIMOTOITE, notes=notes)
    )
    assert refusals == set() and plan.changes


def test_thermal_label_and_consistency_refusals(tmp_path):
    root = _catalog_copy(tmp_path, ["mgo", "akimotoite"])
    thermal = "mgo_speziale_2001_bm3_2"
    refusals, _ = _refusals(
        root, _change(root, "mgo", thermal, thermal_notes="x" * 251)
    )
    assert (thermal, "thermal_length") in refusals
    records = json.loads(
        (root / "peritheos/data/materials/akimotoite.eosmat").read_bytes()
    )["eos_records"]
    first, second = records[0], records[1]
    labels = {
        "akimotoite": {
            first["identifier"]: (first["label"], "L" * 81),
            second["identifier"]: (second["label"], first["label"]),
        }
    }
    refusals, _ = _refusals(root, {}, labels)
    assert (first["identifier"], "label_length") in refusals
    assert (second["identifier"], "label_unique") in refusals
    stale = {"akimotoite": {first["identifier"]: ("Not the label", "New label")}}
    refusals, _ = _refusals(root, {}, stale)
    assert refusals == {(first["identifier"], "label_mismatch")}
    batches = _change(root, "akimotoite", AKIMOTOITE, notes="Short.")
    batches["akimotoite"][AKIMOTOITE]["before"]["notes"] = "Different text."
    refusals, _ = _refusals(root, batches)
    assert refusals == {(AKIMOTOITE, "before_mismatch")}


def test_reverse_procedure_recovers_every_revised_pre_split_file():
    # The lossless split test re-splits each pre-split file; for revised
    # materials the recorded revision is replayed on the split output and
    # must give the committed files byte for byte.
    receipt = json.loads(RECEIPT.read_bytes())
    revised = [
        m for m, b in receipt["materials"].items() if b.get("superseded_record_text")
    ]
    assert revised
    for material in revised:
        current = merged_document(material)
        indices = receipt["materials"][material]["superseded_record_text"]
        for index, record in enumerate(current["eos_records"]):
            assert ("superseded_record_text" in record) == (index in indices)
        restored = json.loads(pre_split_bytes(material))
        for index in indices:
            block = current["eos_records"][index]["superseded_record_text"]
            for name in SUPERSEDED_TEXT_FIELDS:
                if name in block:
                    assert (
                        record_text(restored["eos_records"][index], name) == block[name]
                    )


def test_scientific_snapshot_reads_thermal_notes_at_the_stored_original(tmp_path):
    root = tmp_path
    for directory in ("materials", "audits"):
        (root / "peritheos/data" / directory).mkdir(parents=True)
    document = standalone(thermal={"model": "MGD", "parameters": {}, "notes": "Old."})
    path = root / "peritheos/data/materials/example.eosmat"
    path.write_text(json.dumps(document))
    sidecar = {
        "format": "peritheos.eosmat-audit",
        "format_version": 5,
        "material": "example",
        "records": {},
    }
    audit = root / "peritheos/data/audits/example.audit.json"
    audit.write_text(json.dumps(sidecar))
    before = scientific_sha256(root)
    document["eos_records"][0]["thermal"]["notes"] = "New."
    sidecar["records"]["example_bm3"] = {
        "superseded_record_text": {
            "thermal_notes": "Old.",
            "superseded_on": "2026-11-16",
        }
    }
    path.write_text(json.dumps(document))
    audit.write_text(json.dumps(sidecar))
    assert scientific_sha256(root) == before


def test_record_text_limits_match_the_baseline():
    new, resolved = compare(current_cases(), load_baseline())
    assert new == [], "Record text exceeds docs/record-text-standard.md"
    assert resolved == [], "Run scripts/check_record_text.py --prune"
    assert compare({"m/r": ["notes_length"]}, {}) == (["m/r: notes_length"], [])
    assert compare({}, {"m/r": ["notes_length"]}) == ([], ["m/r: notes_length"])


def test_sentence_count_ignores_abbreviations_and_decimals():
    assert (
        tool.sentence_count("Baty et al. (2024) fit V0 = 14.72(4) Å³. Second one.") == 2
    )
    assert tool.sentence_count("See Fig. 3 and Eq. 2, e.g. the fit") == 1
    assert tool.text_numbers("V0 = 58.868(2) Å³ at 300 K, Z=4, 2026-10-17") == [
        "58.868",
        "300",
        "2026",
        "10",
        "17",
    ]
