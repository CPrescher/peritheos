"""EOSMAT v5 record ``check`` headlines and per-material audit sidecars."""

import copy
import json
import shutil
from dataclasses import replace
from pathlib import Path

import numpy as np
import pytest

import peritheos._eosmat_audit as audit_module
import peritheos_build
from peritheos import (
    Material,
    convert_eosmat,
    get_eos_record,
    get_material,
    get_material_document,
    load_eosmat,
    save_eosmat,
)
from peritheos._eosmat_audit import (
    AUDIT_RECORD_KEYS,
    GENERIC_QUALIFICATIONS,
    GENERIC_VALIDATION_NOTES,
    SUPERSEDED_TEXT_FIELDS,
    check_block,
    merge_audit,
    remove_generic_statements,
    reproduction_headline,
    split_document,
    usage_warning,
)
from peritheos._eosmat_format import format_audit, format_eosmat
from peritheos.eosmat import list_material_documents
from peritheos.errors import EosmatError
from peritheos.materials import _material_from_eosmat
from scripts.apply_record_text import revise_record
from scripts.split_eosmat_audits import (
    RECEIPT,
    merged_document,
    pre_split_bytes,
    record_text,
    restore_record_text,
    set_record_text,
    split_added_materials,
    split_material,
    verify_catalog,
)

ROOT = Path(__file__).resolve().parents[1]
MATERIALS = ROOT / "peritheos/data/materials"
AUDITS = ROOT / "peritheos/data/audits"


def standalone(**record_extra):
    record = {
        "identifier": "example_bm3",
        "model": "BM3",
        "temperature_ref": 300,
        "parameters": {"V0": 74.7, "K0": 160, "K0_prime": 4},
        "reference": "Fictional example; no experimental validation.",
        **record_extra,
    }
    return {
        "format_version": 5,
        "identifier": "example",
        "formula": "MgO",
        "eos_records": [record],
    }


def evidence(status="parity", **validation):
    return {
        "scientific_validation": {
            "status": "primary_source_validated",
            "audit_date": "2026-01-02",
            "primary_source_check": {"finding": "Table 2 values match."},
            **validation,
        },
        "reproduction": {
            "assessment": {
                "status": status,
                "scope": "parameter_agreement",
                "basis": "Independent refit of Table 1.",
            }
        },
    }


@pytest.mark.parametrize(
    ("status", "headline"),
    [
        ("parity", "reproduced"),
        ("similar", "close"),
        ("not_reproduced", "disagrees"),
        ("not_refittable", "untested"),
        ("source_reconstruction", "untested"),
        ("not_assessed", "untested"),
    ],
)
def test_reproduction_headline_is_a_fixed_mapping(status, headline):
    record = standalone(evidence=evidence(status))["eos_records"][0]
    assert reproduction_headline(record) == headline
    document = standalone(evidence=evidence(status))
    loaded = Material.from_eosmat(document).eos_records[0]
    assert loaded.reproduction == headline


def test_a_record_without_assessment_is_untested_never_positive():
    for record in (
        standalone()["eos_records"][0],
        standalone(evidence={"scientific_validation": {"status": "deferred"}})[
            "eos_records"
        ][0],
    ):
        assert reproduction_headline(record) == "untested"
    loaded = Material.from_eosmat(
        standalone(), require_primary_validation=False
    ).eos_records[0]
    assert loaded.reproduction == "untested"
    assert loaded.warning is None
    assert dict(loaded.audit) == {}


def test_check_block_copies_stored_fields_verbatim():
    record = standalone(evidence=evidence("similar"))["eos_records"][0]
    assert check_block(record, audit="a.audit.json#example_bm3") == {
        "reproduction": "close",
        "compared": "published_fit",
        "warning": None,
        "source_check": "verified",
        "source_checked_on": "2026-01-02",
        "source_note": "Table 2 values match.",
        "audit": "a.audit.json#example_bm3",
    }


def test_warning_uses_only_coded_restrictions():
    coded = standalone(
        evidence=evidence(usage_recommendation="not_recommended_for_quantitative_use")
    )["eos_records"][0]
    assert usage_warning(coded) == "Not recommended for quantitative use."
    prose = standalone(
        evidence=evidence(usage_recommendation="DO NOT USE outside 2-10 GPa.")
    )["eos_records"][0]
    assert usage_warning(prose) is None
    named_warning = standalone(evidence=evidence(warning="Possible melting."))[
        "eos_records"
    ][0]
    assert usage_warning(named_warning) is None
    flagged = standalone(
        catalog_access="explicit_selection",
        evidence={**evidence(), "reproduction": {"do_not_use": True}},
    )["eos_records"][0]
    assert usage_warning(flagged) == (
        "Do not use. Never chosen automatically; select it explicitly."
    )
    held = {"recommendation_hold": {"status": "on_hold", "record": "example_bm3"}}
    assert usage_warning(coded, held) == (
        "Not recommended for quantitative use. "
        "Default selection on hold pending review."
    )


def test_catalog_checks_equal_their_derivation_from_the_sidecar():
    counts = dict.fromkeys(("reproduced", "close", "disagrees", "untested"), 0)
    for material in list_material_documents():
        stored = json.loads((MATERIALS / f"{material}.eosmat").read_bytes())
        keys = list(stored)
        if "datasets" in stored:
            assert keys.index("datasets") > keys.index("eos_records")
        full = get_material_document(material, format_version=5)
        for record, complete in zip(stored["eos_records"], full["eos_records"]):
            assert not set(record) & set(AUDIT_RECORD_KEYS), record["identifier"]
            expected = check_block(
                {k: v for k, v in complete.items() if k != "check"},
                full,
                audit=f"../audits/{material}.audit.json#{record['identifier']}",
            )
            assert record["check"] == expected
            counts[record["check"]["reproduction"]] += 1
    assert counts == {
        "reproduced": 358,
        "close": 72,
        "disagrees": 61,
        "untested": 196,
    }


def _revised_indices(material):
    """Records whose text was revised after the split (receipt binding)."""
    binding = json.loads(RECEIPT.read_bytes())["materials"].get(material, {})
    return binding.get("superseded_record_text", [])


def _without_text_revision(document, material):
    """A merged document with the post-split text revision undone."""
    result = copy.deepcopy(document)
    for index in _revised_indices(material):
        restore_record_text(result["eos_records"][index])
    return result


def _replay_text_revision(material, record_bytes, audit_bytes):
    """Apply the bound text revision to freshly split files.

    The revised text is taken from the committed record file; the stored
    originals are recomputed from the split output, so the result equals the
    committed files only if every stored original is the pre-split text.
    """
    indices = _revised_indices(material)
    if not indices:
        return record_bytes, audit_bytes
    current = merged_document(material)
    document, sidecar = json.loads(record_bytes), json.loads(audit_bytes)
    for index in indices:
        revised = current["eos_records"][index]
        block = revised["superseded_record_text"]
        record = document["eos_records"][index]
        new = {f: record_text(revised, f) for f in SUPERSEDED_TEXT_FIELDS if f in block}
        entry = sidecar["records"].setdefault(record["identifier"], {})
        revise_record(record, entry, new, block["superseded_on"])
    order = [r["identifier"] for r in document["eos_records"]]
    entries = sidecar["records"]
    sidecar["records"] = {i: entries[i] for i in order if i in entries}
    return format_eosmat(document).encode(), format_audit(sidecar).encode()


def test_every_material_has_exactly_one_sidecar_and_the_split_is_lossless():
    # Each pre-split file is recovered byte for byte from record file, sidecar
    # and the receipt; re-splitting it and replaying the bound text revision
    # (superseded_record_text) reproduces the committed files exactly.
    assert verify_catalog() == len(list_material_documents())
    revised = 0
    for material in list_material_documents():
        record_bytes, audit_bytes, _ = split_material(pre_split_bytes(material))
        record_bytes, audit_bytes = _replay_text_revision(
            material, record_bytes, audit_bytes
        )
        revised += bool(_revised_indices(material))
        assert record_bytes == (MATERIALS / f"{material}.eosmat").read_bytes()
        assert audit_bytes == (AUDITS / f"{material}.audit.json").read_bytes()
    assert revised


def test_only_documented_statements_are_removed_and_numbers_are_bit_identical():
    for material in list_material_documents():
        before = json.loads(pre_split_bytes(material))
        merged = _without_text_revision(merged_document(material), material)
        defaults = {}
        for old, new in zip(before["eos_records"], merged["eos_records"]):
            new = {k: v for k, v in new.items() if k != "check"}
            flags = {k: old.pop(k) for k in ("default", "default_for") if k in old}
            if flags.get("default") is True or flags.get("default_for"):
                category = (
                    "hugoniot"
                    if old["model"] == "LinearUsUpHugoniot"
                    else "equilibrium"
                )
                defaults[category] = old["identifier"]
            assert remove_generic_statements(old) == new
            for key in ("parameters", "thermal", "parameter_covariance"):
                assert json.dumps(old.get(key)) == json.dumps(new.get(key))
        # The per-record flags became one map naming exactly the same records.
        assert merged.get("default_records", {}) == defaults
        assert {
            k: v
            for k, v in merged.items()
            if k not in ("eos_records", "default_records")
        } == {k: v for k, v in before.items() if k != "eos_records"}
        assert json.dumps(before.get("datasets")) == json.dumps(merged.get("datasets"))
    generic = {*GENERIC_QUALIFICATIONS, *GENERIC_VALIDATION_NOTES}
    for path in AUDITS.glob("*.audit.json"):
        for identifier, entry in json.loads(path.read_bytes())["records"].items():
            evidence = entry["evidence"]
            assessment = evidence["reproduction"].get("assessment", {})
            # Longer record-specific notes and nested historical assessments
            # that merely contain a generic sentence are kept verbatim.
            assert not generic & set(assessment.get("qualifications", []))
            assert evidence["scientific_validation"].get("note") not in generic
            assert not {"report_url", "report_path", "report"} & set(assessment)


def test_sidecar_loads_lazily_and_matches_an_eager_load(monkeypatch):
    calls = []
    original = audit_module.load_bundled_audit

    def counting(material):
        calls.append(material)
        return original(material)

    monkeypatch.setattr(audit_module, "load_bundled_audit", counting)
    records_only = get_material_document(
        "mgo", format_version=None, include_audit=False
    )
    assert all("evidence" not in r for r in records_only["eos_records"])
    lazy = _material_from_eosmat(
        records_only, require_primary_validation=False, bundled_audit=True
    )
    record = lazy.eos_records[0]
    assert calls == []
    assert (
        record.reproduction == records_only["eos_records"][0]["check"]["reproduction"]
    )
    assert record.determination_method == "experimental"
    assert record.pressure(70.0) > 0
    assert calls == []
    assert record.scientific_validation_note
    assert calls == ["mgo"]
    eager = Material.from_eosmat(
        get_material_document("mgo", format_version=5),
        require_primary_validation=False,
    )
    assert calls == ["mgo", "mgo"]
    for a, b in zip(lazy.eos_records, eager.eos_records):
        for name in (
            "identifier",
            "parameter_provenance",
            "pressure_calibration",
            "scientific_validation_status",
            "scientific_validation_note",
            "parameter_errors",
            "parameter_covariance",
            "validity",
            "reference",
            "notes",
            "is_default",
        ):
            assert getattr(a, name) == getattr(b, name), name
        assert a.eos.parameter_values() == b.eos.parameter_values()
        assert dict(a.eosmat_metadata) == dict(b.eosmat_metadata)
        assert dict(a.audit) == dict(b.audit)
        assert (a.reproduction, a.warning) == (b.reproduction, b.warning)
    assert calls == ["mgo", "mgo"]
    assert lazy.to_eosmat() == eager.to_eosmat()


def test_catalog_accessors_return_complete_evidence():
    record = get_eos_record("mgo_sakai_2025_rydberg_stacey_1")
    full = next(
        r
        for r in get_material_document("mgo", format_version=5)["eos_records"]
        if r["identifier"] == record.identifier
    )
    assert record.reproduction == "untested"
    assert record.warning is None
    assert dict(record.audit)["evidence"] == full["evidence"]
    assert dict(record.pressure_calibration) == full["pressure_calibration"]
    validation = record.eosmat_metadata["scientific_validation"]
    assert validation["primary_data_check"]["status"] == "bundled"
    exported = get_material("mgo").to_eosmat()
    by_id = {r["identifier"]: r for r in exported["eos_records"]}
    assert (
        by_id[record.identifier]["evidence"]["source_lineage"]
        == full["evidence"]["source_lineage"]
    )
    assert by_id[record.identifier]["check"] == full["check"]


def _default_flags_as_selection(document):
    """Replace per-record default flags by the selection they express.

    A ``default_records`` map is exported as ``default_for`` on the named
    record; the pre-split files spelled the same selection as ``default:
    true``, ``default_for`` or both, and marked others ``default: false``.
    """
    result = copy.deepcopy(document)
    for record in result["eos_records"]:
        flags = {
            key: record.pop(key) for key in ("default", "default_for") if key in record
        }
        record["selected_default"] = bool(
            flags.get("default") is True or flags.get("default_for")
        )
    return result


def test_compatibility_exports_contain_the_full_evidence():
    for material in list_material_documents():
        full = _without_text_revision(
            get_material_document(material, format_version=5), material
        )
        before = json.loads(pre_split_bytes(material))
        before["eos_records"] = [
            remove_generic_statements(r) for r in before["eos_records"]
        ]
        for version in (3, 4):
            exported = convert_eosmat(full, version=version)
            assert "default_records" not in exported
            # Identical apart from the spelling of the same default selection.
            assert _default_flags_as_selection(exported) == _default_flags_as_selection(
                convert_eosmat(before, version=version)
            )
            assert all("check" not in r for r in exported["eos_records"])
            if "default_records" in full:
                for category, identifier in full["default_records"].items():
                    record = next(
                        r
                        for r in exported["eos_records"]
                        if r["identifier"] == identifier
                    )
                    assert record["default_for"] == category
        dioptas = _without_text_revision(get_material_document(material), material)
        assert _default_flags_as_selection(dioptas) == _default_flags_as_selection(
            convert_eosmat(before, version=3, legacy_uncertainty_units=False)
        )
        for record in dioptas["eos_records"]:
            assert "pressure_calibration" in record
            assert "scientific_validation" in record


def test_material_exports_of_catalog_materials_are_complete():
    material = get_material("palladium")
    full = get_material_document("palladium", format_version=5)
    for version in (3, 4):
        exported = material.to_eosmat(format_version=version)
        evidence_keys = (
            {"evidence"} if version == 4 else {"scientific_validation", "reproduction"}
        )
        for record in exported["eos_records"]:
            assert evidence_keys <= set(record)
            assert "check" not in record
    native = {r["identifier"]: r for r in material.to_eosmat()["eos_records"]}
    for record in full["eos_records"]:
        if record["identifier"] in native:
            assert native[record["identifier"]]["evidence"] == record["evidence"]


def test_bundled_records_evaluate_identical_pressures_before_and_after():
    for material in list_material_documents():
        try:
            current = get_material(material)
        except Exception:
            continue  # Document-only materials have no executable records.
        if current.identifier != material:
            continue
        before = _material_from_eosmat(
            json.loads(pre_split_bytes(material)),
            require_primary_validation=False,
            record_identifiers=[r.identifier for r in current.eos_records],
        )
        for old, new in zip(before.eos_records, current.eos_records):
            assert old.identifier == new.identifier
            if old.is_hugoniot:
                volumes = old.reference_volume * np.array([0.8, 0.9])
                np.testing.assert_array_equal(
                    old.pressure(volumes, check_domain=False),
                    new.pressure(volumes, check_domain=False),
                )
                continue
            volumes = old.reference_volume * np.array([0.85, 0.95, 1.0])
            temperature = old.reference_temperature
            np.testing.assert_array_equal(
                old.pressure(volumes, temperature, check_validity=False),
                new.pressure(volumes, temperature, check_validity=False),
            )
            assert old.parameter_errors == new.parameter_errors
            assert old.parameter_error_confidence == new.parameter_error_confidence
            assert old.parameter_covariance == new.parameter_covariance
            assert old.covariance_parameters == new.covariance_parameters
            if "superseded_record_text" in new.audit:
                # Only the validity notes text may have been revised.
                stored = new.audit["superseded_record_text"]
                notes = stored.get("validity_notes", list(new.validity.notes)) or ()
                assert old.validity == replace(new.validity, notes=tuple(notes))
            else:
                assert old.validity == new.validity
            assert old.is_default == new.is_default
            assert old.equation_kind == new.equation_kind
            assert old.reproduction == new.reproduction
            for key in ("record_kind", "fixed_parameters"):
                assert old.eosmat_metadata.get(key) == new.eosmat_metadata.get(key)
            # Thermal notes are text and may have been revised; the stored
            # original restores them, and every other thermal field is equal.
            thermal = {"thermal": copy.deepcopy(new.eosmat_metadata.get("thermal"))}
            stored = dict(new.audit).get("superseded_record_text", {})
            if "thermal_notes" in stored:
                set_record_text(thermal, "thermal_notes", stored["thermal_notes"])
            assert old.eosmat_metadata.get("thermal") == thermal["thermal"]
            assert old.scientific_validation_status == new.scientific_validation_status
            assert old.parameter_provenance == new.parameter_provenance
            assert old.pressure_calibration == new.pressure_calibration
        for selector in (
            "default_record",
            "default_equilibrium_record",
            "default_hugoniot_record",
        ):
            assert getattr(getattr(before, selector)(), "identifier", None) == getattr(
                getattr(current, selector)(), "identifier", None
            ), (material, selector)


def test_standalone_inline_evidence_needs_no_sidecar(tmp_path):
    document = standalone(evidence=evidence("parity"))
    path = tmp_path / "example.eosmat"
    save_eosmat(path, document)
    assert load_eosmat(path) == document
    assert load_eosmat(path, include_audit=False) == document
    record = Material.from_eosmat(load_eosmat(path)).eos_records[0]
    assert record.reproduction == "reproduced"
    assert (
        dict(record.audit)["evidence"]["reproduction"]["assessment"]["status"]
        == "parity"
    )
    assert record.scientific_validation_status == "primary_source_validated"
    # An explicit check must agree with the evidence it summarizes.
    contradicted = standalone(
        evidence=evidence("parity"), check={"reproduction": "disagrees"}
    )
    with pytest.raises(EosmatError, match="contradicts"):
        Material.from_eosmat(contradicted)
    with pytest.raises(EosmatError, match="check.reproduction"):
        Material.from_eosmat(standalone(check={"reproduction": "great"}))
    with_warning = standalone(check={"warning": "Author-supplied caution."})
    loaded = Material.from_eosmat(with_warning, require_primary_validation=False)
    assert loaded.eos_records[0].warning == "Author-supplied caution."
    assert loaded.eos_records[0].reproduction == "untested"


def test_user_split_pair_merges_on_request_and_missing_sidecar_is_an_error(tmp_path):
    full = get_material_document("larnite", format_version=5)
    (tmp_path / "materials").mkdir()
    (tmp_path / "audits").mkdir()
    path = tmp_path / "materials/larnite.eosmat"
    shutil.copy2(MATERIALS / "larnite.eosmat", path)
    shutil.copy2(AUDITS / "larnite.audit.json", tmp_path / "audits")
    assert load_eosmat(path, include_audit=False) == json.loads(path.read_bytes())
    assert load_eosmat(path) == full
    record_only = Material.from_eosmat(
        load_eosmat(path, include_audit=False)
    ).eos_records[0]
    assert record_only.reproduction == full["eos_records"][0]["check"]["reproduction"]
    with pytest.raises(EosmatError, match="was not loaded"):
        record_only.audit
    (tmp_path / "audits/larnite.audit.json").unlink()
    with pytest.raises(EosmatError, match="is missing") as error:
        load_eosmat(path)
    assert error.value.code == "eosmat.audit_missing"
    exported = convert_eosmat(load_eosmat(path, include_audit=False), version=4)
    assert all("check" in r for r in exported["eos_records"])


def test_missing_bundled_sidecar_is_a_clear_error(monkeypatch):
    audit_module._bundled_audit_text.cache_clear()
    monkeypatch.setattr(audit_module, "audit_file_name", lambda m: f"{m}.absent.json")
    try:
        with pytest.raises(EosmatError, match="audit sidecar for bundled") as error:
            get_material_document("larnite")
        assert error.value.code == "eosmat.audit_missing"
        lazy = _material_from_eosmat(
            get_material_document("larnite", format_version=None, include_audit=False),
            require_primary_validation=False,
            bundled_audit=True,
        )
        record = lazy.eos_records[0]
        assert np.isfinite(record.pressure(0.95 * record.reference_volume, 300.0))
        assert record.reproduction in {"reproduced", "close", "disagrees", "untested"}
        with pytest.raises(EosmatError, match="Reinstall"):
            record.scientific_validation_note
    finally:
        audit_module._bundled_audit_text.cache_clear()


def test_split_and_merge_round_trip_and_reject_conflicts():
    document = standalone(
        evidence=evidence("parity"),
        pressure_calibration={
            "status": "unresolved",
            "methods": [],
            "recalculation": {"status": "not_possible"},
            "audit_date": "2026-01-02",
        },
    )
    derived = audit_module.derived_report_links("example", "example_bm3")
    document["eos_records"][0]["evidence"]["reproduction"]["assessment"].update(derived)
    records, sidecar = split_document(document)
    assert set(records["eos_records"][0]) & set(AUDIT_RECORD_KEYS) == set()
    stored = sidecar["records"]["example_bm3"]["evidence"]["reproduction"]
    assert not set(derived) & set(stored["assessment"])
    merged = merge_audit(records, sidecar)
    without_check = copy.deepcopy(merged)
    del without_check["eos_records"][0]["check"]
    assert without_check == document
    clash = copy.deepcopy(records)
    clash["eos_records"][0]["evidence"] = {}
    with pytest.raises(EosmatError, match="both inline"):
        merge_audit(clash, sidecar)
    stray = copy.deepcopy(sidecar)
    stray["records"]["unknown_record"] = {}
    with pytest.raises(EosmatError, match="unknown records"):
        merge_audit(records, stray)


def test_build_staging_uses_the_same_audit_keys():
    assert peritheos_build.AUDIT_RECORD_KEYS == AUDIT_RECORD_KEYS
    for material in ("mgo", "larnite", "palladium"):
        assert peritheos_build.source_document(
            MATERIALS / f"{material}.eosmat"
        ) == get_material_document(material, format_version=None)


def _split_catalog_copy(tmp_path, materials):
    """A minimal repository root holding some split materials and the receipt."""
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
    return target


def test_added_inline_file_is_split_deterministically_and_recovers(tmp_path):
    receipt_path = _split_catalog_copy(tmp_path, ["larnite", "pyrope"])
    larnite = (tmp_path / "peritheos/data/materials/larnite.eosmat").read_bytes()
    # A curator writes a complete inline-evidence file for a new material and
    # replaces an existing split material with its complete inline document.
    inline = pre_split_bytes("grossular")
    (tmp_path / "peritheos/data/materials/grossular.eosmat").write_bytes(inline)
    pyrope = pre_split_bytes("pyrope")
    (tmp_path / "peritheos/data/materials/pyrope.eosmat").write_bytes(pyrope)
    old_pyrope = json.loads(receipt_path.read_bytes())["materials"]["pyrope"]
    with pytest.raises(ValueError, match="already split"):
        split_added_materials(["larnite"], root=tmp_path, check_snapshot=False)
    assert split_added_materials(root=tmp_path, check_snapshot=False) == [
        "grossular",
        "pyrope",
    ]
    for material, original in (("grossular", inline), ("pyrope", pyrope)):
        assert pre_split_bytes(material, root=tmp_path) == original
        record_bytes, audit_bytes, _ = split_material(original)
        assert (
            tmp_path / f"peritheos/data/materials/{material}.eosmat"
        ).read_bytes() == record_bytes
        assert (
            tmp_path / f"peritheos/data/audits/{material}.audit.json"
        ).read_bytes() == audit_bytes
        # The same input gives the same files as the bundled catalog, once the
        # catalog's bound text revision is replayed.
        record_bytes, audit_bytes = _replay_text_revision(
            material, record_bytes, audit_bytes
        )
        assert record_bytes == (MATERIALS / f"{material}.eosmat").read_bytes()
        assert audit_bytes == (AUDITS / f"{material}.audit.json").read_bytes()
    receipt = json.loads(receipt_path.read_bytes())
    assert list(receipt["materials"]) == ["grossular", "larnite", "pyrope"]
    assert receipt["superseded"] == {"pyrope": [old_pyrope]}
    assert (
        tmp_path / "peritheos/data/materials/larnite.eosmat"
    ).read_bytes() == larnite
    # Nothing is left to split; a second run changes nothing.
    snapshot = receipt_path.read_bytes()
    assert split_added_materials(root=tmp_path, check_snapshot=False) == []
    assert receipt_path.read_bytes() == snapshot
    assert verify_catalog(tmp_path) == 3


def test_complete_document_of_a_split_material_resplits_to_the_same_files(tmp_path):
    receipt_path = _split_catalog_copy(tmp_path, ["larnite"])
    path = tmp_path / "peritheos/data/materials/larnite.eosmat"
    save_eosmat(path, get_material_document("larnite", format_version=5))
    assert split_added_materials(["larnite"], root=tmp_path, check_snapshot=False)
    assert path.read_bytes() == (MATERIALS / "larnite.eosmat").read_bytes()
    assert (tmp_path / "peritheos/data/audits/larnite.audit.json").read_bytes() == (
        AUDITS / "larnite.audit.json"
    ).read_bytes()
    assert "larnite" in json.loads(receipt_path.read_bytes())["superseded"]


def test_added_split_rejects_mixed_and_noncanonical_files(tmp_path):
    _split_catalog_copy(tmp_path, ["larnite"])
    path = tmp_path / "peritheos/data/materials/larnite.eosmat"
    mixed = json.loads(path.read_bytes())
    complete = get_material_document("larnite", format_version=5)
    mixed["eos_records"].append(
        {**complete["eos_records"][0], "identifier": "larnite_inline_copy"}
    )
    mixed["eos_records"][-1].pop("check")
    path.write_text(json.dumps(mixed))
    with pytest.raises(ValueError, match="inline audit blocks while"):
        split_added_materials(root=tmp_path, check_snapshot=False)
    path.write_text(json.dumps(complete))
    with pytest.raises(ValueError, match="canonical"):
        split_added_materials(root=tmp_path, check_snapshot=False)


def test_curator_inline_file_with_record_flags_splits_to_the_bundled_layout(tmp_path):
    receipt_path = _split_catalog_copy(tmp_path, ["mgo"])
    complete = get_material_document("mgo", format_version=5)
    # A curator writes inline evidence with per-record default flags and no
    # check blocks, in the current layout.
    defaults = set(complete.pop("default_records").values())
    assert len(defaults) == 2  # equilibrium and Hugoniot
    for record in complete["eos_records"]:
        record.pop("check")
        if record["identifier"] in defaults:
            record["default"] = True
    path = tmp_path / "peritheos/data/materials/mgo.eosmat"
    save_eosmat(path, complete)
    inline = path.read_bytes()
    assert split_added_materials(["mgo"], root=tmp_path, check_snapshot=False)
    assert path.read_bytes() == (MATERIALS / "mgo.eosmat").read_bytes()
    assert (tmp_path / "peritheos/data/audits/mgo.audit.json").read_bytes() == (
        AUDITS / "mgo.audit.json"
    ).read_bytes()
    binding = json.loads(receipt_path.read_bytes())["materials"]["mgo"]
    assert "layout" not in binding
    assert binding["added_default_records"] is True
    assert pre_split_bytes("mgo", root=tmp_path) == inline
