"""Short reprs, plain-text summaries and forgiving material lookup."""

import copy
import dataclasses
import shutil
from pathlib import Path

import numpy as np
import pytest

import peritheos._eosmat_audit as audit_module
from peritheos import (
    Material,
    get_eos_record,
    get_material,
    get_material_document,
    list_materials,
    load_eosmat,
)
from peritheos.catalog import _catalog_index, _compatibility_materials
from peritheos.errors import MaterialLookupError
from peritheos.materials import _material_from_eosmat

ROOT = Path(__file__).resolve().parents[1]


def lazy_material(identifier):
    """A bundled material whose audit sidecar has not been read."""
    return _material_from_eosmat(
        get_material_document(identifier, format_version=None, include_audit=False),
        require_primary_validation=False,
        bundled_audit=True,
    )


@pytest.fixture
def no_sidecar_reads(monkeypatch):
    calls = []

    def refuse(material):
        calls.append(material)
        raise AssertionError(f"audit sidecar of {material!r} was read")

    monkeypatch.setattr(audit_module, "load_bundled_audit", refuse)
    return calls


# --- repr -----------------------------------------------------------------


def test_material_and_record_reprs_are_one_short_line():
    material = get_material("mgo")
    assert repr(material) == (
        "Material('mgo', formula='MgO', phase='B1 (periclase)', 37 EOS records)"
    )
    assert str(material) == repr(material)
    record = material.get_eos_record("mgo_speziale_2001_bm3_2")
    assert repr(record) == (
        "EOSRecord('mgo_speziale_2001_bm3_2', "
        "model='BM3 + VariableQMieGruneisenDebye', reproduction='reproduced', "
        "default)"
    )
    warned = material.get_eos_record("mgo_shen_2026_vinet_3")
    assert warned.warning
    assert "warning='Not recommended" in repr(warned)
    hugoniot = material.get_eos_record("mgo_b1_duffy_ahrens_1995_hugoniot_5")
    assert repr(hugoniot).startswith("HugoniotRecord('mgo_b1_duffy_ahrens_1995")


def test_reprs_of_every_catalog_object_are_bounded():
    materials = (
        *_catalog_index().materials.values(),
        *_compatibility_materials().values(),
    )
    for material in materials:
        text = repr(material)
        assert "\n" not in text and len(text) < 250, material.identifier
        for record in material.eos_records:
            text = repr(record)
            assert "\n" not in text and len(text) < 250, record.identifier
    # Nested containers no longer expand every record's metadata.
    assert len(repr(list_materials())) < 60_000


def test_repr_and_summary_never_read_the_audit_sidecar(no_sidecar_reads):
    material = lazy_material("mgo")
    repr(material)
    material.summary()
    for record in material.eos_records:
        repr(record)
        record.summary()
    assert no_sidecar_reads == []


def test_repr_survives_a_missing_bundled_sidecar(monkeypatch):
    audit_module._bundled_audit_text.cache_clear()
    monkeypatch.setattr(audit_module, "audit_file_name", lambda m: f"{m}.absent.json")
    try:
        material = lazy_material("larnite")
        record = material.eos_records[0]
        assert repr(record).startswith(f"EOSRecord({record.identifier!r}")
        assert record.identifier in material.summary()
        assert record.identifier in record.summary()
    finally:
        audit_module._bundled_audit_text.cache_clear()


def test_repr_survives_an_unloaded_user_sidecar(tmp_path):
    (tmp_path / "materials").mkdir()
    path = tmp_path / "materials/larnite.eosmat"
    shutil.copy2(ROOT / "peritheos/data/materials/larnite.eosmat", path)
    material = Material.from_eosmat(load_eosmat(path, include_audit=False))
    record = material.eos_records[0]
    assert record.identifier in repr(record)
    assert record.identifier in record.summary()


def test_repr_never_raises_for_a_broken_record():
    record = copy.copy(get_eos_record("mgo_speziale_2001_bm3_2"))
    record.__dict__["_check_headline"] = None
    record.__dict__.pop("_audit_source", None)
    record.__dict__["eosmat_metadata"] = None  # unreadable on purpose
    assert repr(record) == "EOSRecord('mgo_speziale_2001_bm3_2')"


def test_equality_and_replace_are_unchanged():
    record = get_eos_record("gold_dewaele_2004_vinet_5")
    clone = dataclasses.replace(record)
    assert clone == record
    assert repr(clone) == repr(record)
    material = Material.from_eosmat(get_material_document("gold", format_version=5))
    assert dataclasses.replace(material) == material
    fields = {item.name for item in dataclasses.fields(Material)}
    assert {"identifier", "formula", "phase", "eos_records"} <= fields


# --- summaries ------------------------------------------------------------


def test_material_summary_lists_every_record_with_status():
    material = get_material("mgo")
    text = material.summary()
    lines = text.splitlines()
    assert lines[0] == "Periclase B1 (mgo)"
    assert "Formula: MgO" in text and "B1 (periclase)" in text
    assert "4 formula units per conventional unit cell" in text
    assert "angstrom^3/conventional_unit_cell" in text
    default = material.default_record()
    assert f"  Default: {default.name} ({default.identifier})" in lines
    for record in material.eos_records:
        row = next(line for line in lines if f" {record.identifier} " in line)
        assert record.reproduction in row
        assert str(record.reference.year) in row
        assert ("*" in row.split()[0]) is record.is_default
        assert ("!" in row.split()[0]) is bool(record.warning)
    warned = material.get_eos_record("mgo_shen_2026_vinet_3")
    warning_line = lines.index(f"    ! {warned.warning}")
    assert warned.identifier in lines[warning_line + 1]
    # Each distinct warning is stated once.
    assert text.count(warned.warning) == 1
    assert "Speziale et al. 2001" in text
    assert "0-203" in text and "300-3663" in text
    # Unstated ranges are shown as "-", not as sentinel numbers.
    assert "e-324" not in text and "inf" not in text


def test_every_catalog_summary_renders_without_truncating_identifiers():
    for material in list_materials():
        text = material.summary()
        for record in material.eos_records:
            assert f" {record.identifier} " in text
            assert record.identifier in record.summary()


def test_record_summary_shows_coefficients_reference_and_finding():
    record = get_eos_record("mgo_speziale_2001_bm3_2")
    text = record.summary()
    lines = [line.rstrip() for line in text.splitlines()]
    assert lines[0] == "mgo_speziale_2001_bm3_2"
    assert "  Default: the designated equilibrium record" in lines
    assert any(
        line.split() == ["V0", "=", "74.71", "+/-", "0.01", record.volume_unit]
        for line in lines
    )
    assert any(line.split() == ["K0", "=", "160.2", "fixed", "GPa"] for line in lines)
    assert any(line.split() == ["q1", "=", "11.8", "+/-", "0.2"] for line in lines)
    # A fixed coefficient's published error is not shown as a fit error.
    assert "gamma0 = 1.524 fixed (published +/- 0.025)" in text
    assert "confidence level not stated" in text
    assert (
        "  Check: reproduced (refit vs published fit). A refit of the source data "
        "gave the published coefficients."
    ) in lines
    assert "Warning" not in text
    assert (
        "    Validity range (enforced by check_validity=True): P 0-203, T 300-3663"
    ) in lines
    assert "    Experimental data range: P 0-52" in lines
    assert "doi:10.1029/2000JB900318" in text
    assert "  Source check: verified (2026-10-09). Published final BM3" in text
    assert "No documented unit for: q0, q1" in text
    assert "parameter_covariance" not in text


def test_record_summary_converts_molar_volume_errors_to_the_public_unit():
    # Python-constructed compatibility records carry no record-file metadata.
    record = get_material("mgo_b1").get_eos_record("mgo_b1_tange_2009_vinet")
    text = record.summary()
    v0 = next(line for line in text.splitlines() if line.strip().startswith("V0"))
    assert f"{record.reference_volume:.12g}" in v0
    assert "angstrom^3" in v0


def test_record_summary_states_the_confidence_level():
    record = get_eos_record("gold_dewaele_2004_vinet_5")
    text = record.summary()
    assert "95% confidence" in text
    assert any(
        line.split()[:5] == ["K0'", "=", "6", "+/-", "0.02"]
        for line in text.splitlines()
    )


def test_record_summary_shows_the_warning_text():
    record = get_eos_record("mgo_shen_2026_vinet_3")
    assert f"Warning: {record.warning}" in record.summary()


# --- material lookup ------------------------------------------------------


@pytest.mark.parametrize(
    ("query", "identifier"),
    [
        ("periclase", "mgo"),
        ("Periclase", "mgo"),
        ("MgO", "mgo"),
        ("mgo", "mgo"),
        ("MGO", "mgo"),
        ("gold", "gold"),
        ("Gold", "gold"),
        ("GOLD", "gold"),
        ("Au", "gold"),
        ("au", "gold"),
        ("Co", "cobalt_hcp"),
        ("  Gold ", "gold"),
        ("BRIDGMANITE", "bridgmanite"),
        ("au_fcc", "au_fcc"),
    ],
)
def test_forgiving_material_lookup(query, identifier):
    assert get_material(query).identifier == identifier


def test_exact_identifiers_and_aliases_resolve_exactly_as_before():
    index = _catalog_index()
    for identifier, material in index.materials.items():
        assert get_material(identifier) is material
    for alias, material in index.material_aliases.items():
        assert get_material(alias) is material
    for identifier, material in _compatibility_materials().items():
        if identifier not in index.materials:
            assert get_material(identifier) is material


@pytest.mark.parametrize(
    ("query", "matches"),
    [
        ("Fe", ("fe", "iron")),
        ("Si", ("silicon", "silicon_v", "silicon_vii", "silicon_x")),
        ("h2o", ("ice_vi", "ice_vii")),
    ],
)
def test_ambiguous_formula_lists_candidates_and_never_picks(query, matches):
    with pytest.raises(MaterialLookupError, match="Ambiguous material") as error:
        get_material(query)
    assert error.value.context["matches"] == matches
    message = str(error.value)
    for identifier in matches:
        assert identifier in message
    material = get_material(matches[0])
    assert f"{matches[0]} ({' '.join(material.phase.split())})" in message


def test_uppercase_formula_is_case_sensitive():
    # "CO" is not cobalt; only an all-lowercase query matches formulas loosely.
    with pytest.raises(MaterialLookupError, match="Unknown material 'CO'"):
        get_material("CO")


def test_ambiguous_alias_message_is_bounded():
    with pytest.raises(MaterialLookupError, match="Ambiguous material alias") as error:
        get_material("ferropericlase")
    assert len(error.value.context["matches"]) > 20
    assert "more" in str(error.value)
    assert len(str(error.value)) < 700


def test_document_only_material_is_never_redirected():
    for identifier in ("wadsleyite", "fe_fcc", "coesite_v"):
        with pytest.raises(MaterialLookupError, match="no executable EOS record"):
            get_material(identifier)


@pytest.mark.parametrize(
    ("query", "suggestion"),
    [
        ("periclas", "mgo (Periclase B1)"),
        ("zirconum_alpha", "zirconium_alpha"),
        ("quartz", "alpha_quartz"),
        ("ice", "ice_vi"),
    ],
)
def test_unknown_material_suggestions_are_ranked_and_bounded(query, suggestion):
    with pytest.raises(MaterialLookupError, match="Unknown material") as error:
        get_material(query)
    message = str(error.value)
    suggestions = message.split("Did you mean: ", 1)[1].split("?", 1)[0]
    assert suggestions.startswith(suggestion)
    assert suggestions.count(", ") <= 4
    assert len(message) < 400
    assert error.value.operation == "lookup_material"


# --- EOS record lookup within a material ----------------------------------


def test_record_lookup_stays_exact_and_suggests_briefly():
    material = get_material("mgo")
    with pytest.raises(MaterialLookupError) as error:
        material.get_eos_record("speziale")
    message = str(error.value)
    assert "Did you mean: mgo_speziale_2001_bm3_2?" in message
    assert "summary()" in message
    assert "mgo_ye_2017_vinet_300k" not in message
    assert len(message) < 300
    with pytest.raises(MaterialLookupError) as error:
        material.get_eos_record("MGO_SPEZIALE_2001_BM3_2")
    assert "mgo_speziale_2001_bm3_2" in str(error.value)
    with pytest.raises(MaterialLookupError) as error:
        material.get_eos_record("dewaele")
    assert len(str(error.value)) < 500
    assert error.value.context == {"identifier": "dewaele", "material": "mgo"}


def test_record_lookup_lists_a_small_material_completely():
    material = get_material("au_fcc")
    with pytest.raises(MaterialLookupError, match="available") as error:
        material.get_eos_record("missing")
    for record in material.eos_records:
        assert record.identifier in str(error.value)


def test_record_values_are_not_changed_by_display():
    record = get_eos_record("mgo_speziale_2001_bm3_2")
    before = record.eos.parameter_values()
    record.summary()
    repr(record)
    assert record.eos.parameter_values() == before
    assert np.isfinite(record.pressure(70.0, 300.0))
