"""Recommendation gates and source-observation selection, without changing source values."""

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material, get_material_document
from peritheos.errors import MaterialError, MaterialLookupError
from scripts.apply_eos_recommendation_holds import NOT_RECOMMENDED, WADSLEYITE
from scripts.validate_primary_eos_refits import (
    DEWAELE_REFIT_JSON,
    _dewaele_2019_outcome,
    _fit_record,
    _series,
    _solomatova_outcome,
)

ROOT = Path(__file__).resolve().parents[1]


def source(identifier):
    for path in (ROOT / "peritheos/data/materials").glob("*.eosmat"):
        document = json.loads(path.read_text())
        for record in document["eos_records"]:
            if record["identifier"] == identifier:
                return document, record
    raise AssertionError(identifier)


@pytest.mark.parametrize("identifier", NOT_RECOMMENDED)
def test_seven_unreproduced_eos_require_explicit_scientific_acceptance(identifier):
    document, raw = source(identifier)
    record = get_eos_record(identifier)
    assert record.scientific_validation_status == "not_reproduced"
    assert (
        raw["scientific_validation"]["usage_recommendation"]
        == "not_recommended_for_quantitative_use"
    )
    assert raw["catalog_access"] == "explicit_selection"
    assert (
        get_material(document["identifier"]).default_record().identifier != identifier
    )
    strict = Material.from_eosmat(document)
    assert identifier not in {r.identifier for r in strict.eos_records}
    with pytest.raises(MaterialError, match="not_reproduced"):
        Material.from_eosmat(document, record_identifiers=[identifier])
    accepted = Material.from_eosmat(
        document, record_identifiers=[identifier], require_primary_validation=False
    )
    assert accepted.eos_records[0].eosmat_metadata["eos"] == raw["eos"]
    assert accepted.eos_records[0].eosmat_metadata.get("thermal") == raw.get("thermal")


def test_held_wadsleyite_is_document_only_without_a_replacement_default():
    doc = get_material_document("wadsleyite")
    assert doc["recommendation_hold"]["record"] == WADSLEYITE
    assert len(doc["eos_records"]) == 1
    assert doc["eos_records"][0]["scientific_validation"]["status"] == "deferred"
    with pytest.raises(MaterialLookupError):
        get_material("wadsleyite")
    with pytest.raises(MaterialLookupError):
        get_eos_record(WADSLEYITE)
    with pytest.raises(MaterialError, match="deferred"):
        Material.from_eosmat(doc)
    diagnostic = Material.from_eosmat(doc, require_primary_validation=False)
    assert (
        diagnostic.eos_records[0].eosmat_metadata["eos"]["parameters"]
        == doc["eos_records"][0]["eos"]["parameters"]
    )


@pytest.mark.parametrize(
    "identifier,column,count,rms",
    [
        (
            "ferropericlase_fe27_jacobsen_2002_bm3_1",
            "lattice_a_angstrom",
            10,
            0.005738353913,
        ),
        (
            "magnesiowustite_fe56_jacobsen_2002_bm3_1",
            "lattice_a_angstrom",
            11,
            0.012872784351,
        ),
        (
            "magnesiowustite_fe75_jacobsen_2002_bm3_1",
            "lattice_a_angstrom",
            7,
            0.005057308845,
        ),
        (
            "mg073fe027o_jacobsen_2005_bm3_1",
            "mg073fe027o_derived_volume_a3",
            16,
            0.193911160661,
        ),
        (
            "fe093o_b1_jacobsen_2005_bm3_1",
            "fe093o_b1_derived_volume_a3",
            7,
            0.119998266979,
        ),
        (
            "mg0215fe0762vac0023o_finkelstein_2017_bm3_helium_1",
            "cubic_volume_a3",
            17,
            0.277922409860,
        ),
        (
            "mg0215fe0762vac0023o_finkelstein_2017_bm3_neon_cubic_2",
            "cubic_volume_a3",
            5,
            0.118251289423,
        ),
        (
            "mgfe94o_b1_dobrosavljevic_2019_bm3_1",
            "atomic_volume_angstrom3",
            14,
            0.240561470879,
        ),
        (
            "mgfe94o_rhombohedral_dobrosavljevic_2019_bm3_1",
            "atomic_volume_angstrom3",
            12,
            0.728146404152,
        ),
        (
            "e_feooh_hc_low_spin_thompson_2017_bm3_1",
            "unit_cell_volume_a3",
            12,
            0.108193702672,
        ),
    ],
)
def test_ten_audits_use_sample_observable_and_source_branch(
    identifier, column, count, rms
):
    doc, raw = source(identifier)
    dataset = next(d for d in doc["datasets"] if d["identifier"] in raw["fit_datasets"])
    series = _series(doc, raw, dataset)
    assert series.volume_column == column
    assert len(series.pressure) == count
    model = get_eos_record(identifier)
    residuals = model.eos.pressure(series.volume) - series.pressure
    assert float(np.sqrt(np.mean(residuals**2))) == pytest.approx(rms, abs=1e-9)
    if "jacobsen_2002" in identifier:
        assert np.all(series.volume < 80)  # sample cell, not quartz calibrant
        assert series.volume_sigma is not None


def test_gong_does_not_fit_a_reference_isentrope_to_hugoniot_states():
    identifier = "mg092fe008sio3_bridgmanite_gong_2004_adiabatic_bm3"
    doc, raw = source(identifier)
    dataset = next(d for d in doc["datasets"] if d["identifier"] in raw["fit_datasets"])
    outcome = _fit_record(doc, raw, dataset)
    assert outcome["status"] == "not_refittable"
    assert "Hugoniot" in outcome["reason"] and "isentrope" in outcome["reason"]
    assert "published_rmse_gpa" not in outcome


@pytest.mark.parametrize("suffix", ["mao", "dor"])
def test_dewaele_iron_uses_existing_complete_recalibrated_hcp_audit(suffix):
    identifier = f"iron_dewaele_2019_{suffix}_vinet"
    doc, raw = source(identifier)
    refit = json.loads(DEWAELE_REFIT_JSON.read_text())["row_level_refits"][identifier]
    outcome = _dewaele_2019_outcome(doc, raw, refit)
    assert outcome["observations"] == 53
    assert outcome["status"] == "similar"
    assert all(p["within_reported_95pct"] for p in outcome["parameters"])


@pytest.mark.parametrize("branch", ["hs", "ls"])
def test_fp48_uses_coupled_spin_crossover_audit(branch):
    identifier = f"mg0490fe0483ti0027o_solomatova_2016_{branch}_bm3_reference_{1 if branch == 'hs' else 2}"
    _, raw = source(identifier)
    outcome = _solomatova_outcome(raw)
    assert outcome["fit_kind"] == "coupled_bm3_three_level_spin_crossover"
    assert outcome["observations"] == 45
    assert outcome["status"] == "similar"
    assert outcome["rmse_volume_a3"] == pytest.approx(0.1198, abs=5e-5)


@pytest.mark.parametrize("identifier", [*NOT_RECOMMENDED, WADSLEYITE])
def test_generic_transcription_reaudit_cannot_undo_reviewed_recommendation(identifier):
    from scripts.apply_primary_source_audit import audit_record

    document, record = source(identifier)
    restored = audit_record(record, document["identifier"] + ".eosmat")
    assert restored == record
    assert restored is not record
