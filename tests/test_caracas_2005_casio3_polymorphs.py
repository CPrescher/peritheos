"""Primary-source checks for Caracas et al. (2005) CaSiO3 fits."""

from __future__ import annotations

import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_caracas_2005_casio3_polymorphs import PARAMETERS, reproduce

MATERIALS = {
    "Pm-3m": "ca_perovskite",
    "I4/mcm": "ca_perovskite_tetragonal",
    "Imma": "casio3_perovskite_imma",
    "R-3c": "casio3_perovskite_r3c",
    "P4/mbm": "casio3_perovskite_p4mbm",
    "I4/mmm": "casio3_perovskite_i4mmm",
    "Im-3": "casio3_perovskite_im3",
    "P42/nmc": "casio3_perovskite_p42nmc",
    "Pnma": "casio3_perovskite_pnma",
}


@pytest.mark.parametrize(("phase", "material"), MATERIALS.items())
def test_both_source_parameterizations_are_executable(phase, material):
    document = get_material_document(material)
    identifiers = [
        record["identifier"]
        for record in document["eos_records"]
        if record.get("reference", {}).get("doi") == "10.1029/2004GL022144"
    ]
    loaded = Material.from_eosmat(document, record_identifiers=identifiers)
    assert len(loaded.eos_records) == 2
    assert {record.eos.__class__.__name__ for record in loaded.eos_records} == {
        "BM3",
        "BM4",
    }
    for record in loaded.eos_records:
        assert record.eos.pressure(record.eos.V0) == pytest.approx(0.0)
        assert record.eos.bulk_modulus(record.eos.V0) == pytest.approx(record.eos.K0)


def test_table2_contains_nine_complete_bm3_bm4_pairs():
    assert set(PARAMETERS) == set(MATERIALS)
    assert sum(len(pair) for pair in PARAMETERS.values()) == 18


def test_independent_density_statements_select_the_bm3_curves():
    metrics = reproduce()
    assert metrics["Pm-3m"]["bm3_density_0_g_cm3"] == pytest.approx(4.32, abs=0.01)
    assert metrics["Pm-3m"]["bm3_density_130_g_cm3"] == pytest.approx(5.77, abs=0.01)
    assert metrics["I4/mcm"]["bm3_density_0_g_cm3"] == pytest.approx(4.33, abs=0.01)
    assert metrics["I4/mcm"]["bm3_density_130_g_cm3"] == pytest.approx(5.78, abs=0.01)


def test_every_new_phase_card_is_indexing_only_without_invented_coordinates():
    for material in tuple(MATERIALS.values())[2:]:
        document = get_material_document(material)
        assert document["formula_units_per_cell"] is None
        assert "lattice" not in document
        assert "atom_sites" not in document
