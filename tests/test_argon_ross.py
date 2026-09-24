"""Protect precursor identity, observations and honest non-reproduction status."""

import json

import numpy as np
import pytest

from peritheos import Material, get_material_document
from scripts.reproduce_argon_ross import (
    REPORT,
    audit,
    precursor_exp6_energy_kelvin,
)


def test_ross_supporting_study_does_not_fabricate_a_journal_eos():
    document = get_material_document("argon_fcc")
    study = next(
        s
        for s in document["supporting_studies"]
        if s["identifier"] == "argon_ross_1986_supporting_study"
    )
    assert study["reference"]["doi"] == "10.1063/1.451346"
    assert study["reproduction_status"] == "not_reproduced"
    assert study["requested_publication_pdf_status"] == "not_obtained"
    assert study["eos_record_identifiers"] == []
    assert not any("ross" in r["identifier"] for r in document["eos_records"])
    material = Material.from_eosmat(document)
    assert material.to_eosmat()["supporting_studies"] == document["supporting_studies"]
    table = material.get_dataset("argon_ross_1985_precursor_figure1_subset")
    assert table.checksum_verified
    assert table.used_by_eos_records == ()
    assert table.reference["year"] == 1985
    assert "doi" not in table.reference
    pv = table.as_pressure_volume()
    assert pv.pressure_sigma is None and pv.volume_sigma is None
    np.testing.assert_equal(table["temperature_k"], np.full(19, 293))
    # Independent conversion: 1 A3/atom = 0.602214076 cm3/mol of atoms.
    np.testing.assert_allclose(
        pv.volume / 4 * 0.602214076,
        table["molar_volume_cm3_mol"],
        atol=6e-6,
    )


def test_precursor_audit_is_reproducible_without_claiming_an_eos_fit():
    result = audit()
    assert result == json.loads(REPORT.read_text())
    assert result["conversion_max_absolute_roundoff"] < 5e-5
    assert result["eos_fit_performed"] is False
    assert result["pressure_residual_rmse_gpa"] is None
    assert result["matching_journal_pdf_obtained"] is False


@pytest.mark.parametrize("alpha", [13.0, 13.2])
def test_precursor_potential_has_printed_well_depth_and_stationary_minimum(alpha):
    def potential(r):
        return precursor_exp6_energy_kelvin(r, alpha)

    assert potential(3.85) == pytest.approx(-122)
    h = 1e-5
    assert abs((potential(3.85 + h) - potential(3.85 - h)) / (2 * h)) < 1e-6
    assert potential(3.85 - h) > potential(3.85)
    assert potential(3.85 + h) > potential(3.85)
    with pytest.raises(ValueError):
        potential(1.0)
