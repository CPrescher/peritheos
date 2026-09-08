import numpy as np
import pytest

from peritheos import get_material_document
from scripts.reproduce_noguchi_1999_nio import (
    SOURCE_CONFERENCE_K0_ERROR_GPA,
    SOURCE_CONFERENCE_K0_GPA,
    SOURCE_JOURNAL_K0_GPA,
    SOURCE_JOURNAL_K0_PRIME,
    fit_conference_crosscheck,
    fit_journal_isotherm,
)

RECORD_ID = "nickel_oxide_noguchi_1999_bm3_1"


def _source_record():
    document = get_material_document("nickel_oxide")
    return next(
        record
        for record in document["eos_records"]
        if record["identifier"] == RECORD_ID
    )


def test_noguchi_reduction_metadata_preserves_all_source_conventions():
    source = _source_record()
    validation = source["scientific_validation"]["independent_refit"]

    assert source["fit_datasets"] == ["nickel_oxide_noguchi_1999_table1_shock"]
    assert validation["classification"] == "similar"
    assert validation["thermal_parameters"] == {
        "reference_temperature_k": 300.0,
        "debye_temperature_k": 390.0,
        "gamma0": 1.38,
        "gamma_volume_exponent_q": 1.0,
        "atoms_per_formula_unit": 2.0,
        "debye_temperature_law": "integrated_gruneisen",
        "gamma_volume_law": "gamma_over_volume_constant",
    }
    assert validation["hugoniot_parameters"] == {
        "c0_km_s": 5.36,
        "s": 1.19,
        "rho0_g_cm3": 6.781,
    }
    assert "must never be fitted directly as an isotherm" in source["notes"]


def test_noguchi_conference_reduction_reproduces_both_published_checks():
    states, fit = fit_conference_crosscheck()

    assert states.volume_ratio.size == 7
    assert fit.parameters["K0_prime"] == 4.0
    assert fit.parameters["K0"] == pytest.approx(184.2974466, rel=1.0e-8)
    assert abs(fit.parameters["K0"] - SOURCE_CONFERENCE_K0_GPA) < (
        SOURCE_CONFERENCE_K0_ERROR_GPA
    )
    assert states.shock_temperature_k[-1] == pytest.approx(1618.660372, rel=1.0e-8)
    assert states.hugoniot_pressure_gpa[-1] == 132.9


def test_noguchi_journal_refit_uses_thermally_reduced_pressures_only():
    states, fit = fit_journal_isotherm()

    assert states.volume_ratio.size == 8
    assert np.all(states.shock_temperature_k > 300.0)
    assert np.all(states.thermal_pressure_gpa > 0.0)
    assert np.all(states.isothermal_pressure_gpa < states.hugoniot_pressure_gpa)
    assert states.shock_temperature_k[-1] == pytest.approx(1784.715643, rel=1.0e-8)
    assert states.isothermal_pressure_gpa[-1] == pytest.approx(138.6051987, rel=1.0e-8)

    assert fit.parameters["K0"] == pytest.approx(182.4351120, rel=1.0e-8)
    assert fit.parameters["K0_prime"] == pytest.approx(4.12465196, rel=1.0e-8)
    assert np.sqrt(np.mean(fit.residuals**2)) == pytest.approx(2.08431895, rel=1.0e-8)
    assert (
        abs(fit.parameters["K0"] - SOURCE_JOURNAL_K0_GPA) < (fit.standard_errors["K0"])
    )
    assert (
        abs(fit.parameters["K0_prime"] - SOURCE_JOURNAL_K0_PRIME)
        < (fit.standard_errors["K0_prime"])
    )
