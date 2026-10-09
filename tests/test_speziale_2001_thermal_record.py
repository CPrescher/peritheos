"""Published thermal parameters, runtime parity and interchange for Speziale MgO."""

import copy

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_eos_record, get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye, VariableQMieGruneisenDebye
from peritheos.errors import EosValidationError
from scripts.validate_speziale_2001_equations import (
    closed_form_pressure,
    closed_form_theta,
)

IDENTIFIER = "mgo_speziale_2001_bm3_2"


def test_published_parameters_and_author_confirmed_classification():
    doc = get_material_document("mgo")
    raw = next(r for r in doc["eos_records"] if r["identifier"] == IDENTIFIER)
    assert raw["eos"]["parameters"] == {"V0": 74.71, "K0": 160.2, "K0_prime": 3.99}
    assert raw["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 773.0,
        "gamma0": 1.524,
        "q0": 1.65,
        "q1": 11.8,
        "n": 2.0,
    }
    assert raw["thermal"]["parameter_errors"]["q1"] == 0.2
    assert raw["equation_kind"] == "thermal"
    assert raw["scientific_validation"]["status"] == "primary_source_validated"
    assert raw["scientific_validation"]["reproduction_status"] == "partially_reproduced"
    assert (
        raw["scientific_validation"]["related_thermal_reconstruction"][
            "formulation_status"
        ]
        == "author_confirmed"
    )
    assert get_eos_record(IDENTIFIER).is_default
    Draft202012Validator(eosmat_schema()).validate(doc)


def test_runtime_and_python_fallback_against_independent_analytic_audit():
    record = get_eos_record(IDENTIFIER)
    assert isinstance(record.eos, VariableQMieGruneisenDebye)
    fallback = copy.copy(record.eos)
    del fallback._native
    volumes = 74.71 * np.array([0.6, 0.64, 0.7, 0.8, 0.9, 1.0, 1.02])
    temperatures = np.array([300.0, 500.0, 1100.0, 2000.0, 3663.0])
    model_volumes = volumes[:, None] * record.volume_scale
    expected = np.array(
        [[closed_form_pressure(v, t) for t in temperatures] for v in volumes]
    )
    assert record.pressure(volumes[:, None], temperatures) == pytest.approx(
        expected, abs=1e-8
    )
    assert fallback.pressure(model_volumes, temperatures) == pytest.approx(
        expected, abs=1e-8
    )
    assert record.eos.characteristic_temperature(model_volumes[:, 0]) == pytest.approx(
        closed_form_theta(volumes), abs=1e-7
    )
    assert fallback.characteristic_temperature(model_volumes[:, 0]) == pytest.approx(
        closed_form_theta(volumes), abs=1e-9
    )
    # Thermal pressure subtraction preserves every point on the 300 K cold curve.
    assert record.eos.thermal_pressure(model_volumes, 300.0) == pytest.approx(
        np.zeros_like(model_volumes), abs=1e-12
    )


def test_thermal_interchange_and_volume_inversion_preserve_molar_normalization():
    material = Material.from_eosmat(
        get_material_document("mgo"), record_identifiers=[IDENTIFIER]
    )
    doc = material.to_eosmat()
    Draft202012Validator(eosmat_schema()).validate(doc)
    reloaded = Material.from_eosmat(doc).get_eos_record(IDENTIFIER)
    for v, t in [(74.71, 1100.0), (59.768, 2000.0), (47.8144, 3663.0)]:
        expected = closed_form_pressure(v, t)
        assert reloaded.pressure(v, t) == pytest.approx(expected, abs=1e-8)
        assert reloaded.volume(expected, t) == pytest.approx(v, abs=1e-8)
    assert doc["eos_records"][0]["thermal"]["parameters"] == {
        "Tr": 300.0,
        "theta0": 773.0,
        "gamma0": 1.524,
        "q0": 1.65,
        "q1": 11.8,
        "n": 2.0,
    }


@pytest.mark.parametrize("q0,q1", [(1.65, 0.0), (0.0, 11.8), (1.65, 1e-12)])
def test_constant_q_and_constant_gamma_limits(q0, q1):
    reference = BM3(1.124785340449, 160.2, 3.99)
    variable = VariableQMieGruneisenDebye(reference, 300, 773, 1.524, q0, q1, 2)
    constant = MieGruneisenDebye(reference, 300, 773, 1.524, q0, 2)
    volumes = reference.V0 * np.array([0.64, 0.8, 1.0, 1.02])
    for candidate in (variable, copy.copy(variable)):
        if candidate is not variable:
            del candidate._native
        assert candidate.characteristic_temperature(volumes) == pytest.approx(
            constant.characteristic_temperature(volumes), rel=1e-10
        )
        assert candidate.pressure(volumes, 2000) == pytest.approx(
            constant.pressure(volumes, 2000), abs=1e-8
        )


def test_gamma_theta_derivative_and_invalid_states():
    model = get_eos_record(IDENTIFIER).eos
    v, h = model.rt_eos.V0 * 0.8, 1e-5
    inferred_gamma = -(
        np.log(model.characteristic_temperature(v * np.exp(h)))
        - np.log(model.characteristic_temperature(v * np.exp(-h)))
    ) / (2 * h)
    assert inferred_gamma == pytest.approx(model.gruneisen_parameter(v), rel=1e-8)
    for invalid in (0, -1, np.nan, np.inf):
        with pytest.raises(EosValidationError):
            model.characteristic_temperature(invalid)
        with pytest.raises(EosValidationError):
            model.pressure(v, invalid)


def test_native_fitting_recovers_gamma_from_independent_pressure_targets(monkeypatch):
    from peritheos import _rust
    from peritheos.fitting import fit_thermal_eos

    record = get_eos_record(IDENTIFIER)
    volumes = 74.71 * np.array([0.64, 0.8, 0.9, 1.0])
    temperatures = np.array([1100.0, 2000.0, 3000.0, 3663.0])
    targets = np.array(
        [closed_form_pressure(v, t) for v, t in zip(volumes, temperatures)]
    )

    def forbid_callback_solver(*args, **kwargs):
        raise AssertionError("the built-in model must use native fitting")

    monkeypatch.setattr(_rust, "fit_least_squares", forbid_callback_solver)
    result = fit_thermal_eos(
        VariableQMieGruneisenDebye,
        record.eos.rt_eos,
        volumes * record.volume_scale,
        temperatures,
        targets,
        initial={"gamma0": 1.4},
        fixed={"Tr": 300.0, "theta0": 773.0, "q0": 1.65, "q1": 11.8, "n": 2.0},
    )
    assert result.success
    assert result.parameters["gamma0"] == pytest.approx(1.524, abs=1e-7)
