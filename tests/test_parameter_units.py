"""Every stored coefficient in the bundled catalog shows a unit.

The display unit (``show`` and ``overview()``) must state what the code
assumes for the stored value. The checks below compare a few units with
independent facts in the implementation.
"""

from __future__ import annotations

import inspect

import numpy as np
import pytest
from scipy.constants import Avogadro

from peritheos import materials
from peritheos._overview import parameter_unit
from peritheos.eos.thermal.dewaele import Dewaele2006
from peritheos.eos.thermal.dorogokupets_oganov import DorogokupetsOganov2007
from peritheos.materials import get_eos_record, list_materials

MOLAR_VOLUME = "J bar^-1 mol^-1"


def _unit(model: str, name: str) -> str | None:
    return parameter_unit(
        model,
        name,
        volume_unit="angstrom^3/conventional_unit_cell",
        pressure_unit="GPa",
    )


def _thermal_type(record_id: str) -> str:
    return get_eos_record(record_id).overview()["thermal_model"]["identifier"]


def _stored(record_id: str, component: str, name: str) -> float:
    rows = get_eos_record(record_id).overview()["parameters"]
    (value,) = [
        row["value"]
        for row in rows
        if row["component"] == component and row["name"] == name
    ]
    return value


def test_every_catalog_coefficient_with_a_value_has_a_unit():
    missing = set()
    for material in list_materials():
        for record in material.eos_records:
            for row in record.overview()["parameters"]:
                if row["value"] is not None and row["unit"] is None:
                    missing.add((record.identifier, row["component"], row["name"]))
    assert not missing, sorted(missing)


def test_double_debye_volumes_are_molar_and_never_converted():
    # The model is molar, so the loader converts only the reference V0 to
    # J bar^-1 mol^-1; Vp and Ve are used as stored, so they are molar too.
    record_id = "diamond_benedict_2014_dewaele_anchored"
    assert _thermal_type(record_id) == "DoubleDebyeHelmholtz"
    assert "double_debye_helmholtz" in materials._MOLAR_VOLUME_THERMAL_MODELS
    record = get_eos_record(record_id)
    assert record.eos.Vp == _stored(record_id, "thermal_model", "Vp")
    material = materials.get_material("diamond")
    scale = Avogadro * 1.0e-25 / material.formula_units_per_cell
    stored_v0 = _stored(record_id, "model", "V0")
    assert record.eos.rt_eos.V0 == pytest.approx(stored_v0 * scale, rel=1e-12)
    for name in ("Vp", "Ve"):
        assert _unit("DoubleDebyeHelmholtz", name) == MOLAR_VOLUME
    assert _unit("DoubleDebyeHelmholtz", "a_a") == f"({MOLAR_VOLUME})^-1"


def test_alpha_kt_keeps_the_record_volume():
    # AlphaKT is not a molar model: V0 stays in the record unit and the
    # coefficients are per kelvin and GPa, with no volume in their unit.
    record_id = "aragonite_martinez_1996_bm2_2"
    assert _thermal_type(record_id) == "AlphaKT"
    assert "thermal_reference_state" not in materials._MOLAR_VOLUME_THERMAL_MODELS
    record = get_eos_record(record_id)
    assert record.eos.rt_eos.V0 == _stored(record_id, "model", "V0")
    assert _unit("AlphaKT", "alpha0") == "K^-1"
    assert _unit("AlphaKT", "alpha1") == "K^-2"
    assert _unit("AlphaKT", "dK_dT") == "GPa/K"
    assert _unit("AlphaKT", "beta1") == "GPa^-1 K^-1"


def test_dorogokupets_oganov_electronic_e_is_in_micro_per_kelvin():
    # The same electronic pressure term appears in Dewaele2006 with the
    # coefficient in K^-1. DO2007 matches it only when its stored value is
    # read as 10^-6 K^-1.
    record_id = "platinum_dorogokupets_oganov_2007_vinet_4"
    assert _thermal_type(record_id) == "DorogokupetsOganov2007"
    source = get_eos_record(record_id).eos
    names = list(inspect.signature(DorogokupetsOganov2007).parameters)
    values = {name: getattr(source, name) for name in names if name != "rt_eos"}
    electronic_e = 50.0

    def do2007(e: float) -> DorogokupetsOganov2007:
        return DorogokupetsOganov2007(source.rt_eos, **{**values, "electronic_e": e})

    volume = 0.95 * source.rt_eos.V0
    temperature = 2000.0
    do_increment = do2007(electronic_e).thermal_pressure(volume, temperature) - do2007(
        0.0
    ).thermal_pressure(volume, temperature)
    dewaele = Dewaele2006(
        source.rt_eos,
        Tr=values["Tr"],
        theta0=values["theta_E1"],
        gamma0=values["gamma0"],
        gamma_inf=values["gamma_inf"],
        beta=values["beta"],
        anharmonic_a=0.0,
        anharmonic_m=1.0,
        electronic_e=electronic_e * 1.0e-6,
        electronic_g=values["electronic_g"],
        n=values["n"],
    )
    expected = dewaele.electronic_pressure_increment(volume, temperature)
    assert np.isclose(do_increment, expected, rtol=1e-9, atol=0.0)
    assert _unit("DorogokupetsOganov2007", "electronic_e") == "10^-6 K^-1"
    assert _unit("Dewaele2006", "electronic_e") == "K^-1"


def test_multi_oscillator_scaled_coefficients_name_their_factor():
    record_id = "aluminum_sokolova_2013_holzapfel_2"
    assert _thermal_type(record_id) == "Sokolova2016"
    for model in ("Sokolova2016", "MultiOscillatorGruneisen"):
        assert _unit(model, "a_0") == "10^-6 K^-1"
        assert _unit(model, "e_0") == "10^-6 K^-1"
        assert _unit(model, "QE1o") == "K"


def test_linear_hugoniot_slope_is_dimensionless():
    # Us = c0 + s up with Us, up and c0 all in km/s.
    record = get_eos_record("cementite_nonmagnetic_hu_2019_linear_hugoniot")
    hugoniot = record.eos
    up = 2.0
    shock = hugoniot.shock_velocity_from_particle_velocity(up)
    assert shock == pytest.approx(hugoniot.c0 + hugoniot.s * up)
    assert _unit("LinearUsUpHugoniot", "c0") == "km/s"
    assert _unit("LinearUsUpHugoniot", "s") == "1"
