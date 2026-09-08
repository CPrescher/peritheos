"""Regression tests for the Mosenfelder et al. (2009) PPv audit."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
SCRIPT = ROOT / "scripts/reproduce_mosenfelder_2009_mgsio3_post_perovskite.py"


def reproduction_module():
    spec = importlib.util.spec_from_file_location("mosenfelder_reproduction", SCRIPT)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_source_inputs_and_derived_quantities_are_explicitly_distinguished():
    document = get_material_document("mgsio3_post_perovskite")
    datasets = {item["identifier"]: item for item in document["datasets"]}
    static = datasets["mgsio3_post_perovskite_guignot_2007_table1_pvt"]
    shock = datasets["mgsio3_post_perovskite_mosenfelder_2009_table2_shock"]
    derived = datasets[
        "mgsio3_post_perovskite_mosenfelder_2009_derived_reference_isentrope"
    ]

    assert static["resource"]["sha256"] == (
        "16e372df6ecaae23df2b9ec0403522d472d85248eb65c26a9f332780a6095f0c"
    )
    statuses = {
        column["name"]: column.get("measurement_status") for column in shock["columns"]
    }
    assert statuses["shock_velocity_km_s"] == "measured_shock_transit"
    assert statuses["pressure_gpa"] == (
        "derived_by_impedance_matching_and_rankine_hugoniot"
    )
    assert statuses["shock_density_mg_m3"] == (
        "derived_by_rankine_hugoniot_mass_conservation"
    )
    assert statuses["hugoniot_temperature_k"] == (
        "measured_optical_pyrometry_where_reported"
    )
    assert derived["kind"] == "derived_model_grid"
    assert "not a fit input or experimental dataset" in derived["notes"]
    assert derived["resource"]["sha256"] == (
        "f3f98e61c3d6d322cc95d63168410cf47ec2a3005169580ba44b3246f7890626"
    )


def test_record_executes_published_reference_isentrope_and_heat_capacity():
    document = get_material_document("mgsio3_post_perovskite")
    record = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == "mgsio3_post_perovskite_mosenfelder_2009_bm3_1"
    )
    assert record["thermal"]["thermal_pressure_reference"] == "reference_isentrope"
    assert record["thermal"]["parameters"]["Cvmax"] == pytest.approx(1.035 * 100.3875)
    executable = Material.from_eosmat(document).get_eos_record(record["identifier"])
    assert executable.eos.parameter_values(include_reference=False)["Cvmax"] == (
        pytest.approx(103.9010625)
    )
    changed = executable.eos.with_parameters(Cvmax=100.0)
    assert changed.Cvmax == pytest.approx(100.0)
    assert executable.pressure(123.0, 300.0) == pytest.approx(111.01452694098933)
    assert executable.pressure(119.7, 2520.0) == pytest.approx(144.35293739663)


def test_errors_in_variables_reproduction_recovers_table_4_and_chi_square():
    module = reproduction_module()
    static, shock = module.arrays()
    assert static["pressure_gpa"].size == 48
    assert shock["pressure_gpa"].size == 6
    assert np.count_nonzero(np.isfinite(shock["hugoniot_temperature_k"])) == 3

    result = module.fit_orthogonal(static, shock)
    fitted = result.x[:6]
    assert fitted == pytest.approx(
        [224.55806, 4.219879, 2.518345, 2.067256, 1.047205, 998.67485],
        rel=2.0e-5,
    )
    assert 2.0 * result.cost / 48 == pytest.approx(0.18004566, rel=2.0e-5)
    assert module.static_pressure(123.0, 300.0, module.PUBLISHED) == pytest.approx(
        111.01452694098933
    )


def test_shock_temperature_is_not_required_to_close_solid_hugoniot():
    module = reproduction_module()
    _static, shock = module.arrays()
    pressure = module.shock_pressure(
        shock["shock_density_mg_m3"],
        shock["initial_density_mg_m3"],
        shock["transition_energy_j_kg"],
        module.PUBLISHED,
    )
    assert np.all(np.isfinite(pressure))
    assert pressure.size == 6
