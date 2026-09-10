import copy

import numpy as np
import pytest
from scipy.constants import Avogadro, R

from peritheos import get_material_document
from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import SoundVelocityDebyeHelmholtz
from peritheos.materials import Material


def model():
    return SoundVelocityDebyeHelmholtz(
        Vinet(1.12, 169.8, 4.501),
        Tr=300.0,
        molar_mass_g_mol=40.304,
        n=2.0,
        longitudinal_intercept=-2.546415,
        longitudinal_slope=3.43687325,
        shear_intercept=-0.02936402,
        shear_slope=1.69849549,
    )


def test_sound_velocity_debye_implements_appendix_a_relations():
    eos = model()
    volume = 1.0
    density = 40.304 / 10.0
    fitted_longitudinal = -2.546415 + 3.43687325 * density
    fitted_shear = -0.02936402 + 1.69849549 * density
    tau_squared = (fitted_longitudinal / fitted_shear) ** 2
    poisson = (tau_squared - 2.0) / (2.0 * (tau_squared - 1.0))

    assert eos.density(volume) == pytest.approx(density)
    assert eos.poisson_ratio(volume) == pytest.approx(poisson)
    assert eos.longitudinal_velocity(volume) > eos.shear_velocity(volume) > 0.0
    assert eos.characteristic_temperature(volume) > 0.0
    assert eos.gruneisen_parameter(volume) > 0.0


def test_sound_velocity_debye_pressure_is_helmholtz_derivative():
    eos = model()
    volume = 0.9
    temperature = 4000.0
    step = volume * 2.0e-5
    numerical_pressure = -(
        eos.thermal_helmholtz_free_energy(volume + step, temperature)
        - eos.thermal_helmholtz_free_energy(volume - step, temperature)
    ) / (2.0 * step) / 1.0e4

    assert eos.thermal_pressure(volume, temperature) == pytest.approx(
        numerical_pressure, rel=2.0e-5
    )
    assert eos.thermal_pressure_increment(volume, 300.0) == pytest.approx(0.0)
    assert eos.molar_heat_capacity_v(volume, 1.0e6) == pytest.approx(
        3.0 * 2.0 * R, rel=2.0e-4
    )


def test_sound_velocity_debye_hugoniot_energy_solver():
    eos = model()
    initial_volume = 40.304 / 3.590 / 10.0
    volume = 0.75
    pressure = 200.0
    temperature = eos.hugoniot_temperature_from_state(
        volume, pressure, initial_volume=initial_volume
    )
    initial_energy = eos.internal_energy(initial_volume, 300.0)
    target = initial_energy + 0.5 * pressure * (initial_volume - volume) * 1.0e4

    assert temperature > 300.0
    assert eos.internal_energy(volume, temperature) == pytest.approx(target)


def test_sound_velocity_debye_eosmat_registration_and_molar_scaling():
    document = get_material_document("mgo")
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == "mgo_b1_luo_2023_vinet_thermal_5"
    )
    record = copy.deepcopy(source)
    record["thermal"] = {
        "type": "SoundVelocityDebyeHelmholtz",
        "model": "sound_velocity_debye_helmholtz",
        "parameters": {
            "Tr": 300.0,
            "molar_mass_g_mol": 40.304,
            "n": 2.0,
            "longitudinal_intercept": -2.546415,
            "longitudinal_slope": 3.43687325,
            "shear_intercept": -0.02936402,
            "shear_slope": 1.69849549,
        },
    }
    record["volume"]["public_to_model_scale"] = Avogadro * 1.0e-25 / 4.0
    synthetic = copy.deepcopy(document)
    synthetic["eos_records"] = [
        record if item["identifier"] == record["identifier"] else item
        for item in synthetic["eos_records"]
    ]

    loaded = Material.from_eosmat(
        synthetic, record_identifiers=[record["identifier"]]
    ).eos_records[0]
    assert isinstance(loaded.eos, SoundVelocityDebyeHelmholtz)
    assert loaded.eos.rt_eos.V0 == pytest.approx(
        source["eos"]["parameters"]["V0"] * Avogadro * 1.0e-25 / 4.0
    )
    assert np.isfinite(loaded.pressure(60.0, 3000.0))
