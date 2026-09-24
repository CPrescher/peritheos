"""Independent source checks for Xiao's 2025 Helmholtz EOS."""

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.eos.rt import NaturalStrain4
from peritheos.eos.thermal import DebyeAnharmonicHelmholtz
from peritheos.errors import EosValidationError


def model(cls=DebyeAnharmonicHelmholtz):
    k0 = 2.6565
    kp = 2 + 2 * 7.298 / k0
    kpp = (6 * 0.010 / k0 - kp**2 + 3 * kp - 3) / k0
    return cls(
        NaturalStrain4(2.2555, k0, kp, kpp),
        70,
        86.44,
        2.68,
        0.0024,
        0.0128,
        0.388,
        7.85,
    )


def test_authors_workbook_caloric_and_mechanical_state():
    eos = model()
    # Cached source workbook cells at T=70 K and Vm=23 cm³/mol.
    assert eos.pressure(2.3, 70) == pytest.approx(0.0783517375843629, abs=1e-12)
    assert eos.helmholtz_free_energy(2.3, 70) == pytest.approx(
        -998.2296606304144, abs=1e-8
    )
    assert eos.internal_energy(2.3, 70) == pytest.approx(1098.2008906760243, abs=1e-8)
    assert eos.entropy(2.3, 70) == pytest.approx(29.949007875806267, abs=1e-10)
    assert eos.molar_heat_capacity_v(2.3, 70) == pytest.approx(
        22.854917560882143, abs=1e-10
    )
    assert eos.bulk_modulus(2.3, 70) == pytest.approx(
        0.001 / 0.0004448515520965651, abs=1e-8
    )
    assert eos.thermal_pressure_increment(2.3, 70) == 0
    assert eos.thermal_pressure(2.3, 70) > 0


def test_vector_roundtrips_fallback_and_energy_derivatives():
    class Fallback(DebyeAnharmonicHelmholtz):
        pass

    native, fallback = model(), model(Fallback)
    v = np.array([2.3, 1.8, 1.51, 2.2])
    t = np.array([70.0, 300.0, 760.0, 1.0])
    for quantity in [
        "pressure",
        "helmholtz_free_energy",
        "internal_energy",
        "entropy",
        "molar_heat_capacity_v",
        "molar_heat_capacity_p",
    ]:
        assert np.allclose(
            getattr(native, quantity)(v, t),
            getattr(fallback, quantity)(v, t),
            rtol=1e-7,
            atol=1e-7,
        )
    assert np.allclose(native.volume(native.pressure(v, t), t), v, rtol=1e-10)
    dv, dt = 1e-5, 1e-3
    pressure = -(
        native.helmholtz_free_energy(v + dv, t)
        - native.helmholtz_free_energy(v - dv, t)
    ) / (2 * dv * 1e4)
    entropy = -(
        native.helmholtz_free_energy(v, t + dt)
        - native.helmholtz_free_energy(v, t - dt)
    ) / (2 * dt)
    cv = (native.internal_energy(v, t + dt) - native.internal_energy(v, t - dt)) / (
        2 * dt
    )
    assert np.allclose(pressure, native.pressure(v, t), atol=1e-8)
    assert np.allclose(entropy, native.entropy(v, t), atol=1e-7)
    assert np.allclose(cv, native.molar_heat_capacity_v(v, t), atol=1e-7)


def test_invalid_inputs():
    eos = model()
    for v, t in [(-1, 70), (2.3, 0), (2.3, np.nan)]:
        with pytest.raises(EosValidationError):
            eos.pressure(v, t)
    with pytest.raises(EosValidationError):
        DebyeAnharmonicHelmholtz(eos.rt_eos, 70, 86.44, 2.68, 0.0024, 0.0128, -1, 7.85)


def test_eosmat_round_trip():
    document = get_material_document("argon_fcc")
    ids = [
        record["identifier"]
        for record in document["eos_records"]
        if record.get("thermal", {}).get("type") == "DebyeAnharmonicHelmholtz"
    ]
    assert len(ids) == 1
    material = Material.from_eosmat(document, record_identifiers=ids)
    reloaded = Material.from_eosmat(material.to_eosmat())
    for item in [material, reloaded]:
        record = item.eos_records[0]
        cell_volume = 23.0 * 4.0 / 0.602214076
        pressure = record.pressure(cell_volume, 70)
        assert pressure == pytest.approx(0.0783517375843629, abs=1e-12)
        assert record.volume(pressure, 70) == pytest.approx(cell_volume, abs=1e-8)
        assert item.eos_records[0].eos.pressure(2.3, 70) == pytest.approx(
            0.0783517375843629
        )


def test_caloric_limits_and_zero_anharmonic_reduction():
    from peritheos.eos.thermal import MieGruneisenDebye

    eos = model()
    volume = eos.rt_eos.V0
    gas_constant = 8.31451
    # Third-law coefficient includes both Debye and the published T^4 term.
    low_temperature = 0.01
    cubic_cv = (
        12 * np.pi**4 * gas_constant / (5 * eos.theta0**3)
        - 12 * eos.b1 * gas_constant / eos.theta0**3
    )
    assert eos.molar_heat_capacity_v(
        volume, low_temperature
    ) / low_temperature**3 == pytest.approx(cubic_cv, rel=1e-8)
    assert eos.entropy(volume, low_temperature) / low_temperature**3 == pytest.approx(
        cubic_cv / 3, rel=1e-8
    )
    # Mathematical high-T asymptote, outside the empirical validity envelope.
    # The anharmonic Cv grows linearly negative; it must not be silently clipped.
    high_temperature = 1e10
    slope = -2 * eos.b1 * gas_constant / (eos.b2 * eos.theta0)
    assert eos.molar_heat_capacity_v(
        volume, high_temperature
    ) / high_temperature == pytest.approx(slope, rel=1e-6)
    harmonic = DebyeAnharmonicHelmholtz(
        eos.rt_eos, 300, eos.theta0, eos.gamma0, 0, 0, eos.b2, eos.b3
    )
    debye = MieGruneisenDebye(
        eos.rt_eos,
        300,
        eos.theta0,
        eos.gamma0,
        0,
        1,
        Cvmax=3 * gas_constant,
        thermal_pressure_reference="absolute_zero",
    )
    v = np.array([1.6, 2.0, 2.3])
    t = np.array([600.0, 200.0, 70.0])
    assert harmonic.pressure(v, t) == pytest.approx(debye.pressure(v, t), rel=1e-12)
    assert harmonic.molar_heat_capacity_v(v, t) == pytest.approx(
        debye.molar_heat_capacity_v(v, t), rel=1e-8
    )
    rebased = DebyeAnharmonicHelmholtz(
        eos.rt_eos, 300, eos.theta0, eos.gamma0, eos.q, eos.b1, eos.b2, eos.b3
    )
    assert rebased.pressure(v, t) == pytest.approx(eos.pressure(v, t), rel=1e-12)


def test_schema_model_pair_is_explicit():
    import copy

    from jsonschema import Draft202012Validator, ValidationError

    from peritheos import eosmat_schema

    document = get_material_document("argon_fcc")
    validator = Draft202012Validator(eosmat_schema())
    validator.validate(document)
    changed = copy.deepcopy(document)
    record = next(
        record
        for record in changed["eos_records"]
        if record.get("thermal", {}).get("type") == "DebyeAnharmonicHelmholtz"
    )
    record["thermal"]["model"] = "mie_gruneisen_debye"
    with pytest.raises(ValidationError):
        validator.validate(changed)
