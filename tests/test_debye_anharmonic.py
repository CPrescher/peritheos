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
