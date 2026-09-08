import json
from pathlib import Path

import numpy as np
import pytest

from peritheos.eos.rt.vinet import Vinet
from peritheos.eos.thermal import DorogokupetsOganov2007
from peritheos.eosmat import get_material_document
from peritheos.materials import Material
from scripts.audit_dorogokupets_oganov_2007_platinum import build_audit

ROOT = Path(__file__).resolve().parents[1]


def _platinum_eos(model_class=DorogokupetsOganov2007):
    return model_class(
        Vinet(0.9091, 276.07, 5.30),
        Tr=298.15,
        theta_B1=95.2,
        d_B1=8.199,
        m_B1=0.329,
        theta_B2=148.4,
        d_B2=4.005,
        m_B2=0.383,
        theta_E1=214.6,
        m_E1=1.211,
        theta_E2=140.8,
        m_E2=1.077,
        gamma0=2.802,
        gamma_inf=1.538,
        beta=5.550,
        anharmonic_a=160.9,
        anharmonic_m=4.06,
        electronic_e=260.0,
        electronic_g=2.4,
        defect_H=32572.0,
        defect_S=0.631,
        n=1.0,
    )


class _PythonDorogokupetsOganov2007(DorogokupetsOganov2007):
    """Force the supported subclass/fallback path for compatibility testing."""


def test_python_fallback_matches_native_pt_evaluator():
    native = _platinum_eos()
    fallback = _platinum_eos(_PythonDorogokupetsOganov2007)
    volumes = 0.9091 * np.array([1.0, 0.9, 0.8, 0.7])
    temperatures = np.array([298.15, 1000.0, 2000.0, 3000.0])

    assert hasattr(native, "_native")
    assert not hasattr(fallback, "_native")
    assert fallback.pressure(volumes, temperatures) == pytest.approx(
        native.pressure(volumes, temperatures), rel=5.0e-12, abs=5.0e-12
    )
    for quantity in (
        "absolute_thermal_pressure",
        "thermal_helmholtz_free_energy",
        "thermal_internal_energy",
        "thermal_entropy",
        "molar_heat_capacity_v",
        "molar_heat_capacity_p",
        "adiabatic_bulk_modulus",
    ):
        assert getattr(fallback, quantity)(volumes, temperatures) == pytest.approx(
            getattr(native, quantity)(volumes, temperatures),
            rel=5.0e-11,
            abs=5.0e-11,
        )


def test_dorogokupets_oganov_helmholtz_derivatives_are_self_consistent():
    eos = _platinum_eos(_PythonDorogokupetsOganov2007)
    for volume, temperature in ((0.9091, 298.15), (0.82, 1000.0), (0.7, 3000.0)):
        volume_step = volume * 1.0e-6
        pressure_from_free_energy = (
            -(
                eos.thermal_helmholtz_free_energy(volume + volume_step, temperature)
                - eos.thermal_helmholtz_free_energy(volume - volume_step, temperature)
            )
            / (2.0 * volume_step)
            / 1.0e4
        )
        assert pressure_from_free_energy == pytest.approx(
            eos.absolute_thermal_pressure(volume, temperature), rel=2.0e-9
        )

        temperature_step = temperature * 1.0e-5
        heat_capacity_from_energy = (
            eos.thermal_internal_energy(volume, temperature + temperature_step)
            - eos.thermal_internal_energy(volume, temperature - temperature_step)
        ) / (2.0 * temperature_step)
        assert heat_capacity_from_energy == pytest.approx(
            eos.molar_heat_capacity_v(volume, temperature), rel=2.0e-9
        )


def test_dorogokupets_oganov_pt_reproduces_published_standard_entropy():
    eos = _platinum_eos()
    # Text below Table I: calculated S_298 = 41.45 J mol^-1 K^-1.
    assert eos.thermal_entropy(eos.rt_eos.V0, eos.Tr) == pytest.approx(41.45, abs=0.005)


def test_dorogokupets_oganov_pt_reproduces_table_vi_isochors():
    eos = _platinum_eos()
    temperatures = (298.15, 1000.0, 2000.0, 3000.0)
    rows = {
        1.00: (0.0, 5.309, 12.864, 20.349),
        0.95: (16.207, 21.219, 28.485, 35.822),
        0.90: (38.307, 43.111, 50.199, 57.476),
        0.85: (68.389, 73.069, 80.079, 87.376),
        0.80: (109.388, 114.016, 121.043, 128.434),
        0.75: (165.473, 170.119, 177.253, 184.808),
        0.70: (242.676, 247.403, 254.730, 262.523),
    }
    for ratio, expected_pressures in rows.items():
        for temperature, expected_pressure in zip(temperatures, expected_pressures):
            assert eos.pressure(ratio * eos.rt_eos.V0, temperature) == pytest.approx(
                expected_pressure, abs=0.02
            )


def test_bundled_dorogokupets_oganov_pt_record_uses_cell_volume_units():
    material = Material.from_eosmat(
        get_material_document("platinum"),
        record_identifiers=["platinum_dorogokupets_oganov_2007_vinet_4"],
    )
    record = material.eos_records[0]
    assert record.reference_volume == pytest.approx(60.38384263870976)
    assert record.pressure(0.8 * record.reference_volume, 3000.0) == pytest.approx(
        128.434, abs=0.02
    )


def test_dorogokupets_oganov_pt_audit_artifact_is_current_and_qualified():
    expected = build_audit()
    stored = json.loads(
        (
            ROOT / "docs" / "data" / "dorogokupets-oganov-2007-platinum-audit.json"
        ).read_text(encoding="utf-8")
    )
    assert stored == expected
    assert stored["global_objective"]["status"] == "not_exactly_reconstructible"
    assert stored["exact_reconstructions"]["table_vi"]["rows"] == 28
    assert stored["exact_reconstructions"]["table_vi"]["max_abs_gpa"] < 0.02
    assert stored["partial_validations"]["static_diffraction"]["rows"] == 36
    assert stored["partial_validations"]["shock"]["rows"] == 7
