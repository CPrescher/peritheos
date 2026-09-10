import csv
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import Avogadro

from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import Dewaele2006
from peritheos.eosmat import get_material_document
from peritheos.materials import Material

ROOT = Path(__file__).resolve().parents[1]
RECORD_ID = "iron_dewaele_2006_vinet_thermal"


class _PythonDewaele2006(Dewaele2006):
    """Force the Python fallback for native/fallback parity testing."""


def _model(model_class=Dewaele2006):
    molar_atomic_v0 = 11.214 * Avogadro * 1.0e-25
    return model_class(
        Vinet(molar_atomic_v0, 163.4, 5.38),
        Tr=300.0,
        theta0=417.0,
        gamma0=1.875,
        gamma_inf=1.305,
        beta=1.875 / (1.875 - 1.305),
        anharmonic_a=3.7e-5,
        anharmonic_m=1.87,
        electronic_e=1.95e-4,
        electronic_g=1.339,
        n=1.0,
    )


def test_python_fallback_matches_native_dewaele_2006_evaluator():
    native = _model()
    fallback = _model(_PythonDewaele2006)
    volumes = native.rt_eos.V0 * np.array([1.0, 0.9, 0.8, 0.7])
    temperatures = np.array([300.0, 1500.0, 3000.0, 6000.0])

    assert hasattr(native, "_native")
    assert not hasattr(fallback, "_native")
    assert fallback.pressure(volumes, temperatures) == pytest.approx(
        native.pressure(volumes, temperatures), rel=1.0e-11, abs=5.0e-12
    )
    assert fallback.characteristic_temperature(volumes) == pytest.approx(
        native.characteristic_temperature(volumes), rel=5.0e-12
    )


def test_dewaele_2006_terms_are_reference_subtracted_and_additive():
    model = _model()
    volumes = model.rt_eos.V0 * np.array([1.0, 0.85, 0.7])

    assert model.thermal_pressure(volumes, 300.0) == pytest.approx(0.0, abs=1.0e-14)
    total = model.thermal_pressure(volumes, 6000.0)
    terms = (
        model.vibrational_pressure_increment(volumes, 6000.0)
        + model.anharmonic_pressure_increment(volumes, 6000.0)
        + model.electronic_pressure_increment(volumes, 6000.0)
    )
    assert total == pytest.approx(terms, rel=5.0e-13)


def test_bundled_dewaele_2006_record_uses_hcp_cell_volume_basis():
    record = Material.from_eosmat(
        get_material_document("iron"), record_identifiers=[RECORD_ID]
    ).eos_records[0]

    assert record.reference_volume == pytest.approx(22.428)
    assert record.eos.rt_eos.V0 == pytest.approx(11.214 * Avogadro * 1.0e-25)
    assert record.pressure(record.reference_volume, 300.0) == pytest.approx(
        0.0, abs=1.0e-14
    )


def test_dewaele_2006_recalculates_fischer_hcp_fe_pressures():
    record = Material.from_eosmat(
        get_material_document("iron"), record_identifiers=[RECORD_ID]
    ).eos_records[0]
    dataset = ROOT / "peritheos/data/datasets/feo-fischer-2011-table-s1-pvt.csv"
    residuals = []
    with dataset.open(newline="", encoding="utf-8") as stream:
        for row in csv.DictReader(stream):
            if not row["iron_molar_volume_cm3_mol"]:
                continue
            atomic_molar_volume = float(row["iron_molar_volume_cm3_mol"])
            cell_volume = atomic_molar_volume * 1.0e24 / Avogadro * 2.0
            calculated = record.pressure(cell_volume, float(row["temperature_k"]))
            residuals.append(calculated - float(row["pressure_gpa"]))

    assert len(residuals) == 52
    assert max(abs(residual) for residual in residuals) < 0.023
