"""Tests of literal equations, not claims that Table 8 has been reproduced."""

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos.eos.experimental.maltby_2024 import Maltby2024Published, fcc_shells

ROOT = Path(__file__).resolve().parents[1]


def test_geometric_coordination_shells_and_confirmed_source_entries():
    squared, population = fcc_shells(64)
    shells = dict(zip(squared, population))
    assert [shells[m] for m in [1, 2, 3, 4, 5, 6, 7, 8]] == [
        12,
        6,
        24,
        12,
        24,
        8,
        48,
        6,
    ]
    assert shells[49] == 108
    assert shells[64] == 12
    # SI.1 incorrectly lists 48 at m=14. No integer three-square sum equals
    # 2*m=28, so that entry is not part of the geometric fcc lattice.
    assert 14 not in shells
    assert 30 not in shells and 46 not in shells and 56 not in shells


@pytest.mark.parametrize("volume", [1.2, 1.8, 2.397])
@pytest.mark.parametrize("temperature", [0.0, 5.0, 70.0, 300.0])
def test_helmholtz_pressure_identity_and_stable_inverse(volume, temperature):
    model = Maltby2024Published(64)
    h = volume * 1e-5
    numerical_pressure = (
        -(
            model.molar_helmholtz_energy(volume + h, temperature)
            - model.molar_helmholtz_energy(volume - h, temperature)
        )
        / (2 * h)
        * 1e-4
    )
    assert model.pressure(volume, temperature) == pytest.approx(
        numerical_pressure, abs=1e-7
    )
    assert model.volume(
        model.pressure(volume, temperature), temperature
    ) == pytest.approx(volume, abs=1e-10)


def test_independent_si_quadrature_and_retained_literature_mismatch():
    report = json.loads(
        (ROOT / "docs/data/argon-maltby-2024-reproduction.json").read_text()
    )
    model = Maltby2024Published(64)
    for row in report["independent_si_quadrature"]:
        actual = model.pressure(row["volume_cm3_mol"] / 10, row["temperature_k"])
        assert actual == pytest.approx(row["pressure_gpa"], abs=1e-7)
    assert model.volume(0.001, 70) * 10 == pytest.approx(23.8904275342, abs=1e-8)
    assert report["status"] == "not_reproduced"
    assert abs(model.volume(0.001, 70) * 10 - 23.97) > 0.07
    assert model.thermal_pressure_increment(2.0, 300) == 0
    # The published implementation is available without a validation claim.
    records = json.loads(
        (ROOT / "peritheos/data/materials/argon_fcc.eosmat").read_text()
    )["eos_records"]
    record = next(row for row in records if "maltby" in row["identifier"])
    assert record["scientific_validation"]["status"] == "not_reproduced"
    assert record["default"] is False


@pytest.mark.parametrize("cutoff", [0, 30, 257, 1.5, True])
def test_invalid_cutoff(cutoff):
    with pytest.raises(ValueError):
        Maltby2024Published(cutoff)


@pytest.mark.parametrize("v,t", [(0, 70), (-1, 70), (np.nan, 70), (2, -1), (2, np.inf)])
def test_invalid_states(v, t):
    with pytest.raises(ValueError):
        Maltby2024Published(64).pressure(v, t)


def test_invalid_derivative_and_inverse():
    model = Maltby2024Published(64)
    with pytest.raises(ValueError):
        model.bulk_modulus(2, 70, 0)
    with pytest.raises(ValueError):
        model.volume(np.nan, 70)
    with pytest.raises(ValueError):
        model.volume(-100, 70)


def test_unvalidated_catalog_execution_and_export_preserve_qualification():
    from peritheos.catalog import get_eos_record, get_material
    from peritheos.eosmat import get_material_document, validate_eosmat_document
    from peritheos.errors import MaterialError
    from peritheos.materials import Material

    identifier = "argon_fcc_maltby_2024_published"
    record = get_eos_record(identifier)
    assert record.is_thermal and not record.is_default
    assert record.scientific_validation_status == "not_reproduced"
    volumes = np.array([1.2, 1.8, 2.397]) / record.volume_scale
    temperatures = np.array([0, 300, 70])
    pressures = record.pressure(volumes, temperatures)
    assert pressures == pytest.approx(
        [16.50973820028, 2.21438377236, -0.00422597965], abs=1e-7
    )
    assert record.volume(pressures, temperatures) == pytest.approx(volumes, abs=1e-8)
    assert record.thermal_pressure_increment(volumes, 300) == pytest.approx([0, 0, 0])
    material = get_material("argon_fcc")
    exported = material.to_eosmat()
    validate_eosmat_document(exported)
    row = next(
        row for row in exported["eos_records"] if row["identifier"] == identifier
    )
    assert row["scientific_validation"]["status"] == "not_reproduced"
    assert "Unvalidated" in row["reproduction"]["display_label"]
    assert row["eos"]["parameters"] == {"shell_cutoff_squared": 64}
    assert "thermal" not in row
    assert "V0" not in row["eos"]["parameters"]
    assert material.default_record().identifier == "argon_fcc_dewaele_2021_vinet_mgd"
    with pytest.raises(MaterialError, match="not_reproduced"):
        Material.from_eosmat(
            get_material_document("argon_fcc"), record_identifiers=[identifier]
        )
    restored = Material.from_eosmat(exported, require_primary_validation=False)
    assert restored.eos_records[-1].pressure(volumes, temperatures) == pytest.approx(
        pressures
    )


def test_legacy_snapshot_does_not_upgrade_unvalidated_record():
    from peritheos.catalog import get_material
    from peritheos.materials import Material

    material = get_material("argon_fcc")
    snapshot = material.to_snapshot_dict()
    row = snapshot["eos_records"][-1]
    assert row["scientific_validation"]["status"] == "not_reproduced"
    assert "Unvalidated" in row["scientific_validation"]["note"]
    with pytest.raises(ValueError, match="argon_fcc_maltby_2024_published"):
        Material.from_dict(snapshot)
