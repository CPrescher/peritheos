"""Tests of literal equations, not claims that Table 8 has been reproduced."""

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos.eos.experimental.maltby_2024 import Maltby2024Published, fcc_shells

ROOT = Path(__file__).resolve().parents[1]


def test_coordination_shells_match_official_supplement():
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
    # No invented or approximate EOSMAT record has been added.
    records = json.loads(
        (ROOT / "peritheos/data/materials/argon_fcc.eosmat").read_text()
    )["eos_records"]
    assert not any("maltby" in row["identifier"] for row in records)


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
