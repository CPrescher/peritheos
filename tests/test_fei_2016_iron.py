import hashlib
import json

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import DebyeQuadraticThermalPressure
from peritheos.errors import EosValidationError, UnsupportedOperationError
from peritheos.fitting import fit_thermal_eos
from peritheos.materials import Material
from scripts.reproduce_fei_2016_iron import (
    DATA,
    MASS,
    NA,
    RECORDS,
    ROOT,
    V0,
    bm3,
    csv_path,
    reproduce,
    rows,
    thermal_pressure,
)


def model(**changes):
    params = dict(Tr=300, theta0=422, gamma0=1.74, q=0.78, n=1, A=5.78865e-7, m=0.34)
    params.update(changes)
    return DebyeQuadraticThermalPressure(
        BM3(V0 * NA * 1e-25 / 2, 172.7, 4.79), **params
    )


def test_complete_source_tables_and_missing_errors():
    doc = get_material_document("iron")
    for table, count in [("s1", 74), ("s2", 22)]:
        source = rows(table)
        assert len(source) == count
        assert [int(r["source_row"]) for r in source] == list(range(3, count + 3))
        dataset = next(
            d
            for d in doc["datasets"]
            if d["identifier"] == f"iron_fei_2016_table_{table}"
        )
        assert (
            hashlib.sha256(csv_path(table).read_bytes()).hexdigest()
            == dataset["resource"]["sha256"]
        )
        assert list(source[0]) == [c["name"] for c in dataset["columns"]]
        # Independently check the conventional hcp cell against printed a,c.
        for row in source:
            v = (
                np.sqrt(3)
                / 2
                * float(row["a_angstrom"]) ** 2
                * float(row["c_angstrom"])
            )
            assert v == pytest.approx(
                float(row["volume_a3_conventional_cell"]), abs=2e-4
            )
    s1 = rows("s1")
    assert sum(r["pressure_standard"] == "MgO" for r in s1) == 18
    assert sum(r["marker_a_error_angstrom"] == "" for r in s1) == 18
    assert s1[0]["sample_name"] == "c112_009"
    assert float(s1[-1]["pressure_gpa"]) == 204.6
    assert rows("s2")[-1]["sample_name"] == "c118Pt-Fe_292"
    assert float(rows("s2")[0]["temperature_error_k"]) == 18
    assert (DATA / "iron-fei-2016-source.json").exists()


def test_published_parameters_source_states_and_interchange():
    doc = get_material_document("iron")
    material = Material.from_eosmat(doc)
    selected = {
        r["identifier"]: r for r in doc["eos_records"] if r["identifier"] in RECORDS
    }
    assert len(selected) == 3
    assert selected[RECORDS[0]]["eos"]["parameters"] == dict(
        V0=V0, K0=172.7, K0_prime=4.79
    )
    free = selected[RECORDS[1]]
    assert free["eos"]["parameters"]["V0"] == pytest.approx(
        2 * MASS / (NA * 1e-24 * 8.3602)
    )
    assert all(
        r["record_kind"] == "published" and "default_for" not in r
        for r in selected.values()
    )
    record = material.get_eos_record(RECORDS[2])
    volumes = np.array([18.764694, 17.1766, 14.0])
    temperatures = np.array([1208.0, 1795.0, 5000.0])
    expected = [bm3(v) + thermal_pressure(v, t) for v, t in zip(volumes, temperatures)]
    assert record.pressure(volumes, temperatures) == pytest.approx(expected, rel=1e-10)
    assert record.volume(expected, temperatures) == pytest.approx(volumes, rel=1e-8)
    # Source S2 endpoints: tolerate observed residuals, not invented measurement errors.
    assert record.pressure(volumes[0], temperatures[0]) == pytest.approx(53.96, abs=0.1)
    assert record.pressure(volumes[1], temperatures[1]) == pytest.approx(99.69, abs=0.6)
    exported = Material.from_eosmat(material.to_eosmat())
    assert exported.get_eos_record(RECORDS[2]).pressure(
        volumes, temperatures
    ) == pytest.approx(expected)
    assert exported.to_eosmat()["datasets"] == doc["datasets"]
    with pytest.raises(ValueError):
        record.pressure(14, 6000, check_validity=True)


def test_empirical_model_native_python_derivatives_and_limits():
    m = model()
    volume = m.rt_eos.V0 * np.array([[1.0], [0.8], [0.6]])
    temperature = np.array([300.0, 1500.0, 5000.0])
    native = m.pressure(volume, temperature)
    assert m.thermal_pressure(volume, 300) == pytest.approx(np.zeros((3, 1)))
    assert m.pressure(m.rt_eos.V0, 300) == pytest.approx(0, abs=1e-12)
    h = volume * 1e-5
    derivative = (
        -volume
        * (m.pressure(volume + h, temperature) - m.pressure(volume - h, temperature))
        / (2 * h)
    )
    assert m.bulk_modulus(volume, temperature) == pytest.approx(derivative, rel=1e-7)
    assert m.volume(native, temperature) == pytest.approx(
        np.broadcast_to(volume, (3, 3)), rel=1e-8
    )
    assert m.temperature(native, volume) == pytest.approx(
        np.broadcast_to(temperature, (3, 3)), rel=1e-7
    )
    with pytest.raises(UnsupportedOperationError):
        m.molar_heat_capacity_v(volume, temperature)
    del m._native
    assert m.pressure(volume, temperature) == pytest.approx(native, rel=1e-10)
    assert model(A=0).thermal_pressure(volume, temperature) == pytest.approx(
        model()._debye.thermal_pressure(volume, temperature)
    )
    for parameter in ["Tr", "theta0", "n"]:
        with pytest.raises(EosValidationError):
            model(**{parameter: 0})
    for parameter in ["A", "m", "q", "gamma0"]:
        with pytest.raises(EosValidationError):
            model(**{parameter: float("nan")})
    for v, t in [(0, 300), (1, 0), (float("nan"), 300), (1, float("inf"))]:
        with pytest.raises(EosValidationError):
            model().pressure(v, t)


def test_native_fit_recovers_independent_quadratic_pressure():
    source = model()
    v = source.rt_eos.V0 * np.linspace(0.6, 0.95, 12)
    t = np.linspace(1000, 4000, 12)
    cells = v / (NA * 1e-25 / 2)
    p = np.array([bm3(vv) + thermal_pressure(vv, tt) for vv, tt in zip(cells, t)])
    fit = fit_thermal_eos(
        DebyeQuadraticThermalPressure,
        source.rt_eos,
        v,
        t,
        p,
        initial={"A": 5e-7},
        fixed={
            key: value
            for key, value in source.parameter_values(include_reference=False).items()
            if key != "A"
        },
    )
    assert fit.parameters["A"] == pytest.approx(5.78865e-7, rel=1e-6)


def test_partial_refits_and_report_are_reproducible():
    result = reproduce()
    saved = json.loads((ROOT / "docs/data/fei-2016-iron-reproduction.json").read_text())
    assert result["thermal"] == pytest.approx(saved["thermal"])
    assert result["static"]["fixed_density"]["parameters"]["K0"] == pytest.approx(
        174.1450783, rel=1e-7
    )
    assert result["thermal"]["published_pth_rmse_gpa"] == pytest.approx(1.137550586)
    assert result["thermal"]["s2_pth_vs_p_minus_bm3_max_abs_gpa"] == pytest.approx(
        2.065618874
    )
