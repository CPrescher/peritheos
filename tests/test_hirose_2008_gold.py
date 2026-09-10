"""Primary-source reproduction, interchange, and cross-calibration for Hirose Au."""

import copy
import csv
import hashlib

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    eosmat_schema,
    find_pressure_calibration_path,
    get_cross_calibration_edge,
    get_eos_record,
    get_material,
    get_material_document,
    recalculate_xrd_pressure_scale,
)
from peritheos.eos.rt import BM2, BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.fitting import fit_thermal_eos
from scripts.reproduce_hirose_2008_gold import DATA, reproduce, source_pressure

PREFIX = "gold_hirose_2008_bm3_"


def test_complete_paired_primary_table_and_uncertainties():
    dataset = get_material("gold").get_dataset("gold_hirose_2008_table1")
    document = get_material_document("gold")
    raw = next(d for d in document["datasets"] if d["identifier"] == dataset.identifier)
    assert hashlib.sha256(DATA.read_bytes()).hexdigest() == raw["resource"]["sha256"]
    with DATA.open() as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 21
    assert [sum(r["run"] == str(i) for r in rows) for i in (1, 2, 3)] == [9, 4, 8]
    assert float(rows[0]["mgo_a_angstrom"]) == 4.12964
    assert float(rows[8]["pressure_gpa"]) == 140.58
    assert float(rows[14]["temperature_k"]) == 2330
    assert float(rows[19]["gold_fit2_pressure_gpa"]) == 119.7
    assert rows[-1]["gold_a_angstrom_uncertainty"] == ""
    anchor = get_eos_record("mgo_speziale_2001_bm3_2")
    assert (
        max(
            abs(
                anchor.pressure(float(row["mgo_volume_a3"]))
                - float(row["pressure_gpa"])
            )
            for row in rows
            if row["run"] == "1"
        )
        < 0.04
    )  # Inside the smallest printed MgO pressure error (0.08 GPa).
    pv = dataset.as_pressure_volume()
    assert pv.volume[0] == pytest.approx(4.00755**3)
    assert pv.volume_uncertainty[0] == pytest.approx(3 * 4.00755**2 * 0.00046)
    assert pv.pressure_uncertainty[0] == 0.08
    assert pv.volume_sigma is None  # Confidence convention is not published.


def test_table1_predictions_and_diagnostic_refits():
    result = reproduce()
    assert result["fit2"]["table1_output_max_difference_gpa"] < 0.06
    assert result["fit2"]["table1_2070k_pressure_gpa"] == pytest.approx(119.7, abs=0.05)
    assert result["300k"]["refit_parameters"][0] == pytest.approx(5.58, abs=0.02)
    assert result["fit1"]["refit_parameters"] == pytest.approx(
        [-0.028, 3.179e-5, 1.477e-8], rel=0.02
    )
    for metrics in result.values():
        assert metrics["native_equation_max_difference_gpa"] < 1e-10
    assert not result["fit2"]["full_fit_reproduced"]


@pytest.mark.parametrize("suffix", ["fit1", "fit2"])
def test_thermal_native_python_derivatives_and_roundtrips(suffix):
    record = get_eos_record(PREFIX + suffix)
    eos = record.eos
    volumes = np.array([67.85, 58.0, 50.0])[:, None]
    temperatures = np.array([300.0, 1500.0, 2330.0])[None, :]
    expected = source_pressure(volumes, temperatures, fit=1 if suffix == "fit1" else 2)
    assert eos.pressure(volumes, temperatures) == pytest.approx(expected)
    assert eos.thermal_pressure(volumes, 300.0) == pytest.approx(
        np.zeros((3, 1)), abs=1e-12
    )
    fallback = copy.copy(eos)
    del fallback._native
    assert fallback.pressure(volumes, temperatures) == pytest.approx(expected)
    h = 1e-5
    derivative = (
        -volumes
        * (
            eos.pressure(volumes + h, temperatures)
            - eos.pressure(volumes - h, temperatures)
        )
        / (2 * h)
    )
    assert eos.bulk_modulus(volumes, temperatures) == pytest.approx(
        derivative, rel=1e-8
    )
    assert fallback.bulk_modulus(volumes, temperatures) == pytest.approx(
        eos.bulk_modulus(volumes, temperatures)
    )
    pressure = record.pressure([55.0, 51.0], temperature=[1600.0, 2200.0])
    assert record.volume(pressure, temperature=[1600.0, 2200.0]) == pytest.approx(
        [55.0, 51.0]
    )
    assert eos.bulk_modulus(67.85, 300.0) == pytest.approx(167.0)


def test_schema_roundtrip_and_explicit_calibration_scope():
    document = get_material_document("gold")
    Draft202012Validator(eosmat_schema()).validate(document)
    selected = Material.from_eosmat(document, record_identifiers=[PREFIX + "fit2"])
    restored = Material.from_eosmat(selected.to_eosmat())
    assert restored.eos_records[0].pressure(51.0, temperature=2100.0) == pytest.approx(
        selected.eos_records[0].pressure(51.0, temperature=2100.0)
    )
    edge = get_cross_calibration_edge("gold_hirose_2008_to_mgo_speziale_2001_thermal")
    assert edge["executable"] is False
    assert find_pressure_calibration_path(PREFIX + "300k", "mgo_speziale_2001_bm3_2")
    pressure = recalculate_xrd_pressure_scale(
        100.0, PREFIX + "fit2", "gold_fei_2007_vinet_2", temperature_k=2000.0
    )
    assert recalculate_xrd_pressure_scale(
        pressure.target_pressure_gpa,
        "gold_fei_2007_vinet_2",
        PREFIX + "fit2",
        temperature_k=2000.0,
    ).target_pressure_gpa == pytest.approx(100.0)


@pytest.mark.parametrize(
    "updates",
    [
        {"bulk_modulus_law": "unknown"},
        {"beta1": 1e-6},
        {"bulk_modulus_law": "reciprocal_cubic", "dK_dT": -0.01},
        {"kprime_log_coefficient": float("nan")},
    ],
)
def test_temperature_law_invalid_configuration(updates):
    parameters = dict(rt_eos=BM3(67.85, 167.0, 5.58), Tr=300.0, alpha0=0.0, dK_dT=0.0)
    parameters.update(updates)
    with pytest.raises(ValueError):
        ThermalReferenceStateEOS(**parameters)
    with pytest.raises(ValueError):
        ThermalReferenceStateEOS(
            BM2(67.85, 167.0), 300.0, 0.0, 0.0, kprime_log_coefficient=0.001
        )


def test_nonphysical_compressibility_and_invalid_states():
    eos = get_eos_record(PREFIX + "fit2").eos.with_parameters(beta1=-0.01)
    with pytest.raises(ValueError):
        eos.pressure(50.0, 400.0)
    for volume, temperature in (
        (0.0, 300.0),
        (50.0, 0.0),
        (np.nan, 300.0),
        (50.0, np.inf),
    ):
        with pytest.raises(ValueError):
            eos.pressure(volume, temperature)


def test_native_fit_retains_cubic_and_logarithmic_temperature_laws():
    eos = get_eos_record(PREFIX + "fit2").eos
    volumes = np.linspace(49.0, 60.0, 20)
    temperatures = np.linspace(1000.0, 2300.0, 20)
    pressure = eos.pressure(volumes, temperatures)
    fit = fit_thermal_eos(
        ThermalReferenceStateEOS,
        eos.rt_eos,
        volume=volumes,
        temperature=temperatures,
        pressure=pressure,
        initial={"kprime_log_coefficient": 3e-4},
        configuration=eos.configuration_values(),
        fixed={
            name: value
            for name, value in eos.parameter_values(include_reference=False).items()
            if name != "kprime_log_coefficient"
        },
    )
    assert fit.parameters["kprime_log_coefficient"] == pytest.approx(3.61e-4, rel=1e-6)
