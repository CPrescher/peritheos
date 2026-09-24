"""Exact historical Murnaghan/zero-point equation, source data and interchange."""

import csv
import hashlib
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, get_material_document
from peritheos.eos.finger1981 import Finger1981Argon
from peritheos.eos.rt import Murnaghan, SecondOrderMurnaghan
from peritheos.eos.thermal import Dewaele2006
from peritheos.errors import EosValidationError
from peritheos.fitting import fit_thermal_eos
from scripts.reproduce_argon_finger1981 import ID, REPORT, ROOT, reproduce


def record():
    return Material.from_eosmat(get_material_document("argon_fcc")).get_eos_record(ID)


def test_source_transcription_and_reproduction(assert_audit_close):
    doc = get_material_document("argon_fcc")
    assert not list(
        Draft202012Validator(
            json.loads((ROOT / "peritheos/data/eosmat-v3.schema.json").read_text())
        ).iter_errors(doc)
    )
    ds = next(
        d for d in doc["datasets"] if d["identifier"] == "argon_fcc_finger_1981_table1"
    )
    path = ROOT / "peritheos/data" / ds["resource"]["path"]
    assert hashlib.sha256(path.read_bytes()).hexdigest() == ds["resource"]["sha256"]
    with path.open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 19
    for row in rows:
        assert float(row["volume_a3"]) == pytest.approx(
            float(row["a_angstrom"]) ** 3, abs=1e-11
        )
        assert float(row["volume_uncertainty_a3"]) == pytest.approx(
            3 * float(row["a_angstrom"]) ** 2 * float(row["a_uncertainty_angstrom"]),
            abs=1e-11,
        )
    audit = reproduce()
    assert_audit_close(audit, json.loads(REPORT.read_text()))
    assert audit["parameter_reproduction"]["both_within_reported_uncertainties"]
    assert (
        audit["fit_reproduction_status"]
        == "parameters_reproduced_within_reported_uncertainties"
    )
    assert audit["native_max_abs_difference_gpa"] < 1e-9
    assert audit["published_molar_volume_rmse_gpa"] == pytest.approx(
        0.0953552, abs=1e-7
    )
    assert audit["static_offset_plus_zero_point_at_v0_gpa"] == pytest.approx(
        0.00157018, abs=1e-8
    )


def test_murnaghan_equation_and_modulus():
    model = SecondOrderMurnaghan(22.557, 2.3701, 6.97, -0.4, -0.10289)
    v = np.array([14.0, 17.0, 22.557])
    step = 1e-5
    assert model.pressure(model.V0) == pytest.approx(-0.10289)
    assert model.bulk_modulus(v) == pytest.approx(
        -v * (model.pressure(v + step) - model.pressure(v - step)) / (2 * step),
        rel=1e-8,
    )
    assert model.volume(model.pressure(v)) == pytest.approx(v, rel=1e-9)
    # K''=0 reduces to first-order Murnaghan, not Birch-Murnaghan.
    first = SecondOrderMurnaghan(22.557, 2.3701, 6.97, 0, 0)
    assert first.pressure(v) == pytest.approx(
        Murnaghan(22.557, 2.3701, 6.97).pressure(v)
    )
    with pytest.raises(EosValidationError):
        SecondOrderMurnaghan(1, 1, 1, 2)


def test_native_python_fallback_and_interchange():
    r = record()
    v = np.array([100.0, 115.0, 130.0])
    t = np.array([77.0, 293.0, 500.0])
    expected = Finger1981Argon().pressure(v, t)
    assert r.pressure(v, t) == pytest.approx(expected, abs=1e-9)
    assert r.volume(expected, t) == pytest.approx(v, rel=1e-9)

    class Fallback(Dewaele2006):
        pass

    fallback = Fallback(
        r.eos.rt_eos,
        **r.eos.parameter_values(include_reference=False),
        **r.eos.configuration_values(),
    )
    assert fallback.pressure(v * r.volume_scale, t) == pytest.approx(expected, abs=1e-9)
    assert fallback.thermal_pressure_increment(
        v * r.volume_scale, 293
    ) == pytest.approx(0, abs=1e-12)
    assert r.eos.thermal_pressure_increment(v * r.volume_scale, 293) == pytest.approx(
        0, abs=1e-12
    )
    restored = Material.from_eosmat(
        Material.from_eosmat(get_material_document("argon_fcc")).to_eosmat()
    ).get_eos_record(ID)
    assert restored.eos.configuration_values()["zero_point_pressure"] == "included"
    assert restored.pressure(v, t) == pytest.approx(expected, abs=1e-9)
    with pytest.raises(EosValidationError):
        Dewaele2006(
            r.eos.rt_eos,
            **r.eos.parameter_values(include_reference=False),
            zero_point_pressure="included",
        )


def test_native_fit_preserves_absolute_and_zero_point_configuration():
    model = record().eos
    v = np.linspace(1.45, 1.95, 12)
    t = np.linspace(100, 500, 12)
    params = model.parameter_values(include_reference=False)
    params.pop("theta0")
    result = fit_thermal_eos(
        Dewaele2006,
        model.rt_eos,
        v,
        t,
        model.pressure(v, t),
        initial={"theta0": 95},
        fixed=params,
        configuration=model.configuration_values(),
    )
    assert result.model.theta0 == pytest.approx(93.3, abs=1e-4)
    assert result.model.pressure(v, t) == pytest.approx(model.pressure(v, t), abs=1e-7)


def test_zero_point_pressure_is_static_and_not_rebased():
    from scipy.constants import R

    eos = record().eos
    parameters = eos.parameter_values(include_reference=False)
    no_zero = Dewaele2006(
        eos.rt_eos, **parameters, thermal_pressure_reference="absolute_zero"
    )
    v = np.array([1.5, 1.9, 2.2557])
    zp = (
        9
        / 8
        * eos.n
        * R
        * eos.gruneisen_parameter(v)
        * eos.characteristic_temperature(v)
        / v
        / 1e4
    )
    for temperature in [0.001, 77.0, 293.0, 500.0]:
        assert eos.pressure(v, temperature) - no_zero.pressure(
            v, temperature
        ) == pytest.approx(zp, abs=1e-12)
    assert eos.pressure(v, 0.001) == pytest.approx(
        eos.rt_eos.pressure(v) + zp, abs=1e-12
    )
    # Adding zero-point pressure does not add temperature-dependent heat capacity.
    assert eos.thermal_pressure_increment(v, 500.0) == pytest.approx(
        no_zero.thermal_pressure_increment(v, 500.0), abs=1e-12
    )
    parameters["Tr"] = 77.0
    rebased = Dewaele2006(
        eos.rt_eos,
        **parameters,
        thermal_pressure_reference="absolute_zero",
        zero_point_pressure="included",
    )
    assert rebased.pressure(v, 293.0) == pytest.approx(
        eos.pressure(v, 293.0), abs=1e-12
    )


def test_schema_zero_point_requires_static_lattice_baseline():
    import copy

    from jsonschema import Draft202012Validator, ValidationError

    from peritheos import eosmat_schema

    original = get_material_document("argon_fcc")
    validator = Draft202012Validator(eosmat_schema())
    validator.validate(original)
    for configuration in [False, True]:
        document = copy.deepcopy(original)
        thermal = next(
            r["thermal"] for r in document["eos_records"] if r["identifier"] == ID
        )
        if configuration:
            thermal["configuration"] = {
                key: thermal.pop(key)
                for key in ["zero_point_pressure", "thermal_pressure_reference"]
            }
            target = thermal["configuration"]
        else:
            target = thermal
        validator.validate(document)
        target["thermal_pressure_reference"] = "reference_temperature"
        with pytest.raises(ValidationError):
            validator.validate(document)


def test_second_order_murnaghan_inverse_respects_regular_branch():
    # This generic positive-K'' case has a finite-volume compression pole.
    eos = SecondOrderMurnaghan(10.0, 2.0, 4.0, 1.0, 0.2)
    volumes = np.array([5.0, 7.0, 9.0, 10.0, 12.0])
    assert eos.volume(eos.pressure(volumes)) == pytest.approx(volumes, rel=1e-12)
    with pytest.raises(EosValidationError):
        eos.volume(-1.0)
    # Finger's negative-K'' case has an upper pressure limit on its static curve.
    finger = SecondOrderMurnaghan(22.557, 2.3701, 6.97, -0.4, -0.10289)
    with pytest.raises(EosValidationError):
        finger.volume(100.0)
