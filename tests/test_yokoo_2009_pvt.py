"""Au/Pt PVT reconstruction, source holdouts, bounded inversions and interchange."""

import copy
import csv
import json

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    get_eos_record,
    get_material_document,
    search_eos_records,
)
from peritheos.eosmat import validate_eosmat_document
from peritheos.errors import (
    EosmatError,
    EosValidationError,
    MaterialError,
    UnsupportedOperationError,
)
from scripts.fit_yokoo_2009_gold_thermal import check_reconstruction
from scripts.fit_yokoo_2009_platinum_thermal import (
    OUTPUT as PT_OUTPUT,
)
from scripts.fit_yokoo_2009_platinum_thermal import (
    RESIDUALS as PT_RESIDUALS,
)
from scripts.fit_yokoo_2009_platinum_thermal import (
    fit_with_pressure_correction,
)
from scripts.fit_yokoo_2009_platinum_thermal import (
    reconstruct as platinum_report,
)
from scripts.register_yokoo_2009_pvt import (
    OUTPUT,
    RECORDS,
    ROOT,
    reports,
    validate_library,
)


@pytest.mark.parametrize("metal", RECORDS)
def test_complete_source_grid_agrees_with_independent_pressure_and_inverts(metal):
    report = reports()[metal]
    record = get_eos_record(RECORDS[metal])
    rows = report["states"]
    volume = np.array(
        [
            float(r["volume_ratio"]) * report["conventions"]["ambient_cell_a3"]
            for r in rows
        ]
    )
    temperature = np.array([float(r["temperature_k"]) for r in rows])
    expected = np.array([r["model_pressure_gpa"] for r in rows])
    actual = record.pressure(volume, temperature, check_validity=True)
    # Expected values include a fresh optimizer run on each numerical stack.
    assert actual == pytest.approx(expected, abs=5e-8)
    assert record.volume(actual, temperature) == pytest.approx(volume, abs=1e-8)
    assert record.eos.temperature(
        actual, volume * record.volume_scale
    ) == pytest.approx(temperature, abs=1e-6)
    assert record.reference_volume == pytest.approx(
        report["conventions"]["ambient_cell_a3"]
    )
    assert record.eos.rt_eos.V0 / record.volume_scale == pytest.approx(
        report["cold_volume_cell_a3"]
    )
    assert record.reference_volume != pytest.approx(
        report["cold_volume_cell_a3"], abs=0.1
    )


@pytest.mark.parametrize("metal", RECORDS)
def test_diagnostics_require_explicit_acceptance_and_roundtrip_canonical_data(metal):
    record = get_eos_record(RECORDS[metal])
    assert record in search_eos_records(thermal=True)
    document = get_material_document(metal)
    raw = next(
        r for r in document["eos_records"] if r["identifier"] == record.identifier
    )
    assert raw["record_kind"] == "derived"
    assert raw["default"] is False
    assert raw["catalog_access"] == "explicit_selection"
    assert raw["scientific_validation"]["status"] == "not_reproduced"
    assert raw["derivation"]["diagnostic_only"] is True
    assert raw["derivation"]["published_analytical_pvt_reproduced"] is False
    assert raw["derivation"]["experimental_observations"] == 0
    assert raw["derivation"]["author_fit_reproduced"] is False
    assert raw["parameter_covariance"] is None
    assert all(v is None for v in raw["thermal"]["parameter_errors"].values())
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    assert not list(Draft202012Validator(schema).iter_errors(document))
    with pytest.raises(MaterialError, match="not_reproduced"):
        Material.from_eosmat(document, record_identifiers=[record.identifier])
    assert record.identifier not in {
        r.identifier for r in Material.from_eosmat(document).eos_records
    }
    material = Material.from_eosmat(
        document,
        record_identifiers=[record.identifier],
        require_primary_validation=False,
    )
    assert (
        material.to_snapshot_dict()["eos_records"][0]["equation"]["combination"][
            "validated_as_composed"
        ]
        is False
    )
    with pytest.raises(MaterialError, match="Invalid EOS record"):
        Material.from_dict(material.to_snapshot_dict())
    for restored in [
        material.eos_records[0],
        Material.from_eosmat(
            material.to_eosmat(), require_primary_validation=False
        ).eos_records[0],
    ]:
        v = record.reference_volume * np.array([0.65, 0.8, 0.98])[:, None]
        t = np.array([300, 1250, 2750])[None, :]
        assert restored.pressure(v, t) == pytest.approx(record.pressure(v, t), abs=1e-9)
        assert restored.eos.configuration_values() == record.eos.configuration_values()


@pytest.mark.parametrize("metal", RECORDS)
def test_reconstruction_is_bounded_and_roots_cannot_extrapolate(metal):
    record = get_eos_record(RECORDS[metal])
    v0 = record.reference_volume
    for v, t in [
        (v0 * 0.599, 500),
        (v0 * 1.001, 500),
        (v0 * 0.8, -1),
        (v0 * 0.8, 3001),
        (np.nan, 500),
        (v0 * 0.8, np.inf),
    ]:
        with pytest.raises((EosValidationError, ValueError)):
            record.pressure(v, t)
    with pytest.raises(EosValidationError, match="No root"):
        record.volume(float(record.pressure(0.6 * v0, 1000)) + 10, 1000)
    with pytest.raises(EosValidationError, match="No root"):
        record.eos.temperature(
            float(record.pressure(0.8 * v0, 3000)) + 1, 0.8 * v0 * record.volume_scale
        )
    assert record.pressure(v0, 0) == pytest.approx(
        record.eos.rt_eos.pressure(v0 * record.volume_scale)
    )
    # The total model's 300 K isotherm includes absolute phonon and electronic
    # pressure. Its DAC increment must subtract the *full* Tr contribution.
    assert record.thermal_pressure_increment(0.8 * v0, 300) == pytest.approx(
        0, abs=1e-12
    )
    with pytest.raises(UnsupportedOperationError):
        record.eos.molar_heat_capacity_v(0.8 * v0 * record.volume_scale, 1000)


@pytest.mark.parametrize("metal", RECORDS)
def test_intermediate_pressure_surface_is_monotone_and_dac_inversions_use_tr(metal):
    record = get_eos_record(RECORDS[metal])
    eos = record.eos
    v0 = record.reference_volume
    volume = v0 * np.linspace(0.61, 0.99, 19)[:, None]
    temperature = np.array([0, 50, 300, 425, 750, 1250, 1750, 2250, 2750, 3000])[
        None, :
    ]
    pressure = record.pressure(volume, temperature)
    assert np.all(np.diff(pressure, axis=0) < 0)
    assert np.all(np.diff(pressure, axis=1) > 0)
    v = 0.8 * v0 * record.volume_scale
    step = 1e-5
    derivative = -(
        np.log(eos.characteristic_temperature(v * np.exp(step)))
        - np.log(eos.characteristic_temperature(v * np.exp(-step)))
    ) / (2 * step)
    assert derivative == pytest.approx(eos.gruneisen_parameter(v), rel=1e-8)
    # Roundtrip the actual confinement boundary condition, rather than testing
    # an incorrect 0 K cold-curve reference for a 300 K ambient experiment.
    ambient = 0.85 * v0 * record.volume_scale
    reference_pressure = eos.pressure(ambient, 300)
    for f_dac in [0, 0.4]:
        heated = eos.volume_with_dac_confinement(reference_pressure, 1500, f_dac=f_dac)
        assert eos.pressure(heated, 1500) == pytest.approx(
            reference_pressure + f_dac * eos.thermal_pressure_increment(heated, 1500),
            abs=1e-8,
        )
        assert eos.temperature_from_volumes(
            ambient, heated, f_dac=f_dac
        ) == pytest.approx(1500, abs=1e-6)


def test_platinum_correction_is_explicit_and_phase_holdouts_remain_withheld():
    report = platinum_report()
    fit = report["fits"]["primary"]
    assert fit["fit_states"]["states"] == 166
    assert fit["fit_states"]["max_abs_residual_gpa"] < 0.0061
    assert report["volume_interpolation_holdout"]["held_out_states"]["states"] == 79
    assert (
        report["volume_interpolation_holdout"]["held_out_states"][
            "max_abs_residual_gpa"
        ]
        < 0.0061
    )
    assert fit["first_liquid_states"]["states"] == 2
    assert fit["residual_pressure_gpa"][0] == 0
    assert max(fit["residual_pressure_gpa"]) > 0.3
    assert fit["parameters"]["gamma0"] == 2.63
    assert fit["parameters"]["theta0_k"] == 230
    assert (
        report["fits"]["five_parameter_without_pressure_correction"]["fit_states"][
            "max_abs_residual_gpa"
        ]
        > 0.17
    )
    assert "not an electronic-pressure observation" in report["qualification"]
    assert report["adaptive_quad_max_difference_gpa"] < 1e-10
    assert all(
        not row["used_in_primary_fit"]
        for row in report["states"]
        if row["source_phase_annotation"] != "unmarked"
    )
    for start in report["multistart"]:
        assert start["fitted"]["a"] == pytest.approx(fit["parameters"]["a"], abs=1e-6)
        assert start["fitted"]["b"] == pytest.approx(fit["parameters"]["b"], abs=1e-5)
    check_reconstruction(json.loads(PT_OUTPUT.read_text(encoding="utf-8")), report)
    with PT_RESIDUALS.open() as stream:
        archived = list(csv.DictReader(stream))
    assert len(archived) == 168
    assert [r["source_phase_annotation"] for r in archived] == [
        r["source_phase_annotation"] for r in report["states"]
    ]


def test_platinum_withheld_outputs_do_not_define_their_cold_volume_or_correction():
    report = platinum_report()
    rows = report["states"]
    ratio, temp, pressure = [
        np.array([float(r[key]) for r in rows])
        for key in ["volume_ratio", "temperature_k", "pressure_gpa"]
    ]
    unmarked = np.array([r["source_phase_annotation"] == "unmarked" for r in rows])
    withheld = unmarked & np.isin(ratio, np.unique(ratio)[1::2])
    training = unmarked & ~withheld
    fitted, correction = fit_with_pressure_correction(ratio, temp, pressure, training)
    changed = pressure.copy()
    changed[~training] += 100
    control, control_correction = fit_with_pressure_correction(
        ratio, temp, changed, training
    )
    assert fitted == control
    assert correction[1] == pytest.approx(control_correction[1], abs=1e-12)


def test_residual_pressure_table_validation_and_reconstructable_parameters():
    record = get_eos_record(RECORDS["platinum"])
    eos = record.eos
    restored = eos.with_parameters(gamma0=eos.gamma0)
    assert restored.configuration_values() == eos.configuration_values()
    assert restored.pressure(0.8 * eos.ambient_reference_volume, 1250) == pytest.approx(
        eos.pressure(0.8 * eos.ambient_reference_volume, 1250)
    )
    original = get_material_document("platinum")
    for key, value in [
        ("residual_pressure_gpa", [0]),
        ("residual_temperature_k", [0, 0]),
        ("temperature_range_k", [0, 5001]),
    ]:
        document = copy.deepcopy(original)
        raw = next(
            r for r in document["eos_records"] if r["identifier"] == record.identifier
        )
        raw["thermal"]["configuration"][key] = value
        with pytest.raises(EosmatError):
            validate_eosmat_document(document)


def test_saved_library_audit_is_reproducible():
    check_reconstruction(
        json.loads(OUTPUT.read_text(encoding="utf-8")), validate_library()
    )
