"""Density EOS mechanics and preservation of acoustic provenance."""

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, eosmat_schema, get_material_document
from peritheos.eos.rt import DensityPolynomial3
from scripts.reproduce_grimsditch_1986_argon import (
    MASS_FACTOR,
    RECORD,
    bulk,
    pressure,
    reproduce,
)


def test_exact_density_polynomial_and_nonzero_anchor():
    model = DensityPolynomial3(MASS_FACTOR / 2, 2, 12.65, -11.43, 1.5, 0.68)
    assert model.pressure(model.V0) == pytest.approx(1.23)
    rho = np.array([2.017, 2.5, 3.0, 4.0, 4.8])
    v = MASS_FACTOR / rho
    assert model.pressure(v) == pytest.approx(pressure(rho), abs=5e-13)
    assert model.bulk_modulus(v) == pytest.approx(bulk(rho), abs=5e-13)
    step = v * 1e-5
    derivative = -v * (model.pressure(v + step) - model.pressure(v - step)) / (2 * step)
    assert model.bulk_modulus(v) == pytest.approx(derivative, rel=1e-8)
    assert model.volume(model.pressure(v)) == pytest.approx(v, abs=2e-11)
    with pytest.raises(ValueError):
        model.pressure(MASS_FACTOR / 1.0)  # Unstable branch, not an alternate root.
    with pytest.raises(ValueError):
        model.volume(0.0)  # The source polynomial has no ambient solid state.
    with pytest.raises(ValueError):
        DensityPolynomial3(100, 1, 12.65, -11.43, 1.5, 0.68)
    with pytest.raises(ValueError):
        DensityPolynomial3(100, 2, float("nan"), -11.43, 1.5, 0.68)


def test_eosmat_roundtrip_and_provenance():
    doc = get_material_document("argon_fcc")
    Draft202012Validator(eosmat_schema()).validate(doc)
    material = Material.from_eosmat(doc)
    record = material.get_eos_record(RECORD)
    restored = Material.from_eosmat(material.to_eosmat()).get_eos_record(RECORD)
    assert restored.pressure(90.0) == pytest.approx(record.pressure(90.0))
    raw = next(r for r in doc["eos_records"] if r["identifier"] == RECORD)
    assert raw["reproduction"]["fit_status"] == "not_reproduced"
    assert (
        raw["reproduction"]["summary_status"]
        == "approximate_internal_consistency_reproduced"
    )
    summary = raw["reproduction"]["equal_weight_consistency"]
    diagnostic = reproduce()["equal_weight_diagnostics"][summary["selected_diagnostic"]]
    assert summary["included_row_count"] == diagnostic["row_count"] == 74
    assert summary["excluded_source_rows"] == [98]
    assert summary["coefficients_c0_to_c3"] == pytest.approx(
        diagnostic["coefficients_c0_to_c3"]
    )
    assert summary["pressure_rms_gpa"] == pytest.approx(diagnostic["pressure_rms_gpa"])
    assert summary["all_rows_diagnostic_retained"] is True
    assert summary["published_coefficients_retained"] is True
    assert raw.get("fit_datasets", []) == []
    assert raw["temperature_ref_provenance"]["kind"] == "assumed_room_temperature"
    assert raw["eos"]["parameters"]["c3"] == 0.68


def test_source_table_consistency_with_honest_outlier():
    report = reproduce()
    assert report["table1_count"] == 98
    assert report["liquid_count"] == 23
    assert report["fcc_count"] == 75
    assert report["independent_observed_pv_count"] == 0
    assert report["fit_status"] == "not_reproduced"
    assert report["native_max_abs_pressure_error_gpa"] < 1e-12
    assert report["native_max_abs_bulk_error_gpa"] < 1e-12
    assert report["largest_density_discrepancy"]["source_row"] == 98
    assert report["derived_density_consistency_max_abs_gpa"] > 1.2
    for row in report["table2_checks"]:
        assert (
            abs(row["equation_bulk_gpa"] - row["printed_bulk_gpa"])
            < row["printed_bulk_error_gpa"]
        )


def test_stable_expansion_inverse_near_spinodal():
    eos = DensityPolynomial3(MASS_FACTOR / 2, 2, 12.65, -11.43, 1.5, 0.68)
    spinodal = (-3 + np.sqrt(9 + 4 * 2.04 * 11.43)) / (2 * 2.04)
    density = np.array([spinodal * 1.0001, spinodal * 1.01, 1.8, 1.9])
    volume = MASS_FACTOR / density
    assert eos.volume(eos.pressure(volume)) == pytest.approx(volume, rel=1e-10)
    record = Material.from_eosmat(get_material_document("argon_fcc")).get_eos_record(
        RECORD
    )
    assert record.volume(record.pressure(volume)) == pytest.approx(volume, rel=1e-10)
    with pytest.raises(ValueError):
        eos.volume(pressure(spinodal) - 1e-8)


def test_equal_weight_refits_are_only_derived_data_diagnostics():
    report = reproduce()
    fits = report["equal_weight_diagnostics"]
    all_rows = fits["all_75_solid_rows"]
    sensitivity = fits["sensitivity_without_inconsistent_source_row_98"]
    assert all_rows["row_count"] == 75
    assert sensitivity["row_count"] == 74
    assert set(all_rows["source_rows"]) - set(sensitivity["source_rows"]) == {98}
    assert all_rows["coefficients_c0_to_c3"] == pytest.approx(
        [6.72150182, -4.86088512, -0.88317979, 0.96671525], abs=1e-7
    )
    assert sensitivity["coefficients_c0_to_c3"] == pytest.approx(
        [12.57280259, -11.33786970, 1.46413728, 0.68834664], abs=1e-7
    )
    assert sensitivity["pressure_rms_gpa"] < 0.006
    assert report["fit_status"] == "not_reproduced"
    doc = get_material_document("argon_fcc")
    record = next(r for r in doc["eos_records"] if r["identifier"] == RECORD)
    assert [record["eos"]["parameters"][f"c{i}"] for i in range(4)] == [
        12.65,
        -11.43,
        1.5,
        0.68,
    ]
