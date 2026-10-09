"""Pressure-only Au validation keeps source replay and output refits distinct."""

import json

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.audit_gold_pvt import (
    OUTPUT,
    PARAMETERS,
    RECORD,
    audit,
    check_saved,
    fei_pressure,
)


def test_published_fei_gold_equation_and_pvt_inversions():
    record = get_eos_record(RECORD)
    volume = np.array([67.85, 61.0, 55.0])[:, None]
    temperature = np.array([300.0, 1473.0, 2173.0])[None, :]
    expected = fei_pressure(volume, temperature)
    assert record.pressure(volume, temperature) == pytest.approx(expected, abs=1e-9)
    assert fei_pressure(PARAMETERS["V0"], 300.0) == pytest.approx(0, abs=1e-12)
    assert record.volume(expected, temperature) == pytest.approx(
        np.broadcast_to(volume, (3, 3)), rel=1e-8
    )
    assert record.eos.temperature(
        expected, volume * record.volume_scale
    ) == pytest.approx(np.broadcast_to(temperature, (3, 3)), rel=1e-8)
    restored = Material.from_eosmat(
        get_material_document("gold"), record_identifiers=[RECORD]
    ).eos_records[0]
    assert restored.pressure(volume, temperature) == pytest.approx(expected, abs=1e-9)


def test_fei_thermal_errors_and_measured_provenance_are_complete():
    document = get_material_document("gold")
    raw = next(r for r in document["eos_records"] if r["identifier"] == RECORD)
    assert raw["thermal"]["debye_temperature_law"] == "variable_exponent"
    assert raw["thermal"]["parameter_errors"]["gamma0"] == 0.03
    assert raw["thermal"]["parameter_errors"]["q"] == 0.3
    assert raw["thermal"]["parameter_errors"]["theta0"] is None
    assert raw["parameter_error_confidence"] is None
    assert raw["parameter_covariance"] is None
    data_check = raw["scientific_validation"]["primary_data_check"]
    assert data_check["status"] == "bundled"
    assert data_check["digitized_dataset_identifiers"] == [
        "gold_fei_2007_figure1_digitized"
    ]
    linked = data_check["dataset_identifiers"]
    assert set(linked) == {
        "gold_fei_2004_table1",
        "gold_dewaele_2004_table1_compression",
        "gold_hirose_2006_table1",
    }
    for identifier in linked:
        data = next(d for d in document["datasets"] if d["identifier"] == identifier)
        assert RECORD in data["used_by_eos_records"]


def test_measured_comparisons_and_calculated_curves_are_separate():
    report = audit()["fei_2007"]
    assert report["independent_equation"]["states"] == 119
    assert report["independent_equation"]["max_abs_difference_gpa"] < 1e-9
    measured = report["measured_data_comparisons"]
    assert measured["fei_2004_hot"]["states"] == 26
    assert measured["dewaele_2004_cold"]["states"] == 37
    assert measured["fei_2004_hot"]["calibration"] == "Speziale (2001) MgO"
    assert 0.43 < measured["fei_2004_hot"]["rmse_gpa"] < 0.44
    assert 0.31 < measured["dewaele_2004_cold"]["rmse_gpa"] < 0.33
    curves = report["source_curve_comparison"]
    assert curves["kind"] == "calculated_source_curves_not_observations"
    assert [c["temperature_k"] for c in curves["curves"]] == [300, 1473, 2173]
    assert all(
        c["states"] == 37 and c["within_one_stroke_pressure_width"]
        for c in curves["curves"]
    )
    assert report["author_fit_reproduced"] is False


def test_yokoo_pressure_reconstruction_has_qualified_validation():
    report = audit()["yokoo_2009"]
    assert report["full_published_pvt_status"] == "not_reproduced"
    assert (
        report["reconstruction_status"]
        == "numerically_verified_derived_output_reconstruction"
    )
    assert report["primary_fit"]["states"] == 156
    assert report["primary_fit"]["max_abs_residual_gpa"] < 0.022
    assert report["volume_holdout"]["states"] == 74
    assert report["volume_holdout"]["max_abs_residual_gpa"] < 0.019
    assert report["first_liquid_holdout"]["states"] == 6
    assert report["independent_quad_max_difference_gpa"] < 1e-10
    assert report["experimental_observations"] == 0
    assert report["published_coefficients_changed"] is False
    assert report["author_fit_reproduced"] is False
    assert "missing electronic energy is not a PVT blocker" in report["limits"]


def test_gold_pvt_audit_is_reproducible():
    check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), audit())
