"""Recovered primary-method calibration and conditional source regressions."""

import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from peritheos import (
    Material,
    eosmat_schema,
    get_material_document,
    recalculate_ruby_pressure,
    resolve_dataset_pressure,
)
from scripts.audit_kcl_source_gaps import assert_report_equal, audit

ROOT = Path(__file__).resolve().parents[1]
ID = "kcl_campbell_1991_bm2_1"
DATA = "kcl_campbell_1991_table1_compression"


def test_campbell_primary_scale_and_spatial_errors_survive_selection_and_export():
    document = get_material_document("kcl")
    result = resolve_dataset_pressure(document, ID, DATA, check_validity=True)
    assert result.dataset.checksum_verified
    assert len(result.pressure_gpa) == 14
    np.testing.assert_array_equal(result.pressure_gpa, result.reported_pressure)
    np.testing.assert_array_equal(result.row_indices, np.arange(14))
    assert result.source_calibration["identifier"] == "ruby_mao_1978"
    assert result.target_calibration["parameters"]["B"] == 5
    assert not result.source_pressure_extrapolated.any()
    assert result.uncertainty_treatment == "not_propagated"
    for name in ("pressure_uncertainty_gpa", "normalized_volume_uncertainty"):
        assert result.dataset.get_column(name).role == "standard_deviation"
    assert result.dataset["pressure_uncertainty_gpa"][0] == 0.07
    assert "not the standard error" in result.dataset.notes
    selected = Material.from_eosmat(document, record_identifiers=[ID]).to_eosmat()
    jsonschema.validate(selected, eosmat_schema())
    selected_pressure = resolve_dataset_pressure(selected, ID, DATA)
    np.testing.assert_array_equal(result.pressure_gpa, selected_pressure.pressure_gpa)
    assert selected_pressure.reduction == result.reduction
    assert selected["eos_records"][0]["eos"]["parameters"] == {
        "V0": 52.899988,
        "K0": 28.7,
    }
    assert selected["eos_records"][0]["parameter_errors"] == {
        "V0": 0.355452,
        "K0": 0.6,
    }


def test_campbell_mean_coordinate_normalization_uses_1978_not_1986():
    result = resolve_dataset_pressure(get_material_document("kcl"), ID, DATA)
    ratio = (1 + 5 * result.pressure_gpa / 1904) ** (1 / 5)
    expected = 1904 / 7.665 * (ratio**7.665 - 1)
    actual = recalculate_ruby_pressure(
        result.pressure_gpa, "ruby_mao_1978", "ruby_mao_1986"
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13)
    assert actual[-1] > result.pressure_gpa[-1] + 2
    # Nonlinear transformation cannot commute with five-position averaging.
    five = np.array([45.8, 50.8, 55.8, 60.8, 65.8])
    converted = recalculate_ruby_pressure(five, "ruby_mao_1978", "ruby_mao_1986")
    assert abs(converted.mean() - actual[-1]) > 0.01


def test_primary_regression_weights_and_precision_gaps_remain_visible():
    report = audit()
    c = report["campbell_1991"]
    fits = c["fits"]
    weighted = fits["pressure_and_volume_response_weights"]
    assert weighted["V02_over_V01"] == pytest.approx(0.8483, abs=0.00005)
    assert weighted["K02_gpa"] == pytest.approx(28.7, abs=0.05)
    assert weighted["within_printed_parameter_error_widths"]
    assert fits["pressure_only_normalized_stress_weights"][
        "within_printed_parameter_error_widths"
    ]
    assert not fits["equal_normalized_stress_weights"][
        "within_printed_parameter_error_widths"
    ]
    assert c["exact_source_regression_reproduced"] is False
    t = report["tateno_2019"]
    assert t["workbook_package"] == {
        "sheets": ["Sheet1"],
        "formula_cells": 0,
        "external_links": [],
    }
    assert t["outside_holmes_thermal_approximation_excel_rows"] == [
        25,
        27,
        30,
        32,
        34,
        47,
        49,
    ]
    sensitivity = t["coordinate_precision_sensitivity"]
    assert sensitivity["maximum_absolute_pressure_shift_gpa"] == pytest.approx(
        0.1237483309, abs=1e-8
    )
    assert t["exact_source_regression_reproduced"] is False
    assert all(f["solver_success"] for f in sensitivity["fits"].values())
    assert_report_equal(
        report, json.loads((ROOT / "docs/data/kcl-source-gaps-audit.json").read_text())
    )


def test_source_inventory_and_raw_hashes_are_preserved():
    report = audit()
    assert report["resource_sha256"] == {
        "kcl-campbell-1991-table1-compression.csv": "91671aa26aa19fb20450bf708f260a0e6040c56af1e37d7834f7bbec42a0830f",
        "kcl-tateno-2019-official-table-s1.csv": "f039369c2a8d695c6229f35f2ca8b045171f9bee389e6ad07e0c7d87fe9ea39f",
        "kcl-tateno-2019-table-s1-pvt.csv": "a7e3f6307a851c993bd69cbb8a3c5c424f7e0000826dc7739d56a1d84f228ab1",
    }
    manifest = json.loads(
        (ROOT / "peritheos/data/datasets/kcl_variant_sources/manifest.json").read_text()
    )
    content = [
        m for m in manifest["tateno_deposit_inventory"] if m["kind"] != "metadata"
    ]
    assert {m["kind"] for m in content} == {"workbook", "supplemental_figure_pdf"}
    assert len(content) == 2
    assert len(manifest["tateno_pt_pairing_corrections"]) == 110
