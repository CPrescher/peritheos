"""Primary table, pressure-scale and extrapolation checks for Sakai (2011)."""

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import (
    find_pressure_calibration_path,
    get_eos_record,
    get_eos_record_document,
    get_material_document,
    recalculate_eos_pressure_scale,
    recalculate_xrd_pressure_scale,
)
from peritheos.eosmat import validate_eosmat_document
from peritheos.materials import Material
from scripts.compare_sakai_2011_pressure_scales import STANDARDS, compare_pair
from scripts.reproduce_sakai_2011_nacl_b2 import (
    PREFIX,
    ledger_outcome,
    load_rows,
    marker_pressure,
    reference_pressure,
    reproduce,
)


def source_records():
    return [
        record
        for record in get_material_document("nacl_b2")["eos_records"]
        if record["identifier"].startswith(PREFIX)
    ]


def test_primary_transcription_keeps_all_observations_errors_and_markers():
    current = load_rows("nacl-b2-sakai-2011-table1.csv")
    older = load_rows("nacl-b2-sata-2002-table1-sakai-input.csv")
    assert len(current) == 27
    assert len(older) == 29
    assert len({row["run"] for row in current}) == 27
    assert current[0]["run"] == "PNW01_007"
    assert float(current[0]["pt_volume_a3"]) == 52.59
    assert float(current[0]["volume_a3"]) == 24.18
    assert float(current[-1]["pressure_matsui_gpa"]) == 304.4
    assert float(current[-1]["pressure_holmes_gpa"]) == 333.1
    assert float(current[-2]["nacl_a_error_angstrom"]) == 0.0
    assert float(current[-2]["volume_error_a3"]) == 0.02
    assert [r["source_row"] for r in older if r["excluded_by_sakai"] == "1"] == [
        "8",
        "11",
    ]
    assert sum(r["annealed"] == "1" for r in older) == 16
    assert float(older[1]["pt_d111_angstrom"]) == 2.19064
    assert float(older[1]["pt_d111_error_angstrom"]) == 0.00008
    assert float(older[-1]["volume_error_a3"]) == 0.08


def test_source_owned_records_preserve_fixed_extrapolated_volume_and_scales():
    document = get_material_document("nacl_b2")
    validate_eosmat_document(document)
    records = source_records()
    assert len(records) == 8
    assert document["formula_units_per_cell"] == 1
    assert document["space_group_number"] == 221
    assert all(r["record_kind"] == "published" for r in records)
    restored = Material.from_eosmat(document).to_eosmat()
    by_id = {r["identifier"]: r for r in restored["eos_records"]}
    for record in records:
        assert record["fixed_parameters"] == ["V0"]
        assert record["parameter_error_confidence"] is None
        assert record["parameter_covariance"] is None
        assert record["parameter_errors"]["V0"] > 3
        assert (
            by_id[record["identifier"]]["scientific_validation"]
            == record["scientific_validation"]
        )
        assert (
            by_id[record["identifier"]]["pressure_calibration"]
            == record["pressure_calibration"]
        )
    matsui = next(r for r in records if r["identifier"] == PREFIX + "matsui_pt_bm3")
    assert matsui["eos"]["parameters"] == {
        "V0": 37.73,
        "K0": 47.00,
        "K0_prime": 4.10,
    }
    assert matsui["parameter_errors"] == {
        "V0": 4.05,
        "K0": 0.46,
        "K0_prime": 0.02,
    }
    assert matsui["pressure_calibration"]["methods"][0]["reference"]["doi"] == (
        "10.1063/1.3054331"
    )
    assert "nacl_b2_sakai_2014_yokoo_pt_bm3" in by_id


@pytest.mark.parametrize("source", source_records(), ids=lambda r: r["identifier"])
def test_every_scale_executes_inverts_and_reproduces_primary_measurements(source):
    record = get_eos_record(source["identifier"])
    v0, k0, kp = (source["eos"]["parameters"][k] for k in ("V0", "K0", "K0_prime"))
    assert record.pressure(v0) == pytest.approx(0, abs=1e-10)
    assert record.eos.bulk_modulus(v0) == pytest.approx(k0)
    low, high = source["experimental_pressure_range_gpa"]
    p = np.linspace(low + 1, high - 1, 9)
    volumes = record.volume(p, check_validity=True)
    np.testing.assert_allclose(record.pressure(volumes), p, atol=1e-8)
    model = source["identifier"].split("_pt_")[1]
    np.testing.assert_allclose(
        record.pressure(volumes),
        reference_pressure(volumes, v0, k0, kp, model),
        atol=1e-10,
    )
    with pytest.raises(ValueError, match="validity"):
        record.volume(364, check_validity=True)
    with pytest.raises(ValueError, match="isothermal"):
        record.volume(100, 1500, check_validity=True)
    result = reproduce()["records"][source["identifier"]]
    assert result["table1_points_within_three_printed_volume_error_widths"] == 27
    assert result["catalog_max_difference_from_independent_equation_gpa"] < 1e-10
    assert result["fixed_V0_refit"]["observations"] == 54
    assert abs(result["staged_gG_refit"]["V0"] - v0) <= source["parameter_errors"]["V0"]
    outcome = ledger_outcome(source)
    assert outcome["status"] == "parity"
    assert all(p["within_reported_error"] for p in outcome["parameters"])
    assert all(p["within_combined_2sigma"] is None for p in outcome["parameters"])


def test_calibration_reduction_and_independent_published_extrapolation():
    result = reproduce()
    for scale in ("matsui", "fei", "do"):
        assert (
            result["calibrant_checks"][scale]["max_abs_pressure_difference_gpa"] < 0.17
        )
    # This unresolved convention difference must not be hidden by linking to the
    # exact Holmes implementation or silently replacing its coefficients.
    assert result["calibrant_checks"]["holmes"]["max_abs_pressure_difference_gpa"] > 0.6
    assert marker_pressure(60.38 * 0.70, "matsui") == pytest.approx(235.96, abs=0.005)
    comparison = result["independent_published_comparisons"]
    assert comparison["matsui_vinet_minus_bm3_gpa"] == pytest.approx(-5, abs=0.01)
    for scale, source in comparison["source_reported_difference_gpa"].items():
        low, high = comparison["half_last_digit_coefficient_rounding_bounds_gpa"][scale]
        assert low <= source <= high
    assert comparison["bm3_scale_difference_gpa"]["holmes"] > 39
    assert comparison["source_reported_difference_gpa"]["holmes"] == 37
    matsui = result["records"][PREFIX + "matsui_pt_bm3"]
    assert matsui["fixed_V0_refit"]["K0"] == pytest.approx(47.14269657, abs=1e-5)
    assert matsui["volume_error_weighted_refit"]["K0"] > 48


def test_saved_reproduction_retains_qualifications():
    path = (
        Path(__file__).resolve().parents[1]
        / "docs/data/sakai-2011-nacl-b2-reproduction.json"
    )
    stored = json.loads(path.read_text())
    assert set(stored["records"]) == {r["identifier"] for r in source_records()}
    assert (
        stored["independent_published_comparisons"][
            "reported_matsui_model_difference_gpa"
        ]
        == 5
    )


@pytest.mark.parametrize("scale,k0,kp", [("matsui", 273.0, 5.20), ("fei", 277.0, 5.08)])
def test_exact_platinum_reference_isotherms_and_links(scale, k0, kp):
    document = get_material_document("platinum")
    validate_eosmat_document(document)
    source = get_eos_record_document(STANDARDS[scale])
    assert source["eos"]["parameters"] == {"V0": 60.38, "K0": k0, "K0_prime": kp}
    assert source["temperature_ref"] == 300.0
    assert source["parameter_error_confidence"] is None
    assert "thermal" not in source
    restored = Material.from_eosmat(document).to_eosmat()
    record = next(
        r for r in restored["eos_records"] if r["identifier"] == STANDARDS[scale]
    )
    assert record["scientific_validation"] == source["scientific_validation"]
    assert record["pressure_calibration"] == source["pressure_calibration"]
    pt = get_eos_record(STANDARDS[scale])
    volume = pt.volume(np.array([20.0, 60.0, 90.0]), 300.0, check_validity=True)
    np.testing.assert_allclose(
        pt.pressure(volume, 300.0), [20.0, 60.0, 90.0], atol=1e-8
    )
    for model in ("bm3", "vinet"):
        rid = PREFIX + scale + "_pt_" + model
        calibration = get_eos_record_document(rid)["pressure_calibration"]
        assert calibration["methods"][0]["reference_eos_record"] == STANDARDS[scale]
        assert calibration["recalculation"]["status"] == "ready"
    with pytest.raises(ValueError, match="isothermal"):
        pt.pressure(volume, 1000.0)
    with pytest.raises(ValueError, match="validity"):
        pt.volume(310.0, 300.0, check_validity=True)


def test_curve_conversion_uses_pt_state_and_remains_distinct_from_target_fit():
    # Independent Vinet inversion verifies the physical marker transformation,
    # rather than taking the published target NaCl curve as the expected answer.
    from scipy.optimize import brentq

    volumes = np.array([23.0, 20.0, 16.0])
    source_id = PREFIX + "matsui_pt_bm3"
    source_p = reference_pressure(volumes, 37.73, 47.0, 4.10, "bm3")
    pt_volume = np.array(
        [
            brentq(
                lambda v: reference_pressure(v, 60.38, 273, 5.20, "vinet") - p,
                30,
                60.38,
            )
            for p in source_p
        ]
    )
    expected = reference_pressure(pt_volume, 60.38, 277, 5.08, "vinet")
    result = recalculate_eos_pressure_scale(source_id, volumes, STANDARDS["fei"], 300.0)
    assert result.calibration_path == (source_id, STANDARDS["matsui"], STANDARDS["fei"])
    np.testing.assert_allclose(result.target_pressure_gpa, expected, atol=1e-8, rtol=0)
    target_fit = get_eos_record(PREFIX + "fei_pt_bm3").pressure(volumes, 300.0)
    assert np.max(np.abs(result.target_pressure_gpa - target_fit)) > 0.02
    back = recalculate_xrd_pressure_scale(
        result.target_pressure_gpa, STANDARDS["fei"], STANDARDS["matsui"], 300.0
    )
    np.testing.assert_allclose(back.target_pressure_gpa, source_p, atol=1e-8, rtol=0)
    with pytest.raises(ValueError, match="validity"):
        recalculate_eos_pressure_scale(
            source_id, volumes, STANDARDS["fei"], 300.0, check_validity=True
        )


@pytest.mark.parametrize("model", ["bm3", "vinet"])
def test_comparison_reports_extrapolation_and_unresolved_holmes(model):
    fei = compare_pair("matsui", "fei", model)
    assert fei["sample_envelope"]["max_abs_difference_gpa"] < 0.21
    assert fei["states_outside_linked_envelopes"] > 0
    strict = fei["all_linked_envelopes"]
    assert strict["converted_pressure_range_gpa"][1] <= 93.6
    assert strict["max_abs_difference_gpa"] < 0.051
    holmes = compare_pair("matsui", "holmes", model)
    assert holmes["comparison_status"] == "diagnostic_only_unresolved_holmes_variant"
    with pytest.raises(ValueError, match="No executable pressure-calibration path"):
        find_pressure_calibration_path(
            PREFIX + "holmes_pt_" + model, STANDARDS["matsui"]
        )


def test_saved_scale_matrix_retains_all_published_variants_and_controls():
    path = (
        Path(__file__).resolve().parents[1]
        / "docs/data/sakai-2011-pressure-scale-comparison.json"
    )
    result = json.loads(path.read_text())
    assert len(result["comparisons"]) == 20
    assert {
        c["published_target_eos_record"] for c in result["comparisons"].values()
    } == {r["identifier"] for r in source_records()}
    for scale in ("matsui", "fei", "do"):
        for model in ("bm3", "vinet"):
            assert (
                result["comparisons"][f"{scale}_to_{scale}_{model}"]["sample_envelope"][
                    "max_abs_difference_gpa"
                ]
                < 1e-8
            )
    for check in result["calibrant_checks"].values():
        assert check["max_difference_from_independent_vinet_gpa"] < 1e-10
        assert check["max_difference_from_sakai_marker_pressures_gpa"] < 0.17
    assert result["calibrant_checks"]["matsui"][
        "table3_300k_0p7V0_pressure_gpa"
    ] == pytest.approx(235.96, abs=0.005)
