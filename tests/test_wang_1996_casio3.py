"""Distinct Wang (1996) source regressions must keep their masks and weights."""

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import get_material_document
from scripts.reproduce_wang_1996_casio3 import (
    DATASET_ID,
    load_table,
    reproduce,
    thermal_fit,
)

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(scope="module")
def audit():
    return reproduce()


def test_room_temperature_mask_recovers_parameters_within_source_errors(audit):
    fit = audit["room_temperature"]
    assert fit["observations"] == 12
    assert fit["fixed_parameters"] == {"K0_prime": 4.8}
    assert abs(fit["parameters"]["V0_a3"] - 45.58) < 0.04
    assert abs(fit["parameters"]["K0_gpa"] - 232) < 8
    assert fit["pressure_rmse_gpa"] == pytest.approx(0.26188523, abs=1e-7)


def test_equal_weight_mao_mask_reproduces_rounded_parameters_and_errors(audit):
    fit = audit["mao_equal_weight_reanalysis"]
    assert fit["observations"] == 44
    assert fit["objective"] == "unit-weight pressure residuals"
    p, e = fit["parameters"], fit["standard_errors"]
    assert round(p["V0_a3"], 2) == 45.71
    assert round(p["K0_gpa"]) == 244
    assert round(p["K0_prime"], 1) == 4.8
    assert round(e["V0_a3"], 2) == 0.12
    assert round(e["K0_gpa"]) == 12
    assert round(e["K0_prime"], 1) == 0.3


def test_pressure_weighted_mao_diagnostics_keep_missing_error_and_mismatch(audit):
    fit = audit["mao_pressure_weighted_excluding_sub1"]
    assert fit["observations"] == 44
    assert fit["status"] == "not_reproduced_at_published_precision"
    assert fit["parameters"]["K0_gpa"] == pytest.approx(279.1467, abs=1e-3)
    assert abs(fit["parameters"]["V0_a3"] - 45.47) > 0.06
    assert audit["mao_pressure_weighted_all"]["status"] == "not_exactly_reconstructible"


def test_literal_weight_singularity_is_not_regularized_or_replaced(audit):
    _, rows = load_table(DATASET_ID)
    assert rows[34]["source_order"] == "35"
    assert rows[34]["differential_stress_gpa"] == "0.000"
    with pytest.raises(ValueError, match="singular"):
        thermal_fit(rows, weighting="literal_reciprocal_sum")
    literal = audit["thermal_pressure"]["run13_literal_reciprocal_sum"]
    assert literal["observations"] == 34
    assert literal["parameters"]["K0_gpa"] == pytest.approx(207.4627, abs=1e-3)
    assert abs(literal["parameters"]["K0_gpa"] - 229) > 4


def test_thermal_masks_and_proxy_protocols_are_distinct(audit):
    fits = audit["thermal_pressure"]
    all_rows, excluded = fits["unweighted_all_66"], fits["unweighted_exclude_low_64"]
    assert all_rows["observations"] == 66
    assert excluded["observations"] == 64
    assert set(all_rows["source_orders"]) - set(excluded["source_orders"]) == {33, 34}
    weighted2 = fits["weighted2_run13_volume_proxy"]
    assert weighted2["weighting"] == "volume_only_proxy"
    assert weighted2["source_orders"] == list(range(1, 35))
    assert weighted2["parameters"]["K0_gpa"] == pytest.approx(229.6981, abs=1e-3)
    assert weighted2["derived"]["dK_dT_p_gpa_per_k"] == pytest.approx(
        -4.8 * weighted2["parameters"]["dP_dT_v_gpa_per_k"]
    )
    assert abs(all_rows["parameters"]["dP_dT_v_gpa_per_k"] - 0.0071) > 0.0001


def test_artifact_and_metadata_preserve_qualifications_and_covariances(audit):
    saved = json.loads((ROOT / "docs/data/wang-1996-casio3-refit.json").read_text())
    assert audit["published_eos_replaced"] is False
    assert audit["overall_status"] == "partial_reproduction_with_explicit_diagnostics"
    for key in ("room_temperature", "mao_equal_weight_reanalysis"):
        assert saved[key]["parameters"] == pytest.approx(
            audit[key]["parameters"], rel=1e-6
        )
        covariance = np.array(audit[key]["covariance"])
        assert np.all(np.linalg.eigvalsh(covariance) > 0)
    doc = get_material_document("ca_perovskite")
    by_id = {r["identifier"]: r for r in doc["eos_records"]}
    published = by_id["ca_perovskite_wang_1996_mao_equal_weight_bm3"]
    assert published["eos"]["parameters"] == {"V0": 45.71, "K0": 244.0, "K0_prime": 4.8}
    assert (
        published["scientific_validation"]["reproduction_status"]
        == "rounded_parameters_and_errors_reproduced"
    )
    nonparity = by_id["ca_perovskite_wang_1996_mao_pressure_weighted_no_sub1gpa_bm3"]
    assert nonparity["eos"]["parameters"] == {"V0": 45.47, "K0": 268.0, "K0_prime": 4.3}
    assert (
        nonparity["scientific_validation"]["reproduction_status"]
        == "not_reproduced_at_published_precision"
    )
