"""Distinct Wang (1996) source regressions must keep their masks and weights."""

import hashlib
import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_material_document, list_eos_record_documents
from scripts.register_wang_1996_thermal_refit import (
    DATASET_ID as REFIT_DATASET_ID,
)
from scripts.register_wang_1996_thermal_refit import (
    RECORD_ID,
    ledger_outcome,
)
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


def test_registered_thermal_refit_has_exact_source_manifest_and_reproduces(audit):
    doc = get_material_document("ca_perovskite")
    record = next(r for r in doc["eos_records"] if r["identifier"] == RECORD_ID)
    _, source = load_table(DATASET_ID)
    dataset, selected = load_table(REFIT_DATASET_ID)
    assert selected == [r for r in source if int(r["source_order"]) not in (33, 34)]
    assert len(selected) == 64
    assert dataset["used_by_eos_records"] == [RECORD_ID]
    assert record["default"] is False
    assert RECORD_ID in list_eos_record_documents()
    assert record["experimental_pressure_range_gpa"] == [2.66, 12.25]
    assert record["experimental_temperature_range_k"] == [301.0, 1594.0]
    assert (
        ledger_outcome(record)["parity_basis"]
        == "registered_peritheos_refit_reproduced"
    )
    assert record["parameter_covariance"]["parameter_order"] == [
        "rt_eos.V0",
        "rt_eos.K0",
        "alpha_KT",
    ]
    fit = audit["thermal_pressure"]["unweighted_exclude_low_64"]
    assert np.array(record["parameter_covariance"]["matrix"]) == pytest.approx(
        np.array(fit["covariance"]), rel=1e-6
    )


def test_selectable_thermal_record_predictions_inverse_and_covariance():
    doc = get_material_document("ca_perovskite")
    record = Material.from_eosmat(doc, record_identifiers=[RECORD_ID]).get_eos_record(
        RECORD_ID
    )
    raw = next(r for r in doc["eos_records"] if r["identifier"] == RECORD_ID)
    _, rows = load_table(REFIT_DATASET_ID)
    v, t, observed = [
        np.array([float(r[key]) for r in rows])
        for key in ("volume_a3", "temperature_k", "pressure_gpa")
    ]
    parameters = np.array(
        [
            raw["eos"]["parameters"]["V0"],
            raw["eos"]["parameters"]["K0"],
            raw["thermal"]["parameters"]["alpha_KT"],
        ]
    )

    def pressure(x, volumes, temperatures):
        v0, k0, slope = x
        eta = v0 / volumes
        return 1.5 * k0 * (eta ** (7 / 3) - eta ** (5 / 3)) * (
            1 + 0.75 * 0.8 * (eta ** (2 / 3) - 1)
        ) + slope * (temperatures - 300)

    predicted = record.pressure(v, t)
    assert predicted == pytest.approx(pressure(parameters, v, t), abs=1e-10)
    assert np.sqrt(np.mean((predicted - observed) ** 2)) == pytest.approx(0.3562013702)
    assert record.volume(predicted, t) == pytest.approx(v, abs=1e-8)
    gradient = []
    for index, step in enumerate([1e-4, 1e-3, 1e-7]):
        delta = np.zeros(3)
        delta[index] = step
        gradient.append(
            (
                pressure(parameters + delta, 45.0, 1000.0)
                - pressure(parameters - delta, 45.0, 1000.0)
            )
            / (2 * step)
        )
    covariance = np.array(raw["parameter_covariance"]["matrix"])
    expected_error = np.sqrt(np.array(gradient) @ covariance @ np.array(gradient))
    assert record.pressure_with_uncertainty(
        45.0, 1000.0
    ).standard_error == pytest.approx(expected_error, rel=1e-5)


def test_thermal_addition_preserves_all_previous_scientific_content():
    doc = get_material_document("ca_perovskite")
    doc["eos_records"] = [r for r in doc["eos_records"] if r["identifier"] != RECORD_ID]
    doc["datasets"] = [
        d for d in doc["datasets"] if d["identifier"] != REFIT_DATASET_ID
    ]
    for record in doc["eos_records"]:
        record.pop("scientific_validation")
    payload = json.dumps(doc, sort_keys=True, separators=(",", ":")).encode()
    assert (
        hashlib.sha256(payload).hexdigest()
        == "0644a72ce490a9ced97297f40cda35ad31d208622fbb2f012bd740e32b1b4b61"
    )
