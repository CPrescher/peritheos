"""Sharing audit work must preserve independent, mutable test results."""

import copy

import numpy as np
import pytest


def test_reproduction_runs_once_and_each_consumer_gets_an_isolated_copy(reproduce_once):
    calls = []

    def reproduce():
        calls.append(1)
        return {"rows": [{"value": 2.0}], "scope": "independent_refit"}

    first = reproduce_once(reproduce)
    first["rows"][0]["value"] = -1.0
    second = reproduce_once(reproduce)
    assert calls == [1]
    assert second == {"rows": [{"value": 2.0}], "scope": "independent_refit"}
    assert first is not second


def test_distinct_reproduction_functions_keep_separate_results(reproduce_once):
    def first():
        return {"coefficient": 1.0}

    def second():
        return {"coefficient": 2.0}

    assert reproduce_once(first) == {"coefficient": 1.0}
    assert reproduce_once(second) == {"coefficient": 2.0}


def test_failed_calculation_is_not_reused(reproduce_once):
    attempts = []

    def reproduce():
        attempts.append(1)
        if len(attempts) == 1:
            raise ValueError("unavailable input")
        return {"observations": 3}

    with pytest.raises(ValueError, match="unavailable input"):
        reproduce_once(reproduce)
    assert reproduce_once(reproduce) == {"observations": 3}
    assert attempts == [1, 1]


def test_parameter_replay_allows_only_refit_roundoff(assert_parameter_replay_close):
    archived = [
        {
            "published": 174.0,
            "published_error": 6.0,
            "refit": 173.78786,
            "refit_error": 6.19633,
            "difference": -0.21214,
            "similar": True,
        }
    ]
    fresh = copy.deepcopy(archived)
    fresh[0]["refit"] = np.float64(173.78785985)
    fresh[0]["difference"] = np.float64(-0.21214015)
    assert_parameter_replay_close(fresh, archived)
    for key, value in [
        ("published", 174.00000001),
        ("published_error", 6.00000001),
        ("similar", False),
        ("refit", 174.0),
    ]:
        changed = copy.deepcopy(fresh)
        changed[0][key] = value
        with pytest.raises(AssertionError):
            assert_parameter_replay_close(changed, archived)


def test_assessment_replay_preserves_fixed_values_counts_and_decisions(
    assert_assessment_close,
):
    archived = {
        "status": "parity",
        "fit": {
            "observations": 12,
            "fixed_values": {"V0": 155.28},
            "rmse_gpa": 0.257416842,
            "parameter_covariance": [[38.3945126]],
        },
    }
    fresh = copy.deepcopy(archived)
    fresh["fit"]["rmse_gpa"] += 1e-10
    assert_assessment_close(fresh, archived)
    for key, value in [("observations", 13), ("fixed_values", {"V0": 155.28000001})]:
        changed = copy.deepcopy(fresh)
        changed["fit"][key] = value
        with pytest.raises(AssertionError):
            assert_assessment_close(changed, archived)
    fresh["status"] = "not_reproduced"
    with pytest.raises(AssertionError):
        assert_assessment_close(fresh, archived)


def test_record_replay_preserves_published_coefficients_and_workflow_claim(
    assert_record_replay_close,
):
    archived = {
        "parameters": {"K0": {"value": 174.0}},
        "evidence": {
            "reproduction": {
                "assessment": {"status": "parity"},
                "author_workflow_reproduced": False,
                "source_equation": {
                    "predicted_volume_a3": [10.0],
                    "pixel_checkpoints": [[142.5, 131.0]],
                    "tolerance_a3": 0.41,
                },
            }
        },
    }
    fresh = copy.deepcopy(archived)
    fresh["evidence"]["reproduction"]["source_equation"]["predicted_volume_a3"][0] += (
        1e-12
    )
    assert_record_replay_close(fresh, archived)
    for field, value in (
        ("pixel_checkpoints", [[142.500000000001, 131.0]]),
        ("tolerance_a3", 0.41000000001),
    ):
        changed = copy.deepcopy(fresh)
        changed["evidence"]["reproduction"]["source_equation"][field] = value
        with pytest.raises(AssertionError):
            assert_record_replay_close(changed, archived)
    changed = copy.deepcopy(fresh)
    changed["parameters"]["K0"]["value"] += 1e-10
    with pytest.raises(AssertionError):
        assert_record_replay_close(changed, archived)
    fresh["evidence"]["reproduction"]["author_workflow_reproduced"] = True
    with pytest.raises(AssertionError):
        assert_record_replay_close(fresh, archived)
