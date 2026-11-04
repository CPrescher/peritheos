"""Sharing audit work must preserve independent, mutable test results."""

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
