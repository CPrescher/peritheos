"""Shared assertions and isolated results for numerical audit replays."""

import copy

import pytest


@pytest.fixture(scope="session")
def _reproduction_cache():
    return {}


@pytest.fixture
def reproduce_once(_reproduction_cache):
    """Compute an unchanged, no-argument audit once; return an isolated copy.

    Use only for fixed source inputs. Tests that modify inputs, patch a solver,
    or check repeated execution must call the reproduction function directly.
    Keep this cache in tests: production audits must always read current inputs.
    """

    def result(reproduce):
        if reproduce not in _reproduction_cache:
            _reproduction_cache[reproduce] = reproduce()
        return copy.deepcopy(_reproduction_cache[reproduce])

    return result


@pytest.fixture
def assert_audit_close():
    """Compare report structure exactly and fitted floats within solver precision.

    Finite-difference refits vary slightly across NumPy/SciPy and CPU versions.
    Identifiers, counts, flags, and dictionary/list structure remain exact;
    only floating-point leaves receive the explicitly bounded tolerance.
    """

    def compare(actual, expected, *, rel=1e-5, abs=1e-10, path="report"):
        assert type(actual) is type(expected), path
        if isinstance(expected, dict):
            assert actual.keys() == expected.keys(), path
            for key in expected:
                compare(
                    actual[key], expected[key], rel=rel, abs=abs, path=f"{path}.{key}"
                )
        elif isinstance(expected, list):
            assert len(actual) == len(expected), path
            for i, (left, right) in enumerate(zip(actual, expected)):
                compare(left, right, rel=rel, abs=abs, path=f"{path}[{i}]")
        elif isinstance(expected, float):
            assert actual == pytest.approx(expected, rel=rel, abs=abs), path
        else:
            assert actual == expected, path

    return compare
