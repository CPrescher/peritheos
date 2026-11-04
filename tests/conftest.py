"""Shared assertions and isolated results for numerical audit replays."""

import copy
from fnmatch import fnmatchcase

import numpy as np
import pytest

SOURCE_EVIDENCE_FLOAT_PATHS = {
    "source_equation": ("predicted_volume_a3[[]*]", "difference_a3[[]*]"),
    "source_equation_check": (
        "residuals_km_s[[]*]",
        "rms_km_s",
        "max_abs_residual_over_reported_Us_width",
    ),
    "rankine_hugoniot_replay": (
        "rows[[]*].replayed_pressure_gpa",
        "rows[[]*].pressure_difference_gpa",
        "rows[[]*].replayed_density_g_cm3",
        "rows[[]*].density_difference_g_cm3",
        "max_abs_pressure_difference_gpa",
        "max_abs_density_difference_g_cm3",
    ),
}


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

    def compare(
        actual,
        expected,
        *,
        rel=1e-5,
        abs=1e-10,
        path="report",
        ignore_keys=frozenset(),
        float_paths=None,
    ):
        floating = isinstance(expected, (float, np.floating))
        if floating and (
            float_paths is None
            or path in float_paths
            or any(fnmatchcase(path, p) for p in float_paths)
        ):
            assert isinstance(actual, (float, np.floating)), path
            assert actual == pytest.approx(expected, rel=rel, abs=abs), path
            return
        assert type(actual) is type(expected), path
        if isinstance(expected, dict):
            assert actual.keys() == expected.keys(), path
            for key in expected:
                if key in ignore_keys:
                    continue
                compare(
                    actual[key],
                    expected[key],
                    rel=rel,
                    abs=abs,
                    path=f"{path}.{key}",
                    ignore_keys=ignore_keys,
                    float_paths=float_paths,
                )
        elif isinstance(expected, list):
            assert len(actual) == len(expected), path
            for i, (left, right) in enumerate(zip(actual, expected)):
                compare(
                    left,
                    right,
                    rel=rel,
                    abs=abs,
                    path=f"{path}[{i}]",
                    ignore_keys=ignore_keys,
                    float_paths=float_paths,
                )
        else:
            assert actual == expected, path

    return compare


@pytest.fixture
def assert_parameter_replay_close(assert_audit_close):
    """Allow refit variation while keeping published values and decisions exact."""

    def compare(actual, expected):
        assert_audit_close(
            actual,
            expected,
            float_paths=tuple(
                f"report*.{key}"
                for key in ("refit", "refit_error", "difference", "relative_difference")
            ),
        )

    return compare


@pytest.fixture
def assert_assessment_close(assert_audit_close, assert_parameter_replay_close):
    def compare(actual, expected):
        actual, expected = copy.deepcopy(actual), copy.deepcopy(expected)
        if "parameters" in expected:
            assert_parameter_replay_close(
                actual.pop("parameters"), expected.pop("parameters")
            )
        assert_audit_close(
            actual,
            expected,
            float_paths=(
                "report.fit.rmse_*",
                "report.fit.published_rmse_*",
                "report.fit.parameter_covariance*",
                "report.fit.reduced_chi_square",
            ),
        )

    return compare


@pytest.fixture
def assert_record_replay_close(assert_audit_close, assert_assessment_close):
    """Compare computed evidence separately; all remaining record fields are exact."""

    def compare(actual, expected):
        actual, expected = copy.deepcopy(actual), copy.deepcopy(expected)
        left = actual["evidence"]["reproduction"]
        right = expected["evidence"]["reproduction"]
        assert_assessment_close(left.pop("assessment"), right.pop("assessment"))
        for key in (
            "source_equation",
            "source_equation_check",
            "rankine_hugoniot_replay",
        ):
            if key in right:
                assert_audit_close(
                    left.pop(key),
                    right.pop(key),
                    rel=1e-12,
                    abs=1e-10,
                    float_paths=tuple(
                        f"report.{p}" for p in SOURCE_EVIDENCE_FLOAT_PATHS[key]
                    ),
                )
        if "independent_diagnostic_refit" in right:
            assert_audit_close(
                left.pop("independent_diagnostic_refit"),
                right.pop("independent_diagnostic_refit"),
                abs=2e-7,
                float_paths=(
                    "report.parameters_K0_K0_prime*",
                    "report.difference_from_published*",
                    "report.rms_pressure_gpa",
                ),
            )
        assert actual == expected

    return compare


@pytest.fixture
def assert_hugoniot_report_close(assert_audit_close):
    """Only computed shock diagnostics vary; coefficients and normalization do not."""
    paths = [
        f"report.{section}.{p}"
        for section in ("source_equation_check", "rankine_hugoniot_replay")
        for p in SOURCE_EVIDENCE_FLOAT_PATHS[section]
    ]
    paths += [
        f"report.independent_diagnostic_refit.{p}"
        for p in (
            "parameters_c0_s[[]*]",
            "residual_scaled_standard_errors[[]*]",
            "covariance[[]*]",
            "residuals_km_s[[]*]",
            "rms_km_s",
        )
    ]

    def compare(actual, expected):
        assert_audit_close(actual, expected, rel=1e-12, abs=1e-12, float_paths=paths)

    return compare


@pytest.fixture
def assert_pressure_refit_report_close(assert_audit_close):
    """Preserve measured inputs, fixed anchors and decisions in pressure refits."""

    def compare(actual, expected):
        # The lower endpoint is the unchanged Li observation at 27.2 GPa;
        # only the recalculated Sata upper endpoint is a numerical output.
        paths = ["report.diagnostic_input_pressure_range_gpa[1]"]
        for section, fields in (
            (
                "fixed_V0_BM3",
                (
                    "parameters",
                    "conditional_standard_errors",
                    "conditional_covariance",
                    "residual_gpa",
                    "rmse_gpa",
                    "max_abs_gpa",
                ),
            ),
            (
                "reference_mBM3_check",
                (
                    "parameters_Pr_Kr_Kprime",
                    "standard_errors",
                    "covariance",
                    "residual_gpa",
                    "rmse_gpa",
                    "max_abs_gpa",
                ),
            ),
        ):
            for coordinate in ("original_pressures", "recalculated_Sata_unchanged_Li"):
                for field in fields:
                    prefix = f"report.{section}.{coordinate}.{field}"
                    paths.extend((prefix, prefix + "[[]*]"))
            paths.append(
                f"report.{section}.coefficient_shift_recalculated_minus_original[[]*]"
            )
        for model in expected["cross_evaluation"]:
            for coordinate in expected["cross_evaluation"][model]:
                for field in ("rmse_gpa", "max_abs_gpa", "residual_gpa"):
                    prefix = f"report.cross_evaluation.{model}.{coordinate}.{field}"
                    paths.extend((prefix, prefix + "[[]*]"))
        for i, row in enumerate(expected["inputs"]):
            if row["source"] == 0:
                paths += [
                    f"report.inputs[{i}].{key}"
                    for key in ("fit_pressure_gpa", "pressure_shift_gpa")
                ]
        for i, row in enumerate(expected["residuals"]):
            if row["source"] == 0:
                paths.append(f"report.residuals[{i}].fit_pressure_gpa")
            paths += [
                f"report.residuals[{i}].{model}_residual_{coordinate}_gpa"
                for model in (
                    "published_BM3",
                    "original_pressure_diagnostic_BM3",
                    "recalculated_pressure_diagnostic_BM3",
                )
                for coordinate in ("original", "recalculated")
            ]
        assert_audit_close(actual, expected, rel=1e-11, abs=1e-10, float_paths=paths)

    return compare
