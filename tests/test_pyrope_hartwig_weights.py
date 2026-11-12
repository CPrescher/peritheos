"""Weight semantics, unchanged inputs, and conditional solver replay."""

import hashlib
import json

import numpy as np
import pytest

from scripts.audit_pyrope_angel_2022_observations import DATA
from scripts.audit_pyrope_hartwig_weights import fit_rows, weighted_rows
from scripts.refit_pyrope_angel_2022_expanded import EXPANDED_NAMES, expanded_rows
from scripts.refit_pyrope_angel_2022_subset import arrays, effective_sigma


def test_weight_multiplier_changes_only_copied_hartwig_variance():
    rows = expanded_rows(True)
    original = json.dumps(rows, sort_keys=True)
    report = json.loads((DATA / "expanded-refits.json").read_text())
    c = np.array([report["published_parameters"][k] for k in EXPANDED_NAMES])
    changed = weighted_rows(rows, 9)
    assert json.dumps(rows, sort_keys=True) == original
    for before, after in zip(rows, changed):
        assert before is not after
        if before["dataset"] == "hartwig_vt":
            assert after["sY"] == before["sY"] / 3
            assert after["sT"] == before["sT"] / 3
            for key in before.keys() - {"sY", "sT"}:
                assert after[key] == before[key]
        else:
            assert before == after
    ratio = (
        effective_sigma(c, arrays(changed)) ** 2 / effective_sigma(c, arrays(rows)) ** 2
    )
    np.testing.assert_allclose(
        ratio, [1 / 9 if r["dataset"] == "hartwig_vt" else 1 for r in rows], rtol=1e-12
    )


@pytest.mark.parametrize("value", [0, -1, float("nan"), float("inf")])
def test_invalid_multiplier_is_rejected(value):
    with pytest.raises(ValueError):
        weighted_rows([], value)


def test_sensitivity_manifest_and_two_start_convergence():
    report = json.loads((DATA / "hartwig-weight-sensitivity.json").read_text())
    assert report["author_workflow_reproduced"] is False
    for name, digest in report["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    for name, key in [
        ("hartwig-weight-plan.json", "plan_sha256"),
        ("expanded-refits.json", "baseline_fit_sha256"),
    ]:
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == report[key]
    plan = json.loads((DATA / "hartwig-weight-plan.json").read_text())
    baseline = json.loads((DATA / "expanded-refits.json").read_text())
    for count in ("182", "191"):
        cases = report["cases"][count]
        assert [r["weight_multiplier"] for r in cases] == plan["weight_multipliers"]
        for row in cases:
            for key in ("fit_from_published_start", "fit_from_continuation_start"):
                fit = row[key]
                assert fit["numerically_converged"]
                assert fit["scaled_jacobian_rank"] == 7
                assert not fit["at_parameter_bound"]
            delta = row["continuation_minus_published_start"]
            assert abs(delta["theta0_k"]) < 0.001
            assert abs(delta["gamma0"]) < 1e-6
        for key in EXPANDED_NAMES:
            assert (
                cases[0]["fit_from_published_start"]["parameters"][key]
                == baseline["cases"][count]["refit"]["parameters"][key]
            )


@pytest.mark.parametrize("weight", [500, 1000])
def test_independent_weighted_fit_and_fresh_official_receipt(weight):
    report = json.loads((DATA / "hartwig-weight-sensitivity.json").read_text())
    case = next(r for r in report["cases"]["191"] if r["weight_multiplier"] == weight)
    baseline = json.loads((DATA / "expanded-refits.json").read_text())
    fit = fit_rows(
        weighted_rows(expanded_rows(True), weight),
        [baseline["published_parameters"][key] for key in EXPANDED_NAMES],
    )
    np.testing.assert_allclose(
        list(fit["parameters"].values()),
        list(case["fit_from_published_start"]["parameters"].values()),
        rtol=1e-9,
    )
    native = json.loads(
        (DATA / "hartwig-weight-eosfit" / str(weight) / "checks.json").read_text()
    )
    assert native["Hartwig_inverse_variance_multiplier"] == weight
    assert native["numerical_convergence_reported"]
    # The original strict solver-comparison failure must remain visible.
    assert native["within_declared_tolerances"] is False
    assert abs(native["parameters"]["theta0_k"] - fit["parameters"]["theta0_k"]) < 1


def test_same_coefficient_equations_match_official_calculator_across_states():
    report = json.loads(
        (DATA / "hartwig-weight-eosfit/equation-checks/checks.json").read_text()
    )
    assert report["n_coefficient_sets"] == 3
    assert report["n_states_per_coefficient_set"] == 9
    assert report["all_within_tolerances"]
    assert set(report["cases"]) == {"published", "weight_500", "weight_1000"}
    for case in report["cases"].values():
        assert len(case["states"]) == 9
        assert min(s["temperature_k"] for s in case["states"]) == 103
        assert max(s["temperature_k"] for s in case["states"]) == 1700
        for state in case["states"]:
            for key, limit in report["absolute_tolerances"].items():
                assert abs(state["official_minus_independent"][key]) < limit
