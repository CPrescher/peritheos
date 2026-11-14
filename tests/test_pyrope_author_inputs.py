"""Author-input preservation, fit replay and independent uncertainty checks."""

import csv
import hashlib
import json
from collections import Counter

import numpy as np

from scripts.audit_pyrope_author_inputs import DATA, NAMES, PREFIX, check_covariance
from scripts.refit_pyrope_angel_2022_subset import solve


def audit():
    return json.loads((DATA / (PREFIX + "-audit.json")).read_text())


def test_all194_exported_inputs_and_explicit_scale_assignment():
    report = audit()
    with (DATA / (PREFIX + "-inputs.csv")).open() as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 194
    assert Counter(row["dataset"] for row in rows) == report["counts"]
    assert Counter(row["group"] for row in rows) == {"1": 115, "2": 64, "3": 15}
    assert [int(row["fit_row"]) for row in rows if row["group"] == "3"] == list(
        range(118, 133)
    )
    assert all(row["use"] == "1" for row in rows)
    for row in rows:
        for key in ["T", "P", "Y", "sT", "sP", "sY"]:
            assert float(row[key]) == float(row[key + "_token"])
        assert row["scale_index"] == ("4" if row["group"] == "3" else "")
        if row["dataset"] == "hartwig_vt":
            assert float(row["sY"]) == 0.01
        if row["dataset"] == "milani_vt":
            assert float(row["sT"]) == 7
    # Old conditional reconstruction sources remain byte-identical.
    old = json.loads((DATA / "expanded-refits.json").read_text())
    for name, digest in old["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest


def test_fit_replay_covariance_and_residual_rescaling():
    report = audit()
    a, b = report["independent_fits"]
    assert a["numerically_converged"] and b["numerically_converged"]
    assert not a["at_parameter_bound"] and not b["at_parameter_bound"]
    np.testing.assert_allclose(
        list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-8
    )
    rows = [
        {
            key: value
            for key, value in row.items()
            if not key.startswith(("calculated_", "effective_", "standardized_"))
        }
        for row in a["rows"]
    ]
    replay = solve(
        rows,
        [170, 4.5, 800, 1.2, 0.98],
        list(range(5)),
        names=NAMES,
        scaling=np.array([170, 4.5, 770, 1.2, 1]),
        lower=np.array([50, 2, 50, 0.05, 0.9]),
        upper=np.array([350, 10, 3000, 4, 1.1]),
    )
    np.testing.assert_allclose(
        list(replay["parameters"].values()), list(a["parameters"].values()), rtol=1e-8
    )
    check = check_covariance(rows, a)
    assert check["covariance_positive_definite"]
    assert check["degrees_of_freedom"] == 194 - 5
    for step in check["finite_difference_checks"]:
        assert step["maximum_width_relative_difference"] < 1e-7
        assert step["scaled_inverse_information_maximum_error"] < 1e-10
    assert abs(check["correlation"][2][3] - 0.861078483) < 1e-8
    assert (
        abs(
            check["standard_error_multiplier"]
            - np.sqrt(a["conditional_reduced_chi_square"])
        )
        < 1e-12
    )


def test_native_snapshot_comparison_and_property_checks():
    report = json.loads((DATA / (PREFIX + "-eosfit.json")).read_text())
    assert report["eosfit_version"] == "7.60"
    assert report["eosfit_replay"]["converged"]
    assert report["final194_observation_error_and_group_tokens_match_author_log"]
    assert not report["native_exact_final_log_value_match"]
    for row in report["comparison"]:
        assert row["eosfit_difference_over_author_esd"] < 0.01
        assert row["independent_difference_over_author_esd"] < 0.01
        assert abs(row["independent_error_relative_to_eosfit"]) < 0.005
    native = report["eosfit_replay"]
    c = np.array(native["saved_rescaled_covariance"])
    assert np.linalg.eigvalsh(c).min() > 0
    np.testing.assert_allclose(
        np.sqrt(np.diag(c)), list(native["saved_covariance_widths"].values())
    )
    for name in NAMES:
        r = native["saved_covariance_rounding_checks"][name]
        assert (
            r["width_lower"] - 5.1e-6
            <= native["residual_scaled_esds"][name]
            <= r["width_upper"] + 5.1e-6
        )
    properties = report["fixed_parameter_property_check"]
    assert properties["all_within_tolerance"] and len(properties["states"]) == 12
    for row in properties["states"]:
        for key, tolerance in properties["tolerances"].items():
            assert abs(row["differences"][key]) <= tolerance
