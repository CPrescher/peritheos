"""Estimated widths stay separate from source errors and full-fit claims."""

import csv
import hashlib
import itertools
import json
import math
from collections import Counter

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import (
    DATA,
    GAS_CONSTANT,
    observations,
)
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_full_sensitivity import (
    BASELINE,
    FACTORS,
    PLAN,
    estimated_isaak_width,
    prepared_rows,
    scenarios,
)
from scripts.refit_grossular_angel_2022_normalized import (
    prepared_rows as reference_rows,
)
from scripts.refit_grossular_angel_2022_subset import arrays, pressure, solve


def test_complete_non_anchor_selection_and_explicit_estimates_preserve_originals():
    rows = prepared_rows(scenarios()[0])
    assert len(rows) == 219 and rows[:193] == reference_rows()
    assert Counter(r["dataset"] for r in rows[193:]) == {
        "isaak_ks_estimated": 19,
        "thieblot_vt_estimated": 7,
    }
    missing = {
        int(r["source_row"])
        for r in observations("isaak-1992-table-3b.csv")
        if not r["ks_gpa_reported_error"]
    }
    assert {r["row"] for r in rows if r["dataset"] == "isaak_ks_estimated"} == missing
    assert missing.isdisjoint({1, 11, 21})
    np.testing.assert_allclose(estimated_isaak_width(350, "interpolated"), 0.741)
    np.testing.assert_allclose(estimated_isaak_width(1350, "interpolated"), 1.218)
    for r in rows[193:212]:
        assert r["sT"] is None and r["KS"] and r["scale"] is None
    for r, original in zip(rows[212:], observations("figure4-thieblot-selected.csv")):
        assert r["T"] == float(original["temperature_k"])
        assert r["Y"] == float(original["molar_volume_cm3_mol"])
        assert r["row"] == int(original["drawing_index"])
        np.testing.assert_allclose(r["sY"], math.sqrt(0.03**2 + 0.003**2 / 3))
        np.testing.assert_allclose(r["sT"], math.sqrt((0.01 * r["T"]) ** 2 + 4 / 3))
    for case in scenarios():
        variant = prepared_rows(case)
        assert variant[:193] == rows[:193]
        for old, new in zip(rows, variant):
            assert (old["dataset"], old["row"], old["P"], old["T"], old["Y"]) == (
                new["dataset"],
                new["row"],
                new["P"],
                new["T"],
                new["Y"],
            )
    # Blank source cells are still blank; estimates exist only in diagnostics.
    assert (
        sum(
            not r["ks_gpa_reported_error"]
            for r in observations("isaak-1992-table-3b.csv")
        )
        == 19
    )


def test_predeclared_grid_shapes_convergence_and_parameter_envelopes():
    cases = scenarios()
    assert len(cases) == 29 and cases[0]["id"] == BASELINE
    assert {tuple(c[k] for k in FACTORS) for c in cases[:27]} == set(
        itertools.product([0.5, 1, 2], repeat=3)
    )
    assert {c["isaak_profile"] for c in cases[27:]} == {"constant_low", "constant_high"}
    report = json.loads((DATA / "full-sensitivity-refits.json").read_text())
    assert report["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert report["scientific_cases"] == 29 and report["fitted_rows"] == 219
    assert report["all_starts_converged"] and not report["any_start_at_bound"]
    assert (
        not report["author_workflow_reproduced"]
        and not report["original_errors_recovered"]
    )
    assert report["maximum_scaled_start_difference"] < 1e-6
    for case in report["cases"]:
        a, b = case["fits"]
        for fit in [a, b]:
            assert fit["degrees_of_freedom"] == 213 and len(fit["irls_trace"]) <= 30
            assert np.linalg.eigvalsh(fit["conditional_unscaled_covariance"]).min() > 0
        np.testing.assert_allclose(
            list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
        )
    for name, limits in report["all_case_envelope"].items():
        values = [c["fits"][0]["parameters"][name] for c in report["cases"]]
        assert limits == {"minimum": min(values), "maximum": max(values)}
    for name in ["K0_prime", "theta0_k"]:
        assert (
            report["all_case_envelope"][name]["minimum"]
            > report["comparison"][name]["published"]
        )


def test_baseline_and_grid_corner_replay_and_independent_forward_properties():
    plan = json.loads(PLAN.read_text())
    report = json.loads((DATA / "full-sensitivity-refits.json").read_text())
    for case in [report["cases"][0], report["cases"][1]]:
        replay = solve(
            prepared_rows(case),
            plan["starts"][1],
            names=plan["free_parameters"],
            scaling=np.array(plan["parameter_scaling"]),
            lower=np.array(plan["bounds"]["lower"]),
            upper=np.array(plan["bounds"]["upper"]),
        )
        np.testing.assert_allclose(
            list(replay["parameters"].values()),
            list(case["fits"][0]["parameters"].values()),
            rtol=1e-6,
        )
    p, _ = source("grossular")
    c = list(report["cases"][0]["fits"][0]["parameters"].values())
    fitted = p | dict(zip([2, 3, 11, 18], c[:4]))
    for parameters, coefficients in [
        (p, [p[2], p[3], p[11], p[18], p[53], p[54]]),
        (fitted, c),
    ]:
        rows = []
        for r in prepared_rows(scenarios()[0]):
            prop = properties(parameters, r["P"], r["T"], GAS_CONSTANT)
            multiplier = coefficients[r["scale"]] if r["scale"] is not None else 1
            rows.append(
                r | {"Y": multiplier * prop["ks_gpa" if r["KS"] else "volume_cm3_mol"]}
            )
        np.testing.assert_allclose(
            pressure(coefficients, arrays(rows)), [r["P"] for r in rows], atol=1e-8
        )


def test_case_artifacts_and_historical_evidence_are_preserved():
    report = json.loads((DATA / "full-sensitivity-refits.json").read_text())
    for name, digest in report["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    for name in [
        "artifact-manifest.json",
        "followup-artifact-manifest.json",
        "expanded-artifact-manifest.json",
        "normalized-artifact-manifest.json",
    ]:
        for item in json.loads((DATA / name).read_text())["files"]:
            assert (
                hashlib.sha256((DATA / item["file"]).read_bytes()).hexdigest()
                == item["sha256"]
            )
    with (DATA / "full-sensitivity-residuals.csv").open(newline="") as stream:
        residuals = list(csv.DictReader(stream))
    assert Counter(r["case"] for r in residuals) == {c["id"]: 219 for c in scenarios()}
    by_id = {c["id"]: c for c in report["cases"]}
    for name in by_id:
        group = [r for r in residuals if r["case"] == name]
        objective = sum(
            (float(r["pressure_residual_gpa"]) / float(r["effective_error_gpa"])) ** 2
            for r in group
        )
        np.testing.assert_allclose(
            objective, by_id[name]["fits"][0]["objective_sum_squares"]
        )
    with (DATA / "full-baseline-native-comparisons.csv").open(newline="") as stream:
        native = list(csv.DictReader(stream))
    assert len(native) == 219
    for name, summary in report["baseline_native_residuals"].items():
        group = [r for r in native if r["dataset"] == name]
        for curve in ["source", "fit"]:
            np.testing.assert_allclose(
                np.sqrt(
                    np.mean(
                        [(float(r["observed"]) - float(r[curve])) ** 2 for r in group]
                    )
                ),
                summary[curve]["rms"],
            )
