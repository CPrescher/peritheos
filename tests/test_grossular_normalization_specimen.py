"""Controlled input changes separate normalization from the Isaak specimen."""

import csv
import hashlib
import json
import math
from collections import Counter

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, GAS_CONSTANT
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_full_sensitivity import (
    prepared_rows as full_rows,
)
from scripts.refit_grossular_angel_2022_full_sensitivity import scenarios
from scripts.refit_grossular_angel_2022_normalization_specimen import (
    BASELINE,
    DIFFRACTION,
    PLAN,
    cases,
    denominators,
    isaak_curves,
    normalization_offsets,
    prepared_case,
)
from scripts.refit_grossular_angel_2022_subset import arrays, pressure, solve


def test_factorial_design_preserves_background_and_isolates_changes():
    design = cases()
    assert len(design) == 20 and design[0]["id"] == BASELINE
    assert len({(c["normalization"], c["isaak"]) for c in design}) == 20
    baseline = full_rows(scenarios()[0])
    assert prepared_case(design[0])[0] == baseline
    for case in design:
        rows, starts, config = prepared_case(case)
        assert len(rows) == (204 if case["isaak"] == "omit_expansion" else 219)
        expected_free = (
            5
            if case["normalization"] == "free_diffraction"
            else (2 if case["normalization"] == "free_milani_du" else 0)
        )
        assert len(config["names"]) == 6 + expected_free
        assert all(len(s) == len(config["names"]) for s in starts)
        original = {(r["dataset"], r["row"]): r for r in baseline}
        for r in rows:
            old = original[(r["dataset"], r["row"])]
            # All moduli and the compression experiment are exactly unchanged.
            if r["KS"] or r["dataset"] == "milani_pv":
                assert r == old
            assert r["P"] == old["P"] and r["T"] == old["T"] and r["sT"] == old["sT"]
            if r["dataset"] in DIFFRACTION and case["normalization"] in [
                "measured_room",
                "linear_298",
            ]:
                np.testing.assert_allclose(r["Y"] / old["Y"], r["sY"] / old["sY"])
            if r["scale"] is not None and r["scale"] >= 6:
                assert r["dataset"] in DIFFRACTION and not r["KS"]
        assert (
            sum(r["dataset"] in ["isaak_ks", "isaak_ks_estimated"] for r in rows) == 22
        )


def test_derived_reference_and_calculated_endmember_curve_are_not_observations():
    d = denominators("linear_298")
    np.testing.assert_allclose(d["milani_vt"], 1663.99875)
    np.testing.assert_allclose(d["du_vt"], 1664.2350409836065)
    np.testing.assert_allclose(d["bosenick_vt"], 125.30104)
    curves = isaak_curves()
    assert len(curves) == 15
    ratio = 1 + 0.00013 * 5 / 7
    for r in curves:
        np.testing.assert_allclose(r["gr2b_ratio_298_to_293"], ratio)
        np.testing.assert_allclose(
            r["observed_gr2b_reference_298"], r["observed_gr2b_as_plotted"] / ratio
        )
    # Independent exact integration of the first linear alpha segment.
    a298 = (19.2 - 2 * (21.3 - 19.2) / 50) * 1e-6
    integral300 = 0.5 * (a298 + 19.2e-6) * 2
    np.testing.assert_allclose(
        curves[0]["calculated_endmember_from_alpha"], 125.35 * math.exp(integral300)
    )
    integral350 = integral300 + 0.5 * (19.2e-6 + 21.3e-6) * 50
    np.testing.assert_allclose(
        curves[1]["calculated_endmember_from_alpha"], 125.35 * math.exp(integral350)
    )
    report = json.loads((DATA / "normalization-specimen-refits.json").read_text())
    assert not report["calculated_endmember_curves_are_measurements"]
    assert (
        not report["physical_validation"] and not report["author_workflow_reproduced"]
    )


def test_all_cases_converge_and_replay_free_offset_calculated_curve_control():
    report = json.loads((DATA / "normalization-specimen-refits.json").read_text())
    assert report["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert report["scientific_cases"] == 20 and report["all_starts_converged"]
    assert not report["any_start_at_bound"]
    old = json.loads((DATA / "full-sensitivity-refits.json").read_text())["cases"][0][
        "fits"
    ][0]["parameters"]
    np.testing.assert_allclose(
        list(report["cases"][0]["fits"][0]["parameters"].values()),
        list(old.values()),
        rtol=1e-6,
    )
    for case in report["cases"]:
        a, b = case["fits"]
        np.testing.assert_allclose(
            list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
        )
        assert all(
            np.linalg.eigvalsh(f["conditional_unscaled_covariance"]).min() > 0
            for f in [a, b]
        )
        for f in [a, b]:
            assert f["degrees_of_freedom"] == case["rows"] - len(
                case["free_parameters"]
            )
        assert a["parameters"]["K0_prime"] > case["effects"]["K0_prime"]["published"]
    case = report["cases"][-1]
    rows, starts, config = prepared_case(case)
    replay = solve(rows, starts[1], **config)
    np.testing.assert_allclose(
        list(replay["parameters"].values()),
        list(case["fits"][0]["parameters"].values()),
        rtol=1e-6,
    )
    c = list(case["fits"][0]["parameters"].values())
    p, _ = source("grossular")
    fitted = p | dict(zip([2, 3, 11, 18], c[:4]))
    calculated = []
    for r in rows:
        prop = properties(fitted, r["P"], r["T"], GAS_CONSTANT)
        multiplier = c[r["scale"]] if r["scale"] is not None else 1
        calculated.append(
            r | {"Y": multiplier * prop["ks_gpa" if r["KS"] else "volume_cm3_mol"]}
        )
    np.testing.assert_allclose(
        pressure(c, arrays(calculated)), [r["P"] for r in calculated], atol=1e-8
    )


def test_saved_row_receipts_effects_and_historical_evidence_match():
    report = json.loads((DATA / "normalization-specimen-refits.json").read_text())
    for name, digest in report["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    for manifest in [
        "artifact-manifest.json",
        "followup-artifact-manifest.json",
        "expanded-artifact-manifest.json",
        "normalized-artifact-manifest.json",
        "full-sensitivity-artifact-manifest.json",
    ]:
        for item in json.loads((DATA / manifest).read_text())["files"]:
            assert (
                hashlib.sha256((DATA / item["file"]).read_bytes()).hexdigest()
                == item["sha256"]
            )
    with (DATA / "normalization-specimen-residuals.csv").open(newline="") as f:
        residuals = list(csv.DictReader(f))
    assert Counter(r["case"] for r in residuals) == {
        r["id"]: r["rows"] for r in report["cases"]
    }
    for case in report["cases"]:
        group = [r for r in residuals if r["case"] == case["id"]]
        objective = sum(
            (float(r["pressure_residual_gpa"]) / float(r["effective_error_gpa"])) ** 2
            for r in group
        )
        np.testing.assert_allclose(objective, case["fits"][0]["objective_sum_squares"])
    params = {r["id"]: r["fits"][0]["parameters"] for r in report["cases"]}
    baseline = params[BASELINE]
    for case in report["cases"]:
        for name, effect in case["effects"].items():
            expected = (
                params[case["id"]][name]
                - params[case["normalization"] + "__as_plotted"][name]
                - params["plotted__" + case["isaak"]][name]
                + baseline[name]
            )
            np.testing.assert_allclose(effect["two_factor_interaction"], expected)
    with (DATA / "normalization-specimen-native-comparisons.csv").open(newline="") as f:
        native = list(csv.DictReader(f))
    alpha = [
        r
        for r in native
        if r["input_kind"] == "calculated_from_estimated_endmember_alpha"
    ]
    assert len(alpha) == 5 * 15 and all(r["dataset"] == "isaak_vt" for r in alpha)


def test_free_scale_denominator_direction_and_native_export_consistency():
    report = json.loads((DATA / "normalization-specimen-refits.json").read_text())
    offset = next(
        r
        for r in normalization_offsets(report["cases"])
        if r["case"] == "free_milani_du__as_plotted" and r["dataset"] == "du_vt"
    )
    multiplier = offset["calculated_volume_multiplier"]
    original = offset["original_denominator"]
    equivalent = offset["equivalent_denominator"]
    np.testing.assert_allclose(equivalent, original * multiplier)
    np.testing.assert_allclose(
        125.35 * 1666 / original / multiplier, 125.35 * 1666 / equivalent
    )
    assert -0.81 < offset["denominator_change"] < -0.78
    with (DATA / "normalization-specimen-native-comparisons.csv").open(newline="") as f:
        native = list(csv.DictReader(f))
    for case in report["cases"]:
        group = [r for r in native if r["case"] == case["id"]]
        for dataset, summary in case["native_residuals"].items():
            subset = [r for r in group if r["dataset"] == dataset]
            for curve in ["source", "fit"]:
                residual = [float(r["observed"]) - float(r[curve]) for r in subset]
                np.testing.assert_allclose(
                    np.sqrt(np.mean(np.square(residual))), summary[curve]["rms"]
                )
