"""Source normalization, covariance bounds and the qualified 193-row replay."""

import csv
import hashlib
import json
from collections import Counter

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, GAS_CONSTANT
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_normalized import (
    PLAN,
    conrad_width_interval,
    prepared_rows,
)
from scripts.refit_grossular_angel_2022_subset import arrays, pressure, solve


def test_all_73_channels_have_an_evidence_based_disposition():
    audit = json.loads((DATA / "remaining-channel-audit.json").read_text())
    assert sum(g["nominal_channels"] for g in audit["groups"]) == 73
    outcome = audit["outcome"]
    assert outcome["newly_usable_diagnostic_rows"] == 46
    assert outcome["remaining_unweighted_rows"] == 26
    assert outcome["normalization_anchor_not_independent"] == 1
    assert not outcome["full_author_input_recovered"]
    # A new evidence receipt upgrades verification; historical bytes stay intact.
    for filename in [
        "artifact-manifest.json",
        "followup-artifact-manifest.json",
        "expanded-artifact-manifest.json",
    ]:
        for item in json.loads((DATA / filename).read_text())["files"]:
            assert (
                hashlib.sha256((DATA / item["file"]).read_bytes()).hexdigest()
                == item["sha256"]
            )


def test_plot_normalization_is_traceable_and_does_not_repair_primary_values():
    receipt = json.loads((DATA / "normalization-recovery.json").read_text())
    with (DATA / receipt["csv"]).open(newline="") as stream:
        points = list(csv.DictReader(stream))
    assert len(points) == 77
    for name, summary in receipt["groups"].items():
        assert (
            hashlib.sha256((DATA / summary["source_table"]).read_bytes()).hexdigest()
            == summary["source_table_sha256"]
        )
        group = [r for r in points if r["dataset"] == name]
        denominator = summary["inferred_normalization_denominator"]
        residual = np.array(
            [
                float(r["digitized_molar_volume_cm3_mol"])
                - 125.35 * float(r["original_volume"]) / denominator
                for r in group
            ]
        )
        np.testing.assert_allclose(
            np.sqrt(np.mean(residual**2)),
            summary["profiled_plot_residual_cm3_mol"]["rms"],
        )
        if name != "rodehorst":
            assert abs(residual).max() < 0.003
    bose = receipt["groups"]["bosenick"]
    assert bose["source_rows_plotted"] == list(range(2, 21))
    assert abs(bose["inferred_normalization_denominator"] - 125.293) < 0.003
    rode = [r for r in points if r["dataset"] == "rodehorst"]
    for r in rode:
        residual = (
            float(r["digitized_molar_volume_cm3_mol"])
            - 125.35 * float(r["original_volume"]) / 11.850**3
        )
        assert (
            abs(residual) < 0.003
            if float(r["original_temperature_k"]) != 156
            else 0.03 < residual < 0.04
        )
    source156 = next(
        r for r in prepared_rows() if r["dataset"] == "rodehorst_vt" and r["T"] == 156
    )
    np.testing.assert_allclose(source156["Y"], 125.35 * (11.842 / 11.850) ** 3)


def test_new_selection_widths_and_conrad_covariance_bounds():
    rows = prepared_rows()
    assert len(rows) == 193
    assert Counter(r["dataset"] for r in rows) == {
        "milani_vt": 35,
        "du_vt": 13,
        "milani_pv": 7,
        "kono_ks": 36,
        "isaak_ks": 3,
        "bass_ks": 1,
        "gwanmesia_ks": 52,
        "bosenick_vt": 19,
        "rodehorst_vt": 9,
        "isaak_vt": 15,
        "conrad_ks": 3,
    }
    for r in rows:
        if not r["KS"]:
            assert r["scale"] is None
        if r["dataset"] == "kono_ks":
            assert r["scale"] == 4
        if r["dataset"] == "gwanmesia_ks":
            assert r["scale"] == 5
        if r["dataset"] == "bosenick_vt":
            assert r["sT"] == (3 if r["T"] >= 100 else 5)
        if r["dataset"] == "isaak_vt":
            assert r["sT"] == 2 and r["sY"] == 0.01
        if r["dataset"] == "rodehorst_vt":
            assert r["T"] != 298 and r["sT"] is None
    for s11, s12 in [(0.03, 0.03), (0.05, 0.11), (0.04, 0.05)]:
        lower, upper = conrad_width_interval(s11, s12)
        for rho in [-1, 0, 1]:
            width = np.sqrt(s11**2 + 4 * s12**2 + 4 * rho * s11 * s12) / 3
            assert lower - 1e-12 <= width <= upper + 1e-12
    c = [r for r in rows if r["dataset"] == "conrad_ks"]
    np.testing.assert_allclose([r["sY"] for r in c], [0.03, 0.09, 0.14 / 3])
    assert all(r["sP"] == 0.1 and r["T"] == 298 for r in c)


def test_normalized_fit_replay_independent_properties_and_scope():
    plan = json.loads(PLAN.read_text())
    report = json.loads((DATA / "normalized-refits.json").read_text())
    assert report["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert not report["author_workflow_reproduced"]
    assert not report["full_220_channel_fit"] and not report["physical_validation"]
    assert report["mixed_composition_isaak_expansion_included"]
    for name, digest in report["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    a, b = report["fits"]
    for fit in (a, b):
        assert fit["numerically_converged"] and not fit["at_parameter_bound"]
        assert fit["degrees_of_freedom"] == 187
        assert np.linalg.eigvalsh(fit["conditional_unscaled_covariance"]).min() > 0
    np.testing.assert_allclose(
        list(a["parameters"].values()), list(b["parameters"].values()), rtol=1e-6
    )
    replay = solve(
        prepared_rows(),
        plan["starts"][1],
        names=plan["free_parameters"],
        scaling=np.array(plan["parameter_scaling"]),
        lower=np.array(plan["bounds"]["lower"]),
        upper=np.array(plan["bounds"]["upper"]),
    )
    np.testing.assert_allclose(
        list(replay["parameters"].values()), list(a["parameters"].values()), rtol=1e-6
    )
    p, _ = source("grossular")
    c = list(a["parameters"].values())
    fitted = p | dict(zip([2, 3, 11, 18], c[:4]))
    for parameters, coefficients in [
        (p, [p[2], p[3], p[11], p[18], p[53], p[54]]),
        (fitted, c),
    ]:
        cases = []
        for r in prepared_rows():
            prop = properties(parameters, r["P"], r["T"], GAS_CONSTANT)
            multiplier = coefficients[r["scale"]] if r["scale"] is not None else 1
            y = multiplier * prop["ks_gpa" if r["KS"] else "volume_cm3_mol"]
            cases.append(r | {"Y": y})
        np.testing.assert_allclose(
            pressure(coefficients, arrays(cases)), [r["P"] for r in cases], atol=1e-8
        )
    with (DATA / report["native_comparison_csv"]).open(newline="") as stream:
        native = list(csv.DictReader(stream))
    assert len(native) == 193
    for name, summary in report["native_residuals"].items():
        group = [r for r in native if r["dataset"] == name]
        for curve in ("source", "fit"):
            residual = [float(r["observed"]) - float(r[curve]) for r in group]
            np.testing.assert_allclose(
                np.sqrt(np.mean(np.square(residual))), summary[curve]["rms"]
            )
