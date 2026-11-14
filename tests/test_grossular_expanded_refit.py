"""The Gwanmesia extension preserves selection, equations and diagnostic scope."""

import csv
import hashlib
import json
from collections import Counter

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, GAS_CONSTANT
from scripts.audit_pyrope_angel_2022_observations import properties
from scripts.refit_grossular_angel_2022_expanded import PLAN, prepared_rows
from scripts.refit_grossular_angel_2022_subset import (
    arrays,
    pressure,
    solve,
)
from scripts.refit_grossular_angel_2022_subset import (
    prepared_rows as previous_rows,
)


def test_extension_changes_only_declared_rows_scales_and_weight_assumptions():
    rows = prepared_rows()
    assert rows[:95] == previous_rows()
    assert Counter(r["dataset"] for r in rows) == {
        "milani_vt": 35,
        "milani_pv": 7,
        "du_vt": 13,
        "kono_ks": 36,
        "isaak_ks": 3,
        "bass_ks": 1,
        "gwanmesia_ks": 52,
    }
    assert rows[95]["T"] == 303.15 and rows[-1]["T"] == 1000.15
    for r in rows[95:]:
        assert r["KS"] and r["scale"] == 7 and r["sT"] is None
        assert r["sY"] == 0.008 * r["Y"] and r["sP"] == 0.007 * r["P"]
    plan = json.loads(PLAN.read_text())
    assert "not recovered independent one-sigma" in plan["gwanmesia_coordinate_widths"]
    assert plan["selection"]["total"] == 147
    assert plan["free_parameters"][-1] == "scale_gwanmesia_ks"


def test_new_scale_and_ks_inversion_agree_with_independent_forward_quadrature():
    p, _ = source("grossular")
    report = json.loads((DATA / "expanded-refits.json").read_text())
    for c in (
        [p[2], p[3], p[11], p[18], 1, 1, p[53], p[54]],
        list(report["fits"][0]["parameters"].values()),
    ):
        pp = p | dict(zip([2, 3, 11, 18], c[:4]))
        # All 52 distinct P/T states, with the new non-unit scale at index 7.
        rows = []
        for r in prepared_rows()[95:]:
            ks = properties(pp, r["P"], r["T"], GAS_CONSTANT)["ks_gpa"]
            rows.append(r | {"Y": c[7] * ks})
        np.testing.assert_allclose(
            pressure(c, arrays(rows)), [r["P"] for r in rows], atol=1e-8
        )


def test_expanded_fit_replays_from_declared_start_without_rewriting_sources():
    plan = json.loads(PLAN.read_text())
    report = json.loads((DATA / "expanded-refits.json").read_text())
    assert report["plan_sha256"] == hashlib.sha256(PLAN.read_bytes()).hexdigest()
    assert report["channels"] == 147
    assert not report["source_coefficients_changed"]
    assert not report["catalog_record_written"]
    assert not report["author_workflow_reproduced"]
    for name, digest in report["input_sha256"].items():
        assert hashlib.sha256((DATA / name).read_bytes()).hexdigest() == digest
    assert (
        report["official_source_sha256"]
        == "279139b3835582e97e5419283271151a9b10491f7eef89e986ce4e84e9d3dca5"
    )
    a, b = report["fits"]
    for fit in (a, b):
        assert fit["numerically_converged"] and not fit["at_parameter_bound"]
        assert fit["degrees_of_freedom"] == 139
        assert len(fit["irls_trace"]) <= 30
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
    comparison = report["parameter_comparison"]
    assert comparison["K0_gpa"]["absolute_distance_decreased"]
    assert comparison["K0_prime"]["absolute_distance_decreased"]
    assert comparison["gamma0"]["absolute_distance_decreased"]
    assert not comparison["theta0_k"]["absolute_distance_decreased"]
    with (DATA / "expanded-native-comparisons.csv").open(newline="") as stream:
        native = list(csv.DictReader(stream))
    assert len(native) == 147
    for name, summary in report["native_residual_comparison"]["groups"].items():
        group = [r for r in native if r["dataset"] == name]
        for curve in ("source", "fit"):
            residuals = np.array(
                [float(r["observed"]) - float(r[curve]) for r in group]
            )
            np.testing.assert_allclose(
                np.sqrt(np.mean(residuals**2)), summary[curve]["rms"]
            )
