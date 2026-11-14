"""Validate the constructed subset's equation, inversion and fit scope."""

import json

import numpy as np

from scripts.audit_garnet_angel_2022 import independent_pressure, source
from scripts.audit_pyrope_angel_2022_observations import DATA, properties
from scripts.refit_pyrope_angel_2022_subset import (
    GAS_CONSTANT,
    NAMES,
    arrays,
    cold,
    effective_sigma,
    ks_volume,
    prepared_rows,
    pressure_at_observations,
    solve,
    thermal,
)


def test_fast_equation_and_ks_inversion_against_independent_quadrature():
    p, _ = source("pyrope")
    for theta in (50, p[11], 3000):
        pp = dict(p, **{})
        pp[11] = theta
        c = np.array([pp[i] for i in (2, 3, 11, 18)] + [1.0, 1.0])
        t = np.array([103.0, 298.0, 1000.0, 1700.0])
        v = p[1] * np.array([1.003, 1.0, 0.95, 0.90])
        pth, cv = thermal(t, c)
        expected = [
            independent_pressure(pp, vi, ti, GAS_CONSTANT) for vi, ti in zip(v, t)
        ]
        np.testing.assert_allclose(
            cold(v, c)[0] + pth, expected, atol=1e-11, rtol=1e-11
        )
        ks = cold(v, c)[1] + c[3] ** 2 * cv * t * v / (1000 * p[1] ** 2)
        np.testing.assert_allclose(ks_volume(ks, t, c, cv), v, atol=1e-9, rtol=0)
    for pressure, temperature in ((0, 298), (10, 1000), (19.86, 1700)):
        c = np.array([p[i] for i in (2, 3, 11, 18)] + [1.0, 1.0])
        expected = properties(p, pressure, temperature, GAS_CONSTANT)
        d = arrays(
            [
                {
                    "T": temperature,
                    "P": pressure,
                    "Y": expected["ks_gpa"],
                    "sT": None,
                    "sP": None,
                    "sY": 1.0,
                    "kind": "KS",
                    "scale_index": None,
                }
            ]
        )
        np.testing.assert_allclose(pressure_at_observations(c, d), pressure, atol=1e-9)


def test_partial_inputs_and_effective_variances_are_explicit():
    rows = prepared_rows()
    assert len(rows) == 92
    assert {
        name: sum(r["dataset"] == name for r in rows)
        for name in ("milani_vt", "hartwig_vt", "zou_pv", "zou_ks")
    } == {"milani_vt": 33, "hartwig_vt": 12, "zou_pv": 7, "zou_ks": 40}
    assert all(r["sT"] is None and r["sP"] is None for r in rows[:33])
    assert all(r["sY"] == 0.03 for r in rows[33:45])
    assert rows[43]["Y"] == 114.121  # Source conflict preserved, not a^3 repair.
    p, _ = source("pyrope")
    c = np.array([p[i] for i in (2, 3, 11, 18)] + [1.001, 1.0003])
    assert np.all(effective_sigma(c, arrays(rows)) > 0)


def test_joint_irls_replays_and_agrees_with_official_constructed_subset_fit():
    saved = json.loads((DATA / "subset-refits.json").read_text())
    expected = saved["cases"]["joint_recovered"]["refit"]
    seed = saved["cases"]["joint_recovered"]["source_eos_scales_only"]["parameters"]
    result = solve(prepared_rows(), [seed[n] for n in NAMES], list(range(6)))
    assert result["numerically_converged"]
    assert result["scaled_jacobian_rank"] == 6
    assert not result["at_parameter_bound"]
    np.testing.assert_allclose(
        [result["parameters"][n] for n in NAMES],
        [expected["parameters"][n] for n in NAMES],
        atol=1e-4,
        rtol=1e-7,
    )
    oracle = json.loads((DATA / "subset-eosfit/checks.json").read_text())
    assert oracle["within_declared_comparison_tolerances"]
    assert oracle["n_channels"] == 92
    for n in NAMES:
        assert (
            abs(result["parameters"][n] - oracle["parameters"][n])
            < 0.01 * result["conditional_unscaled_widths"][n]
            + oracle["saved_parameter_half_last_digit"][n]
        )
    assert saved["author_workflow_reproduced"] is False
    assert (
        abs(result["parameters"]["theta0_k"] - saved["source_parameters"]["theta0_k"])
        > 250
    )
