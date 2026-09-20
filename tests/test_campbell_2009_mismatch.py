"""Independent invariants for the Campbell mismatch diagnostics."""

import json

import numpy as np
import pytest

from scripts.diagnose_campbell_2009_mismatch import (
    OUTPUT,
    fit,
    thermal_variant,
)
from scripts.reconstruct_seagle_2008_fcc import CELL_TO_MOLAR, V0_1400, fcc_pressure
from scripts.reproduce_campbell_2009_buffers import SOURCE, pressure


def test_thermal_baseline_matches_original_and_all_anchors_preserve_reference_isotherm():
    cells = np.array([38.704, 42.87, 43.317, 45.0])
    temperatures = np.array([2570.0, 1400.0, 1900.0, 1800.0])
    for delta in [6, 6.5, 7]:
        assert thermal_variant(cells, temperatures, delta=delta) == pytest.approx(
            fcc_pressure(cells, temperatures, delta), abs=1e-8
        )
    r = cells * CELL_TO_MOLAR / V0_1400
    bm3 = 180 * (r ** (-7 / 3) - r ** (-5 / 3)) * (1 + 0.75 * (r ** (-2 / 3) - 1))
    for anchor in [300, 1400, "observed"]:
        # At 1400 K no thermal remapping is needed, whatever alpha convention.
        assert thermal_variant(cells, 1400, anchor) == pytest.approx(bm3, abs=1e-8)


def test_diagnostic_fit_recovers_synthetic_coefficients_with_and_without_fixed_q():
    # Independent range of temperature and volume avoids the pressure/temperature
    # confounding in the experimental rows. Missing errors do not affect equal weights.
    source = np.array(SOURCE["feo"])
    v, t = np.meshgrid(np.linspace(9, 12, 6), np.linspace(1000, 2400, 5))
    v, t = v.ravel(), t.ravel()
    data = np.array(
        [
            v,
            t,
            pressure(v, t, source),
            np.full(v.size, np.nan),
            np.full(v.size, np.nan),
            np.full(v.size, np.nan),
        ]
    )
    for q in [None, source[5]]:
        result = fit("feo", data, fixed_q=q)
        assert result["solver_success"]
        assert result["rmse_gpa"] < 1e-9
        assert list(result["parameters"].values()) == pytest.approx(source[[1, 4, 5]])


def test_diagnostic_partitions_account_for_every_observation_and_residual():
    report = json.loads(OUTPUT.read_text())
    entries = [*report["materials"].values(), report["feo_dewaele_control"]]
    for entry in entries:
        rows = entry["rows"]
        ids = {r["row_id"] for r in rows}
        total_sse = sum(r["published_residual_gpa"] ** 2 for r in rows)
        for partition in ["study_and_marker", "temperature"]:
            groups = list(entry["groups"][partition].values())
            seen = [key for group in groups for key in group["row_ids"]]
            assert len(seen) == len(set(seen)) == len(ids)
            assert set(seen) == ids
            assert sum(g["published"]["sse_gpa2"] for g in groups) == pytest.approx(
                total_sse
            )
            assert sum(g["published_sse_fraction"] for g in groups) == pytest.approx(1)
        for groups in entry["groups"].values():
            for group in groups.values():
                assert group["without_group"]["observations"] + len(
                    group["row_ids"]
                ) == len(ids)
                assert group["without_group"]["solver_success"]
        for profile in entry["q_profiles"].values():
            for fit_result in profile["fixed_q_fits"]:
                assert fit_result["solver_success"]
                assert (
                    fit_result["sse_gpa2"] >= profile["free_q_fit"]["sse_gpa2"] - 1e-6
                )
    for variant in report["thermal_variants"].values():
        assert len(variant["checkpoint_differences_gpa"]) == 8
        assert all(f["solver_success"] for f in variant["fits"].values())
