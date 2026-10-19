"""Native report coordinates replay stored audits without performing new fits."""

import copy
import json

import numpy as np
import pytest

from peritheos import get_material_document
from scripts.reconstruction_diagnostics import ROOT, DiagnosticComparison, ledger
from scripts.reconstruction_figures import current_curve, curve_svg
from scripts.reconstruction_plot_data import curve_comparison
from scripts.report_numbers import rounded_metadata
from scripts.report_plot_themes import dark_svg


def replay(material, rid):
    document = get_material_document(material)
    before = copy.deepcopy(document)
    record = next(r for r in document["eos_records"] if r["identifier"] == rid)
    comparison, reason = curve_comparison(document, record)
    assert isinstance(comparison, DiagnosticComparison), reason
    assert document == before
    return comparison


@pytest.mark.parametrize(
    "material,rid,counts,units,partial",
    [
        (
            "gold",
            "gold_dorfman_2012_tange_mgo_k0_fixed_vinet",
            [160, 241, 241],
            ["GPa", "GPa", "√GPa"],
            False,
        ),
        ("mgo", "mgo_li_2006_bm3_absolute_acoustic", [11, 11], ["km/s", "km/s"], False),
        ("kcl", "kcl_b2_ma_2024_bm3_mgd", [11, 8], ["GPa", "GPa"], False),
        (
            "ca_perovskite",
            "ca_perovskite_fu_2023_bm3_mgd_refit",
            [140, 34, 34, 34],
            ["GPa"] * 4,
            False,
        ),
        (
            "mgsio3_liquid",
            "mgsio3_liquid_akins_2004_adiabatic_bm3",
            [3],
            ["GPa"],
            False,
        ),
        (
            "mgsio3_post_perovskite",
            "mgsio3_post_perovskite_mosenfelder_2009_bm3_1",
            [48, 6],
            ["GPa", "GPa"],
            False,
        ),
        (
            "mgo",
            "mgo_b1_luo_2023_vinet_thermal_5",
            [5, 4, 3],
            ["GPa", "K", "km/s"],
            False,
        ),
        ("iron", "iron_brown_2000_linear_hugoniot", [37], ["km/s"], False),
        ("mgo", "mgo_b1_duffy_ahrens_1995_hugoniot_5", [4], ["km/s"], False),
        (
            "nickel_oxide",
            "nickel_oxide_noguchi_1999_linear_hugoniot_2",
            [8],
            ["km/s"],
            False,
        ),
        ("mgo", "mgo_b1_tange_2009_vinet", [61, 17, 24], ["GPa"] * 3, True),
        ("gold", "gold_zhu_2025_pvt", [10], ["GPa"], True),
        ("platinum", "platinum_zhu_2025_pvt", [17], ["GPa"], True),
    ],
)
def test_native_inputs_units_and_partial_scope_without_regression(
    monkeypatch, material, rid, counts, units, partial
):
    def no_fit(*args, **kwargs):
        raise AssertionError("Report generation must not run regression")

    for name in ("least_squares", "minimize", "curve_fit"):
        monkeypatch.setattr("scipy.optimize." + name, no_fit)
    monkeypatch.setattr("scripts.validate_primary_eos_refits.fit_rt_eos", no_fit)
    c = replay(material, rid)
    assert [len(p.x) for p in c.panels] == counts
    assert [p.unit for p in c.panels] == units
    assert c.partial_replay is partial
    if partial:
        assert c.note.startswith("PARTIAL REPLAY:")
    for p in c.panels:
        assert np.isfinite([p.rms("source"), p.rms("refit")]).all()
        for _, x, a, b in p.curves:
            assert np.shape(x) == np.shape(a) == np.shape(b)
            assert np.isfinite([x, a, b]).all()


@pytest.mark.parametrize("mode", ["fixed", "free"])
def test_dorfman_joint_objective_and_calculated_baseline_are_distinct(mode):
    c = replay("gold", f"gold_dorfman_2012_tange_mgo_k0_{mode}_vinet")
    report = json.loads(
        (ROOT / "docs/data/dorfman-2012-cocompression-refit.json").read_text(
            encoding="utf-8"
        )
    )["fits"][mode]
    baseline, pairs, normalized = c.panels
    assert "Calculated" in baseline.observed_label
    assert baseline.rms("source") == 0
    assert pairs.observed_label == "Equal-pressure constraint"
    np.testing.assert_array_equal(normalized.observed, np.zeros(241))
    assert np.sum(normalized.source**2) == pytest.approx(
        report["published_objective"], rel=1e-10
    )
    assert np.sum(normalized.refit**2) == pytest.approx(report["objective"], rel=1e-10)


def test_native_replays_match_saved_audit_statistics_and_row_scope():
    rid = "mgsio3_liquid_akins_2004_adiabatic_bm3"
    c = replay("mgsio3_liquid", rid)
    assert c.panels[0].rms("refit") == pytest.approx(ledger()[rid]["rmse_gpa"])
    rid = "kcl_b2_ma_2024_bm3_mgd"
    c = replay("kcl", rid)
    assert c.panels[0].rms("refit") == pytest.approx(ledger()[rid]["kt_rmse_gpa"])
    assert c.panels[1].rms("refit") == pytest.approx(ledger()[rid]["rmse_gpa"])
    rid = "mgo_li_2006_bm3_absolute_acoustic"
    c = replay("mgo", rid)
    for p, key in zip(
        c.panels, ("rmse_p_wave_velocity_km_s", "rmse_s_wave_velocity_km_s")
    ):
        assert np.sqrt(np.mean((p.refit[1:] - p.observed[1:]) ** 2)) == pytest.approx(
            ledger()[rid][key]
        )
    assert "excluding the ambient anchor" in c.note
    rid = "mgo_b1_luo_2023_vinet_thermal_5"
    c = replay("mgo", rid)
    chi_square = sum(np.sum(((p.refit - p.observed) / p.sigma) ** 2) for p in c.panels)
    assert chi_square == pytest.approx(ledger()[rid]["chi_square"])
    c = replay("ca_perovskite", "ca_perovskite_fu_2023_bm3_mgd_refit")
    report = json.loads(
        (ROOT / "docs/data/fu-2023-casio3-refit-audit.json").read_text(encoding="utf-8")
    )
    fit = next(
        f
        for f in report["sensitivity_fits"]
        if f["objective"] == "unweighted_absolute_gpa"
    )
    assert sum(np.sum((p.refit - p.observed) ** 2) for p in c.panels) == pytest.approx(
        fit["residual_sum_squares"]
    )


def test_every_spin_endmember_report_uses_the_same_complete_coupled_trajectory():
    from scripts.reproduce_solomatova_2016_ferropericlase import bm3_pressure
    from scripts.validate_primary_eos_refits import SOLOMATOVA_RECORD_SOURCES

    by_source = {}
    for rid, source in SOLOMATOVA_RECORD_SOURCES.items():
        if not ledger()[rid].get("parameters"):
            continue
        material = ledger()[rid]["material_identifier"]
        c = replay(material, rid)
        assert len(c.panels) == 1
        p = c.panels[0]
        assert "Coupled" in p.title
        assert len(p.x) == ledger()[rid]["observations"]
        if source in by_source:
            np.testing.assert_array_equal(p.source, by_source[source].source)
            np.testing.assert_array_equal(p.refit, by_source[source].refit)
        else:
            by_source[source] = p
        if "_hs_" in rid:
            v0, k0, kp = [
                next(
                    r["published"]
                    for r in ledger()[rid]["parameters"]
                    if r["parameter"] == n
                )
                for n in ("V0", "K0", "K0_prime")
            ]
            assert not np.allclose(p.source, bm3_pressure(p.x, v0, k0, kp))
    assert len(by_source) == 6


def test_specialized_figure_is_reproducible_and_theme_keeps_native_labels():
    c = replay("mgo", "mgo_li_2006_bm3_absolute_acoustic")
    svg = curve_svg(c)
    assert svg == curve_svg(c) and current_curve(svg, c)
    assert b"km/s" in svg and b"Density (g/cm" in svg
    dark = dark_svg(svg)
    assert b"fill: #ffffff" not in dark and b"fill: #1e2129" in dark
    assert b"km/s" in dark
    c.panels[0].refit[0] += 0.1
    assert not current_curve(svg, c)


def test_shared_plot_limits_keep_large_source_residuals_visible(monkeypatch):
    from matplotlib.figure import Figure

    c = replay("ca_perovskite", "ca_perovskite_fu_2023_bm3_mgd_refit")
    saved = []
    original = Figure.savefig

    def capture(figure, *args, **kwargs):
        original(figure, *args, **kwargs)
        saved.append(figure)

    monkeypatch.setattr(Figure, "savefig", capture)
    curve_svg(c)
    axes = np.array(saved[0].axes).reshape(2 * len(c.panels), 2)
    for row, p in enumerate(c.panels):
        for col, side in enumerate(("source", "refit")):
            lo, hi = axes[2 * row + 1, col].get_ylim()
            residuals = getattr(p, side) - p.observed
            assert lo < residuals.min() <= residuals.max() < hi
        assert axes[2 * row + 1, 0].get_ylim() == axes[2 * row + 1, 1].get_ylim()


def test_auxiliary_acoustic_parameters_follow_the_uncertainty_rounding_rule():
    row = ledger()["mgo_li_2006_bm3_absolute_acoustic"]["acoustic_parameters"][0]
    before = copy.deepcopy(row)
    display = rounded_metadata({"acoustic_parameters": [row]})["acoustic_parameters"][0]
    assert display["published"] == "163.5" and display["published_error"] == "1.1"
    assert display["refit"] == "163.5" and display["refit_error"] == "1.4"
    assert row == before


def test_plot_cache_accepts_roundoff_but_rejects_changed_values_and_labels():
    from scripts.report_plot_cache import same_inputs

    assert same_inputs([0.0, 5.0], [2e-11, 5.0 + 2e-10])
    assert not same_inputs([0.0, 5.0], [1e-6, 5.0])
    assert not same_inputs([0.0], [0.0, 0.0])
    assert not same_inputs([True], [1])
    assert not same_inputs(["5"], [5])
    assert same_inputs([None], [None])
    c = replay("mgo", "mgo_li_2006_bm3_absolute_acoustic")
    svg = curve_svg(c)
    c.panels[0].refit[0] += 1e-11
    assert current_curve(svg, c)
    c.panels[0].refit[0] += 1e-5
    assert not current_curve(svg, c)
    c.panels[0].refit[0] -= 1e-5
    c.panels[0].unit = "Changed unit"
    assert not current_curve(svg, c)
