"""Removing source artwork must preserve numerical evidence and usable plots."""

import json
from pathlib import Path

import matplotlib.axes
import numpy as np
import pytest

from scripts import fei_curve_readouts
from scripts.audit_fei_2007_thermal_parts import reproduce as thermal_reproduce
from scripts.plot_fei_2004_2007_platinum_refits import main as plot_pt_comparison
from scripts.plot_fei_2007_figure2_validation import main as plot_pt_validation
from scripts.reproduce_fei_2007_neon import SOURCE, curve_check, plot_validation

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data"


def test_neon_replay_preserves_historical_pixels_and_independent_predictions():
    result = curve_check(json.loads(SOURCE.read_text(encoding="utf-8")))
    historical = json.loads(
        (DATA / "fei-2007-neon-reproduction.json").read_text(encoding="utf-8")
    )["curve_check"]
    assert result["image_sha256"] == historical["image_sha256"]
    assert result["axis_calibration_pixels"] == historical["axis_calibration_pixels"]
    assert len(result["points"]) == len(historical["points"]) == 18
    for actual, expected in zip(result["points"], historical["points"]):
        for key in ("temperature_k", "x_pixel", "y_pixel", "stroke_x_span"):
            assert actual[key] == expected[key]
        for key in (
            "volume_a3",
            "figure_pressure_gpa",
            "stroke_plus_one_pixel_axis_sensitivity_gpa",
        ):
            assert actual[key] == pytest.approx(expected[key], rel=0, abs=1e-10)
        for branch in ("vinet", "bm3"):
            assert actual["calculated_pressure_gpa"][branch] == pytest.approx(
                expected["calculated_pressure_gpa"][branch], abs=1e-9
            )
    assert result["replay_basis"] == (
        "preserved_numerical_pixel_readouts_source_image_not_bundled"
    )


def test_changed_pixel_readouts_are_rejected(tmp_path, monkeypatch):
    data = json.loads(fei_curve_readouts.READOUTS.read_bytes())
    data["neon_2007_figure5"]["points"][0]["x_pixel"] += 0.001
    path = tmp_path / "readouts.json"
    path.write_text(json.dumps(data), encoding="utf-8")
    monkeypatch.setattr(fei_curve_readouts, "READOUTS", path)
    with pytest.raises(ValueError, match="differ from the preserved audit"):
        fei_curve_readouts.curve_readouts("neon_2007_figure5")


def test_thermal_replay_separates_verified_inputs_from_historical_image_hashes():
    import hashlib

    result = thermal_reproduce()
    historical = json.loads((DATA / "fei-2007-thermal-comparison.json").read_bytes())
    assert result["source_images"]
    for path, metadata in result["source_images"].items():
        assert path not in result["input_sha256"]
        assert metadata["historical_sha256"] == historical["input_sha256"][path]
        assert metadata["verification"] == "not_performed_source_image_not_bundled"
    for path, checksum in result["input_sha256"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == checksum


def test_all_three_plots_render_without_source_images(tmp_path, monkeypatch):
    def reject_source_image(*args, **kwargs):
        raise AssertionError("Fei comparison plots must not embed source rasters")

    monkeypatch.setattr(matplotlib.axes.Axes, "imshow", reject_source_image)
    pt = tmp_path / "pt.png"
    comparison = tmp_path / "comparison.png"
    neon = tmp_path / "neon.png"
    plot_pt_validation(output=pt)
    result = plot_pt_comparison(output=comparison)
    historical = json.loads(
        (DATA / "fei-2004-2007-platinum-refit-curve-comparison.json").read_bytes()
    )
    for actual, expected in zip(
        result["fei2004_figure4"]["samples"], historical["fei2004_figure4"]["samples"]
    ):
        assert actual == pytest.approx(expected, rel=0, abs=1e-9)
    report = json.loads((DATA / "fei-2007-neon-reproduction.json").read_bytes())
    plot_validation(report, output=neon)
    for path in (pt, comparison, neon):
        assert path.read_bytes().startswith(b"\x89PNG\r\n\x1a\n")
        assert path.stat().st_size > 10000


def test_removed_source_rasters_are_absent_and_explicitly_excluded():
    policy = json.loads((DATA / "license-inventory-scopes.json").read_bytes())
    config = (ROOT / "pyproject.toml").read_text(encoding="utf-8")
    prefix = "docs/data/fei-2007-thermal-source-figures/"
    removed = [
        p for p in policy["local_distribution_exclusions"] if p.startswith(prefix)
    ]
    assert len(removed) == 4
    assert f'"{prefix}**"' in config
    assert all(not (ROOT / path).exists() for path in removed)
    rows = fei_curve_readouts.curve_readouts("platinum_2004_figure4")
    np.testing.assert_array_equal(
        [r["temperature_k"] for r in rows["points"]],
        [300, 1473, 1673, 1873] * 3,
    )
