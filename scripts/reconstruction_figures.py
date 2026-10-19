"""Deterministic P–V and pressure-residual figures for recorded reconstructions."""

from __future__ import annotations

import hashlib
import inspect
import io
from dataclasses import asdict

import numpy as np

from scripts.diagnostic_figures import diagnostic_signature, diagnostic_svg
from scripts.reconstruction_diagnostics import DiagnosticComparison
from scripts.reconstruction_plot_data import CurveComparison
from scripts.report_numbers import number
from scripts.report_plot_cache import current_svg, finish_svg


def plot_arrays(comparison: CurveComparison):
    """Isotherms are model evaluations, never lines joining mixed-T observations."""
    c = comparison
    temperatures = (
        np.array([c.source.reference_temperature])
        if c.series.temperature is None
        else np.unique(c.series.temperature)
    )
    if temperatures.size > 5:
        temperatures = temperatures[[0, len(temperatures) // 2, -1]]
    volume = np.linspace(float(min(c.series.volume)), float(max(c.series.volume)), 160)
    curves = []
    for model in (c.source, c.refit):
        curves.append(
            [
                np.asarray(model.pressure(volume, temperature), dtype=float)
                for temperature in temperatures
            ]
        )
    return volume, temperatures, curves


def plot_inputs(c: CurveComparison | DiagnosticComparison):
    if isinstance(c, DiagnosticComparison):
        return asdict(c)
    volume, temperatures, curves = plot_arrays(c)
    arrays = [
        c.series.pressure,
        c.series.volume,
        c.series.temperature,
        c.series.pressure_sigma,
        c.series.volume_sigma,
        c.source_pressure,
        c.refit_pressure,
        volume,
        temperatures,
        *curves,
    ]
    return [arrays, c.volume_unit, c.independent_refit, c.source_rmse, c.refit_rmse]


def signature(c: CurveComparison | DiagnosticComparison) -> str:
    if isinstance(c, DiagnosticComparison):
        return diagnostic_signature(c)
    payload = (
        inspect.getsource(curve_svg)
        + inspect.getsource(plot_arrays)
        + inspect.getsource(number)
    )
    return hashlib.sha256(payload.encode()).hexdigest()


def current_curve(payload: bytes, c: CurveComparison | DiagnosticComparison) -> bool:
    return current_svg(payload, plot_inputs(c), signature(c))


def curve_svg(c: CurveComparison | DiagnosticComparison) -> bytes:
    if isinstance(c, DiagnosticComparison):
        return diagnostic_svg(c)
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    volume, temperatures, curves = plot_arrays(c)
    style = {
        "svg.hashsalt": "peritheos-curves-v1",
        "svg.fonttype": "none",
        "font.family": "DejaVu Sans",
        "font.size": 10,
    }
    with plt.rc_context(style):
        fig, axes = plt.subplots(
            2,
            2,
            figsize=(10, 7),
            sharex="col",
            sharey="row",
            height_ratios=[2, 1],
            layout="constrained",
        )
        titles = (
            "Stored Peritheos refit (target)"
            if c.independent_refit
            else "Original / source fit",
            "Reconstruction / refit",
        )
        for index, (predicted, rmse, color) in enumerate(
            zip(
                (c.source_pressure, c.refit_pressure),
                (c.source_rmse, c.refit_rmse),
                ("#245b8a", "#ae4c14"),
            )
        ):
            top, residual = axes[:, index]
            top.errorbar(
                c.series.pressure,
                c.series.volume,
                xerr=c.series.pressure_sigma,
                yerr=c.series.volume_sigma,
                fmt="o",
                ms=3,
                color="#4e5660",
                ecolor="#aab2bc",
                elinewidth=0.6,
                alpha=0.8,
                label="Audit input rows",
            )
            for pressure, temperature, linestyle in zip(
                curves[index], temperatures, ("-", "--", ":", "-.", (0, (5, 2)))
            ):
                top.plot(
                    pressure,
                    volume,
                    color=color,
                    ls=linestyle,
                    lw=1.5,
                    label=f"{number(float(temperature), 4)} K",
                )
            top.set_title(titles[index], loc="left", fontweight="bold")
            top.legend(fontsize=8, frameon=False)
            residual.scatter(
                c.series.pressure, predicted - c.series.pressure, s=13, color=color
            )
            residual.axhline(0, color="#4e5660", lw=0.7)
            residual.set_title(
                f"Pressure RMS error = {number(rmse, 3)} GPa", loc="left", fontsize=10
            )
            residual.set_xlabel("Observed / calibrated pressure (GPa)")
            for ax in (top, residual):
                ax.grid(color="#e7eaee", linewidth=0.6)
                ax.spines[["top", "right"]].set_visible(False)
                ax.ticklabel_format(axis="y", useOffset=False)
        unit = c.volume_unit.replace("angstrom^3", "Å³").replace("_", " ")
        axes[0, 0].set_ylabel(f"Volume ({unit})")
        axes[1, 0].set_ylabel("P(model) − P(input) (GPa)")
        fig.suptitle(
            "Stored curves and residuals on the same audit input rows", fontsize=12
        )
        buffer = io.BytesIO()
        fig.savefig(
            buffer,
            format="svg",
            metadata={
                "Date": None,
                "Creator": "Peritheos reconstruction reports",
                "Description": "Stored coefficients; no new optimization. Pressure residuals evaluated at observed volume and temperature. RMS is descriptive, not parameter uncertainty.",
            },
        )
        plt.close(fig)
    payload = buffer.getvalue().replace(
        b"'DejaVu Sans'", b"'DejaVu Sans', Arial, sans-serif"
    )
    return finish_svg(payload, plot_inputs(c), signature(c))
