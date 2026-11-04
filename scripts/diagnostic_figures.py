"""Native-coordinate curves and residuals for specialized reconstruction audits."""

from __future__ import annotations

import hashlib
import inspect
import io
import textwrap
from dataclasses import asdict

from scripts.reconstruction_diagnostics import DiagnosticPanel
from scripts.report_numbers import number
from scripts.report_plot_cache import finish_svg


def diagnostic_signature(comparison):
    return hashlib.sha256(
        (
            inspect.getsource(diagnostic_svg)
            + inspect.getsource(DiagnosticPanel.rms)
            + inspect.getsource(number)
        ).encode()
    ).hexdigest()


def diagnostic_svg(c):
    import matplotlib

    matplotlib.use("Agg")
    from matplotlib import pyplot as plt

    with plt.rc_context(
        {
            "svg.hashsalt": "peritheos-diagnostic-v1",
            "svg.fonttype": "none",
            "font.family": "DejaVu Sans",
            "font.size": 10,
        }
    ):
        fig, axes = plt.subplots(
            2 * len(c.panels),
            2,
            figsize=(10, 5 * len(c.panels)),
            squeeze=False,
            height_ratios=[2, 1] * len(c.panels),
            layout="constrained",
        )
        # Link axes before drawing so autoscaling includes both fits. Linking
        # afterward can inherit just the refit's limits and hide source residuals.
        for row in range(len(c.panels)):
            top_left = axes[2 * row, 0]
            for offset, col in ((0, 1), (1, 0), (1, 1)):
                axes[2 * row + offset, col].sharex(top_left)
            axes[2 * row, 1].sharey(top_left)
            axes[2 * row + 1, 1].sharey(axes[2 * row + 1, 0])
        for row, p in enumerate(c.panels):
            for col, (side, color, title) in enumerate(
                zip(
                    ("source", "refit"),
                    ("#245b8a", "#ae4c14"),
                    (
                        "Stored Peritheos refit (target)"
                        if c.independent_refit
                        else c.source_title,
                        "Reconstruction / refit",
                    ),
                )
            ):
                top, residual = axes[2 * row, col], axes[2 * row + 1, col]
                prediction = getattr(p, side)
                xx, yy = (p.observed, p.x) if p.swap_axes else (p.x, p.observed)
                top.errorbar(
                    xx,
                    yy,
                    xerr=p.sigma if p.swap_axes else p.x_sigma,
                    yerr=p.x_sigma if p.swap_axes else p.sigma,
                    fmt="-" if p.observed_label.startswith("Calculated") else "o",
                    ms=3,
                    color="#4e5660",
                    ecolor="#aab2bc",
                    elinewidth=0.6,
                    alpha=0.8,
                    label=p.observed_label,
                )
                if not p.curves:
                    px, py = (prediction, p.x) if p.swap_axes else (p.x, prediction)
                    top.scatter(
                        px,
                        py,
                        s=14,
                        color=color,
                        marker="x",
                        label="Model at input coordinates",
                    )
                for line, (label, x, source, refit) in enumerate(p.curves):
                    y = source if side == "source" else refit
                    top.plot(
                        y if p.swap_axes else x,
                        x if p.swap_axes else y,
                        color=color,
                        ls=("-", "--", ":", "-.")[line % 4],
                        lw=1.4,
                        label=label,
                    )
                top.set_title(
                    textwrap.fill(p.title, width=48) + "\n" + title,
                    fontsize=10,
                    loc="left",
                )
                top.set_xlabel(p.ylabel if p.swap_axes else p.xlabel)
                top.set_ylabel(p.xlabel if p.swap_axes else p.ylabel)
                top.legend(fontsize=8, frameon=False)
                residual.scatter(
                    p.observed if p.swap_axes else p.x,
                    prediction - p.observed,
                    color=color,
                    s=13,
                )
                residual.axhline(0, color="#4e5660", lw=0.7)
                residual.set_xlabel(p.ylabel if p.swap_axes else p.xlabel)
                residual.set_ylabel(f"{p.residual_label} ({p.unit})")
                residual.set_title(
                    f"{p.statistic_label} = {number(p.rms(side), 3)} {p.unit}",
                    loc="left",
                    fontsize=10,
                )
                for ax in (top, residual):
                    ax.grid(color="#e7eaee", linewidth=0.6)
                    ax.spines[["top", "right"]].set_visible(False)
                    ax.ticklabel_format(axis="y", useOffset=False)
        buffer = io.BytesIO()
        fig.savefig(
            buffer,
            format="svg",
            metadata={
                "Date": None,
                "Creator": "Peritheos reconstruction reports",
                "Description": c.note,
            },
        )
        plt.close(fig)
    payload = buffer.getvalue().replace(
        b"'DejaVu Sans'", b"'DejaVu Sans', Arial, sans-serif"
    )
    return finish_svg(payload, asdict(c), diagnostic_signature(c))
