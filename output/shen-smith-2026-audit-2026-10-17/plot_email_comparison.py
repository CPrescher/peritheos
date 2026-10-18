"""Export the two-phase email figure from the independently checked audit."""

import csv
import hashlib
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
AUDIT = Path(__file__).with_name("results.json")
OUT = ROOT / "output/pdf"
OUT.mkdir(parents=True, exist_ok=True)
STEM = OUT / "shen-smith-2026-mgo-nacl-b2-comparison"
source = ROOT / "peritheos/data/datasets/shen-smith-2026-table-s1-volumes.csv"
audit = json.loads(AUDIT.read_text())
assert hashlib.sha256(source.read_bytes()).hexdigest() == audit["source_csv_sha256"]
rows = list(csv.DictReader(source.open()))
cu = {(r["experiment"], r["run_number"]): r for r in rows if r["phase"] == "Cu"}


def copper(v):
    from scipy.constants import Avogadro

    v0 = 4 * 63.546 / (8.939 * Avogadro) * 1e24
    y = 1 - (np.asarray(v) / v0) ** (1 / 3)
    return 3 * 133.6 * y / (1 - y) ** 2 * np.exp(6.29 * y + 2.06 * y**2 + 1.65 * y**3)


def vinet(v, v0, k0, kp):
    x = (np.asarray(v) / v0) ** (1 / 3)
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def invert(p, v0, parameters):
    return np.array(
        [
            brentq(
                lambda v: (
                    vinet(v, v0, parameters["K0_gpa"], parameters["K0_prime"])
                    - pressure
                ),
                0.1 * v0,
                v0,
                xtol=1e-12,
            )
            for pressure in p
        ]
    )


plt.rcParams.update(
    {
        "font.family": "DejaVu Sans",
        "font.size": 10,
        "axes.titlesize": 11,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.edgecolor": "#4e5660",
        "axes.linewidth": 0.8,
        "xtick.labelsize": 9,
        "ytick.labelsize": 9,
        "pdf.fonttype": 42,
        "ps.fonttype": 42,
        "svg.hashsalt": "shen-smith-2026-mgo-nacl-b2-comparison",
        "savefig.facecolor": "white",
    }
)
PUBLISHED = "#b74330"
REFIT = "#176ca4"
DATA = "#62676d"
fig, axes = plt.subplots(
    2, 2, figsize=(10.6, 8.5), sharex="col", gridspec_kw={"height_ratios": [1.6, 1]}
)
fig.subplots_adjust(
    left=0.085, right=0.98, bottom=0.205, top=0.82, wspace=0.24, hspace=0.31
)
fig.suptitle(
    "MgO and NaCl-B2: published EOS and independent refits", y=0.984, fontsize=15
)
fig.text(
    0.5,
    0.947,
    "Shen and Smith (2026), Table II and Supplemental Table S1",
    ha="center",
    fontsize=11,
)
handles = [
    Line2D(
        [],
        [],
        marker="o",
        linestyle="none",
        color=DATA,
        markersize=4,
        label="Measured volume",
    ),
    Line2D([], [], color=PUBLISHED, linewidth=2, label="Published EOS (Table II)"),
    Line2D(
        [], [], color=REFIT, linewidth=2, linestyle="--", label="Our diagnostic refit"
    ),
    Line2D(
        [],
        [],
        marker="x",
        linestyle="none",
        color=DATA,
        markersize=5,
        label="Excluded B2 rows",
    ),
]
fig.legend(
    handles=handles,
    loc="upper center",
    bbox_to_anchor=(0.52, 0.923),
    ncol=4,
    frameon=False,
    fontsize=9,
    handlelength=2.3,
)
summary = []
for col, phase in enumerate(["MgO", "NaCl-B2"]):
    checked = next(r for r in audit["records"] if r["phase"] == phase)
    raw = [r for r in rows if r["phase"] == phase and r["experiment"] == "DAC-2"]
    p = np.array(
        [
            copper(float(cu[(r["experiment"], r["run_number"])]["unit_cell_volume_a3"]))
            for r in raw
        ]
    )
    v = np.array([float(r["unit_cell_volume_a3"]) for r in raw])
    sigma = np.array([float(r["unit_cell_volume_standard_error_a3"]) for r in raw])
    selected = np.ones(len(raw), dtype=bool) if phase == "MgO" else p > 35
    assert int(selected.sum()) == checked["observations"]
    v0 = checked["fixed_V0_a3"]
    published, refit = checked["published"], checked["refit"]
    rp = vinet(v, v0, published["K0_gpa"], published["K0_prime"]) - p
    rr = vinet(v, v0, refit["K0_gpa"], refit["K0_prime"]) - p
    rms_p = float(np.sqrt(np.mean(rp[selected] ** 2)))
    rms_r = float(np.sqrt(np.mean(rr[selected] ** 2)))
    assert abs(rms_p - checked["published_residuals"]["rms_gpa"]) < 1e-9
    assert abs(rms_r - checked["refit_residuals"]["rms_gpa"]) < 1e-9
    order = np.argsort(p)
    p, v, sigma, selected, rp, rr = [x[order] for x in (p, v, sigma, selected, rp, rr)]
    upper, lower = axes[:, col]
    pgrid = np.linspace(float(p[selected].min()), float(p[selected].max()), 400)
    upper.errorbar(
        p[selected],
        v[selected],
        yerr=sigma[selected],
        fmt="o",
        markersize=2.6,
        markeredgewidth=0,
        color=DATA,
        ecolor="#9da1a6",
        elinewidth=0.6,
        alpha=0.7,
        zorder=2,
    )
    upper.plot(
        pgrid, invert(pgrid, v0, published), color=PUBLISHED, linewidth=2, zorder=3
    )
    upper.plot(
        pgrid,
        invert(pgrid, v0, refit),
        color=REFIT,
        linewidth=2,
        linestyle="--",
        zorder=4,
    )
    upper.set_title(f"({chr(97 + col)}) {phase}", loc="left", fontweight="bold", pad=27)
    upper.text(
        0,
        1.025,
        f"{selected.sum()} fitted observations; fixed $V_0$ = {v0:g} "
        + r"$\mathrm{\AA}^3$",
        transform=upper.transAxes,
        fontsize=9.5,
    )
    upper.set_ylabel(r"Unit-cell volume ($\mathrm{\AA}^3$)")
    upper.text(
        0.97,
        0.97,
        f"Published: $K_0$ = {published['K0_gpa']:g} GPa, $K'_0$ = {published['K0_prime']:g}\nRefit: $K_0$ = {refit['K0_gpa']:.3f} GPa, $K'_0$ = {refit['K0_prime']:.5f}",
        transform=upper.transAxes,
        ha="right",
        va="top",
        fontsize=9,
        bbox={"facecolor": "white", "edgecolor": "none", "alpha": 0.9, "pad": 3},
    )
    lower.scatter(p[selected], rp[selected], s=10, color=PUBLISHED, alpha=0.8, zorder=3)
    lower.scatter(
        p[selected], rr[selected], s=10, marker="^", color=REFIT, alpha=0.8, zorder=3
    )
    lower.axhline(0, color="#4e5660", linewidth=0.9, linestyle=":", zorder=1)
    lower.set_title(
        f"RMS: published {rms_p:.3f} GPa; refit {rms_r:.3f} GPa",
        loc="left",
        pad=9,
        fontsize=9.5,
    )
    lower.set_ylabel(r"$P_{\mathrm{EOS}}(V_{\mathrm{obs}})-P_{\mathrm{Cu}}$ (GPa)")
    lower.set_xlabel("Cu-derived pressure (GPa)")
    for ax in (upper, lower):
        ax.grid(alpha=0.15, linewidth=0.6)
        ax.set_axisbelow(True)
        ax.tick_params(direction="out", length=3.5)
    lower.set_ylim(
        np.floor((min(rp.min(), rr.min()) - 0.15) * 2) / 2,
        np.ceil((max(rp.max(), rr.max()) + 0.15) * 2) / 2,
    )
    if phase == "NaCl-B2":
        upper.errorbar(
            p[~selected],
            v[~selected],
            yerr=sigma[~selected],
            fmt="x",
            markersize=4,
            color=DATA,
            elinewidth=0.6,
            zorder=2,
        )
        lower.scatter(
            p[~selected],
            rp[~selected],
            marker="x",
            s=24,
            linewidths=0.8,
            color=PUBLISHED,
        )
        lower.scatter(
            p[~selected], rr[~selected], marker="x", s=24, linewidths=0.8, color=REFIT
        )
        for ax in (upper, lower):
            ax.axvspan(30, 35, color="#777777", alpha=0.08, zorder=0)
        upper.set_xlim(30, 112)
        upper.set_ylim(19.7, 28.9)
    else:
        upper.set_xlim(0, 112)
        upper.set_ylim(float(np.min(v - sigma)) - 0.5, 76)
    summary.append(
        dict(
            phase=phase,
            shown=len(raw),
            fitted=int(selected.sum()),
            excluded=int((~selected).sum()),
            fixed_V0_a3=v0,
            published=published,
            refit=refit,
            published_rms_gpa=rms_p,
            refit_rms_gpa=rms_r,
        )
    )

caption = [
    "Cu pressure: Fratanduono et al. (2020), 298 K EOS; each volume paired with Cu from the same DAC and run.",
    "Refits: ordinary Vinet, unweighted pressure residuals, fixed published V0. Error bars: reported sample-volume standard errors only.",
    "Cu-derived pressures are held fixed in fitting. Residual RMS uses fitted rows; crosses mark six NaCl-B2 rows at P(Cu) <= 35 GPa.",
    "Fits are shown within the fitted pressure intervals. Diagnostic refits; original author weights and exact fit inputs remain unconfirmed.",
    "Source: Shen and Smith, Phys. Rev. B 113, 144113 (2026), doi:10.1103/fxgq-96sg; Cu: Phys. Rev. Lett. 124, 015701 (2020).",
]
fig.text(
    0.085,
    0.124,
    "\n".join(caption),
    va="top",
    fontsize=8,
    linespacing=1.55,
    color="#41464d",
)
fig.savefig(
    STEM.with_suffix(".pdf"),
    metadata={
        "Title": "MgO and NaCl-B2: published EOS and independent refits",
        "Subject": "Comparison of Shen and Smith (2026) Table II Vinet EOSs with Supplemental Table S1 observations",
        "Author": "Clemens Prescher",
    },
)
fig.savefig(STEM.with_suffix(".png"), dpi=300)
fig.savefig(STEM.with_suffix(".svg"))
STEM.with_suffix(".json").write_text(
    json.dumps(
        {
            "source_csv_sha256": audit["source_csv_sha256"],
            "residual_convention": "P_EOS(V_observed) - P_Cu",
            "records": summary,
        },
        indent=2,
    )
    + "\n"
)
print(json.dumps(summary, indent=2))
