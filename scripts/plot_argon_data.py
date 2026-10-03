"""Plot the collected argon tables without fitting or altering source coordinates.

Requires matplotlib and plotly. Outputs offline HTML, PNG/PDF figures, a
normalized P-V table, and a source inventory. Run from any working directory.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import subprocess
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import plotly.graph_objects as go
from plotly.offline import get_plotlyjs

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
NA = 6.02214076e23
MOLAR_TO_ATOMIC = 1e24 / NA
MASS_PER_ATOM = 39.948 * MOLAR_TO_ATOMIC
COLORS = {
    "Dewaele 2021": "#007F86",
    "Finger 1981": "#DE8F05",
    "Ross 1986": "#4774B5",
    "Ono 2020": "#9C4DA2",
    "Errandonea 2006": "#D34D43",
    "Chen 2010": "#766247",
    "Wittlinger 1997": "#222222",
    "Anderson 1972/1975": "#6B8F23",
    "Xiao 2025": "#146FB1",
    "Barker 1987": "#CA8135",
    "Grimsditch 1986": "#72949E",
    "Ross 1985 precursor": "#A9BBD4",
    "Maltby 2024": "#E377A1",
}
DOIS = {
    "Dewaele 2021": "10.1038/s41598-021-93995-y",
    "Finger 1981": "10.1063/1.92597",
    "Ross 1986": "10.1063/1.451346",
    "Ono 2020": "10.1038/s41598-020-58252-8",
    "Errandonea 2006": "10.1103/PhysRevB.73.092106",
    "Chen 2010": "10.1103/PhysRevB.81.144110",
    "Anderson 1972/1975": "10.1016/0022-3697(75)90004-9",
    "Xiao 2025": "10.1007/s10765-024-03469-2",
    "Barker 1987": "10.1063/1.452187",
    "Grimsditch 1986": "10.1103/PhysRevB.33.7192",
    "Maltby 2024": "10.1063/5.0237497",
}


def number(row, name):
    value = row.get(name, "")
    return float(value) if value != "" else None


def load():
    tables = {}
    inventory = []
    for path in sorted(DATA.glob("argon*.csv")):
        rows = list(csv.DictReader(path.open()))
        tables[path.name] = rows
        inventory.append(
            {
                "file": path.name,
                "rows": len(rows),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            }
        )
    records = []

    def append(
        file,
        source,
        role,
        volume,
        pressure="pressure_gpa",
        temperature=None,
        phase="fcc",
        pscale=1,
        vscale=1,
        skip=None,
        note="",
        sp=None,
        sv=None,
    ):
        for i, row in enumerate(tables[file], 2):
            if skip and skip(row):
                continue
            v = volume(row) if callable(volume) else number(row, volume)
            p = pressure(row) if callable(pressure) else number(row, pressure)
            if v is None:
                continue
            t = number(row, "temperature_k")
            records.append(
                dict(
                    source=source,
                    dataset=file,
                    csv_line=i,
                    role=role,
                    phase=row.get("phase", phase),
                    pressure_gpa=p * pscale if p is not None else None,
                    volume_atomic_a3=v * vscale,
                    temperature_k=t if t is not None else temperature,
                    pressure_error_gpa=sp(row) if sp else None,
                    volume_error_atomic_a3=sv(row) if sv else None,
                    fit_selected=row.get("selected_for_room_temperature_fit", "")
                    == "1",
                    notes=note
                    + ("; " + row["source_note"] if row.get("source_note") else ""),
                    original=row,
                )
            )

    append(
        "argon-fcc-dewaele-2021-supplement.csv",
        "Dewaele 2021",
        "measurement",
        lambda r: number(r, "argon_a_angstrom") ** 3 / 4,
    )
    append(
        "argon-fcc-finger-1981-table1.csv",
        "Finger 1981",
        "measurement",
        "volume_a3",
        vscale=0.25,
        sp=lambda r: number(r, "pressure_uncertainty_kbar") / 10,
        sv=lambda r: number(r, "volume_uncertainty_a3") / 4,
    )
    append(
        "argon-ross-1986-table1.csv",
        "Ross 1986",
        "measurement",
        "volume_a3",
        vscale=0.25,
        sp=lambda r: number(r, "pressure_error_gpa"),
        note="Printed molar volume retained, including anomalous source rows; pressure errors have unspecified confidence.",
    )
    append(
        "argon-fcc-ono-2020-table1.csv",
        "Ono 2020",
        "measurement",
        "volume_a3",
        vscale=0.25,
        temperature=300,
        sp=lambda r: number(r, "pressure_sigma_gpa"),
        sv=lambda r: number(r, "volume_sigma_a3") / 4,
    )
    append(
        "argon-fcc-errandonea-2006-figure5-digitized.csv",
        "Errandonea 2006",
        "digitized measurement",
        "volume_a3",
        vscale=0.25,
        note="Room temperature (numeric value not tabulated here); digitization bounds are not experimental sigma.",
    )
    append(
        "argon-fcc-chen-2010-figure5-brillouin.csv",
        "Chen 2010",
        "acoustically inferred measurement",
        "volume_a3",
        vscale=0.25,
        note="Brillouin-derived density digitized from Figure 5; graphical halfwidths retained in original row.",
    )
    append(
        "argon-hcp-wittlinger-1997-figure3-digitized.csv",
        "Wittlinger 1997",
        "digitized measurement",
        "volume_a3_conventional_cell",
        vscale=0.5,
        phase="hcp",
        note="Room temperature; hcp cell contains 2 atoms. Published fit not reproduced; points retained.",
    )
    append(
        "argon-anderson-1972-appendix-a.csv",
        "Anderson 1972/1975",
        "measurement",
        "volume_a3",
        "pressure_kbar",
        pscale=0.1,
        vscale=0.25,
        skip=lambda r: r["usable"] != "1",
        note="Raw thesis observations: actual T; small-holder length/common-temperature corrections not applied.",
    )
    append(
        "argon-fcc-xiao-2025-table5.csv",
        "Xiao 2025",
        "measurement",
        "molar_volume_cm3_mol",
        pressure=lambda r: None,
        vscale=MOLAR_TO_ATOMIC,
        sv=lambda r: (
            number(r, "molar_volume_cm3_mol")
            * MOLAR_TO_ATOMIC
            * number(r, "molar_volume_relative_sigma")
        ),
        note="Near sublimation; sample pressure not tabulated. Shown on T-V axes, never assigned a measured P=0.",
    )
    append(
        "argon-ross-1985-precursor-figure1-subset.csv",
        "Ross 1985 precursor",
        "earlier overlapping measurement",
        "volume_a3",
        vscale=0.25,
        note="Earlier subset; potentially overlaps Ross 1986. Hidden initially.",
    )
    append(
        "argon-ross-1986-table3.csv",
        "Ross 1986",
        "calculated solid isotherm",
        "volume_a3",
        vscale=0.25,
    )
    append(
        "argon-ross-1986-table2.csv",
        "Ross 1986",
        "calculated liquid Hugoniot",
        "volume_a3",
        vscale=0.25,
        phase="liquid",
    )
    append(
        "argon-fcc-chen-2010-figure5-curve.csv",
        "Chen 2010",
        "digitized published fit",
        "volume_a3",
        vscale=0.25,
    )
    for kind, role in [
        ("theoretical", "digitized theoretical solid curve"),
        ("experimental", "redigitized Ross measurement"),
    ]:
        append(
            "argon-barker-1987-figure5-selected.csv",
            "Barker 1987",
            role,
            "volume_atomic_a3",
            phase="solid (phase unspecified)",
            skip=lambda r, k=kind: r["determination_method"] != k,
            note="Experimental circles originate from Ross 1986; not new independent Barker measurements."
            if kind == "experimental"
            else "Digitized BFW+AT calculation, not independently regenerated.",
        )
    append(
        "argon-grimsditch-1986-table1.csv",
        "Grimsditch 1986",
        "adopted density (not independent P-V)",
        lambda r: MASS_PER_ATOM / number(r, "density_g_cm3"),
        note="Room temperature unspecified numerically; solid density from earlier diffraction EOS, liquid density from reference source.",
    )
    for kind in ["experimental", "theoretical"]:
        append(
            "argon-barker-1987-table1-liquid.csv",
            "Barker 1987",
            "liquid " + kind,
            "molar_volume_cm3_mol",
            pressure=lambda r: (
                number(r, "pv_over_nkt")
                * 8.31446261815324
                * number(r, "temperature_k")
                / number(r, "molar_volume_cm3_mol")
                / 1000
            ),
            vscale=MOLAR_TO_ATOMIC,
            phase="liquid",
            skip=lambda r, k=kind: r["determination_method"] != k,
            note="Pressure converted from pV/NkT using SI R; six theory alternatives share T and V.",
        )
    maltby = json.loads((ROOT / "curation/argon/maltby-2024.json").read_text())[
        "source_table8"
    ]
    records.append(
        dict(
            source="Maltby 2024",
            dataset="curation/argon/maltby-2024.json:source_table8",
            csv_line=None,
            role="calculated reference checkpoint",
            phase="fcc",
            pressure_gpa=maltby["pressure_mpa"] / 1000,
            volume_atomic_a3=maltby["volume_cm3_mol"] * MOLAR_TO_ATOMIC,
            temperature_k=70,
            pressure_error_gpa=None,
            volume_error_atomic_a3=None,
            fit_selected=False,
            notes="Published calculated checkpoint; literal model has documented reproduction discrepancy. Not experimental data.",
            original=maltby,
        )
    )
    return tables, inventory, records


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    out = args.output.resolve()
    out.mkdir(parents=True, exist_ok=True)
    tables, inventory, records = load()
    pv = [r for r in records if r["pressure_gpa"] is not None]
    roles = {
        "measurement",
        "digitized measurement",
        "acoustically inferred measurement",
    }
    measured = [r for r in pv if r["role"] in roles]
    room = [
        r
        for r in measured
        if r["temperature_k"] is None or 285 <= r["temperature_k"] <= 305
    ]
    cold = [
        r
        for r in measured
        if r["temperature_k"] is not None and r["temperature_k"] < 285
    ]
    neutron = [r for r in records if r["source"] == "Xiao 2025"]
    # Scientific invariants and conversion checks guard against common overlay errors.
    assert len([r for r in pv if r["source"] == "Dewaele 2021"]) == 280
    assert (
        sum(r["fit_selected"] for r in records if r["source"] == "Dewaele 2021") == 95
    )
    assert len([r for r in measured if r["source"] == "Anderson 1972/1975"]) == 379
    assert len(neutron) == 22 and all(r["pressure_gpa"] is None for r in neutron)
    assert (
        len([r for r in records if r["role"] == "digitized theoretical solid curve"])
        == 21
    )
    assert (
        len([r for r in records if r["role"] == "redigitized Ross measurement"]) == 11
    )
    assert all(5 < r["volume_atomic_a3"] < 70 for r in records)

    plt.rcParams.update(
        {
            "font.family": "DejaVu Sans",
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "axes.titleweight": "bold",
            "axes.labelcolor": "#344454",
            "grid.color": "#DDE3E8",
            "grid.linewidth": 0.6,
        }
    )
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    markers = ["o", "s", "^", "D", "v", "P", "X"]
    for source, marker in zip(dict.fromkeys(r["source"] for r in room), markers):
        rs = [r for r in room if r["source"] == source]
        label = f"{source}{' · hcp' if source == 'Wittlinger 1997' else ''} ({len(rs)})"
        for ax in axs[0]:
            ax.scatter(
                [r["pressure_gpa"] for r in rs],
                [r["volume_atomic_a3"] for r in rs],
                s=21 if source != "Wittlinger 1997" else 40,
                marker=marker,
                color=COLORS[source],
                alpha=0.7,
                linewidths=0.5,
                label=label,
            )
    for ax, xmax in zip(axs[0], [145, 12]):
        ax.set(
            xlim=(0, xmax),
            ylim=(11.5, 36),
            xlabel="Pressure (GPa)",
            ylabel="Volume (Å³ / atom)",
        )
        ax.grid(alpha=0.6)
    axs[0, 0].set_title("A   Room-temperature compression · full range", loc="left")
    axs[0, 1].set_title("B   Room-temperature compression · low-P detail", loc="left")
    axs[0, 1].set_ylim(21, 35)
    axs[0, 0].legend(fontsize=8, frameon=True, facecolor="white", loc="upper right")
    for source, marker in [("Anderson 1972/1975", "o"), ("Dewaele 2021", "^")]:
        rs = [r for r in cold if r["source"] == source]
        sc = axs[1, 0].scatter(
            [r["pressure_gpa"] for r in rs],
            [r["volume_atomic_a3"] for r in rs],
            c=[r["temperature_k"] for r in rs],
            cmap="viridis",
            vmin=0,
            vmax=285,
            marker=marker,
            s=20,
            alpha=0.8,
            linewidths=0,
            label=f"{source} ({len(rs)})",
        )
    axs[1, 0].set(
        xscale="log",
        xlim=(0.02, 25),
        ylim=(16, 40.5),
        xlabel="Pressure (GPa, log scale)",
        ylabel="Volume (Å³ / atom)",
    )
    axs[1, 0].set_title("C   Cryogenic compression · actual temperatures", loc="left")
    axs[1, 0].legend(fontsize=8, loc="lower left")
    fig.colorbar(sc, ax=axs[1, 0], label="Temperature (K)", fraction=0.045, pad=0.02)
    axs[1, 0].grid(alpha=0.6)
    axs[1, 1].errorbar(
        [r["temperature_k"] for r in neutron],
        [r["volume_atomic_a3"] for r in neutron],
        xerr=0.02,
        yerr=[r["volume_error_atomic_a3"] for r in neutron],
        fmt="o",
        ms=4,
        color=COLORS["Xiao 2025"],
        ecolor="#97B7D5",
        capsize=2,
        label="Xiao 2025 · 22 neutron observations",
    )
    axs[1, 1].set(xlabel="Temperature (K)", ylabel="Volume (Å³ / atom)", xlim=(5, 53))
    axs[1, 1].set_title("D   Near-sublimation thermal expansion", loc="left")
    axs[1, 1].legend(fontsize=8, loc="upper left")
    axs[1, 1].text(
        0.03,
        0.75,
        "Pressure not tabulated; no P = 0 assigned.\nError bars: source standard uncertainties (k = 1).",
        transform=axs[1, 1].transAxes,
        fontsize=8,
        color="#536270",
    )
    axs[1, 1].grid(alpha=0.6)
    fig.suptitle(
        "ARGON  |  Collected experimental data",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=20,
        weight="bold",
    )
    fig.text(
        0.075,
        0.937,
        "Reported coordinates retained • fcc volumes ÷ 4; hcp volumes ÷ 2 • no pressure-scale recalibration",
        color="#536270",
        fontsize=11,
    )
    fig.text(
        0.075,
        0.022,
        "Chen: acoustically inferred density. Anderson: raw thesis volumes; unresolved small-holder corrections.\nCalculated curves, adopted densities and overlapping Ross re-digitizations appear separately in the interactive archive.",
        fontsize=9,
        color="#536270",
    )
    fig.subplots_adjust(
        left=0.075, right=0.96, top=0.88, bottom=0.105, wspace=0.25, hspace=0.33
    )
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out / f"argon-experimental-overview.{ext}", dpi=190)
    plt.close(fig)

    # Companion figure: distinct physical quantities and explicitly calculated states.
    fig, axs = plt.subplots(2, 2, figsize=(14, 10))
    for source, role in [
        ("Ross 1986", "calculated solid isotherm"),
        ("Barker 1987", "digitized theoretical solid curve"),
        ("Chen 2010", "digitized published fit"),
    ]:
        rs = sorted(
            [
                r
                for r in pv
                if r["source"] == source and r["role"] == role and r["pressure_gpa"] > 0
            ],
            key=lambda r: r["pressure_gpa"],
        )
        axs[0, 0].plot(
            [r["pressure_gpa"] for r in rs],
            [r["volume_atomic_a3"] for r in rs],
            "--",
            color=COLORS[source],
            label=f"{source} · {role}",
        )
    axs[0, 0].scatter(
        [r["pressure_gpa"] for r in room],
        [r["volume_atomic_a3"] for r in room],
        color="#AAB2B9",
        s=8,
        alpha=0.5,
        label="Room-T observations",
    )
    axs[0, 0].set(
        xscale="log", xlabel="Pressure (GPa, log scale)", ylabel="Volume (Å³ / atom)"
    )
    axs[0, 0].legend(fontsize=7, loc="upper right")
    axs[0, 0].set_title("A   Published curves and calculated states", loc="left")
    liquid = [r for r in pv if r["phase"] == "liquid"]
    for source, role in dict.fromkeys((r["source"], r["role"]) for r in liquid):
        rs = [r for r in liquid if r["source"] == source and r["role"] == role]
        axs[0, 1].scatter(
            [r["pressure_gpa"] for r in rs],
            [r["volume_atomic_a3"] for r in rs],
            s=28,
            marker="x" if "calculated" in role or "theoretical" in role else "o",
            color=COLORS[source],
            label=f"{source} · {role}",
        )
    axs[0, 1].set_xscale("symlog", linthresh=0.03)
    axs[0, 1].set(
        xlabel="Pressure (GPa, symmetric log; linear below 0.03)",
        ylabel="Volume (Å³ / atom)",
    )
    axs[0, 1].set_title("B   Liquid states · temperatures differ", loc="left")
    axs[0, 1].legend(fontsize=7, loc="upper right")
    for phase, marker in [("liquid", "o"), ("fcc", "^")]:
        rs = [
            r for r in tables["argon-grimsditch-1986-table1.csv"] if r["phase"] == phase
        ]
        axs[1, 0].scatter(
            [float(r["pressure_gpa"]) for r in rs],
            [float(r["brillouin_shift_cm_inverse"]) for r in rs],
            s=23,
            marker=marker,
            label=phase,
        )
    axs[1, 0].set(xlabel="Pressure (GPa)", ylabel="Brillouin shift (cm⁻¹)")
    axs[1, 0].set_title("C   Grimsditch 1986 · 98 acoustic observations", loc="left")
    axs[1, 0].legend(fontsize=8)
    elast = tables["argon-fcc-grimsditch-1986-table2.csv"]
    for name, err, label in [
        ("bulk_modulus_gpa", "bulk_modulus_error_gpa", "B (adopted)"),
        ("cstar_gpa", "cstar_error_gpa", "C* (envelope)"),
        ("c11_upper_bound_gpa", "c11_error_gpa", "C11 upper bound"),
        ("c44_gpa", "c44_error_gpa", "C44 (derived)"),
        ("c12_lower_bound_gpa", "c12_error_gpa", "C12 lower bound"),
    ]:
        axs[1, 1].errorbar(
            [float(r["pressure_gpa"]) for r in elast],
            [float(r[name]) for r in elast],
            yerr=[float(r[err]) for r in elast],
            fmt="o--",
            ms=3,
            lw=0.8,
            capsize=2,
            label=label,
        )
    axs[1, 1].set(xlabel="Pressure (GPa)", ylabel="Elastic modulus or bound (GPa)")
    axs[1, 1].set_title("D   Grimsditch 1986 · elastic constraints", loc="left")
    axs[1, 1].legend(fontsize=8)
    for ax in axs.flat:
        ax.grid(alpha=0.6)
    fig.suptitle(
        "ARGON  |  Theory, liquid states and acoustic constraints",
        x=0.075,
        y=0.975,
        ha="left",
        fontsize=17,
        weight="bold",
    )
    fig.text(
        0.075,
        0.937,
        "Dashed curves are published graphical/tabulated results, not new fits or independent model reproductions.",
        color="#536270",
        fontsize=10,
    )
    fig.text(
        0.075,
        0.025,
        "Ross liquid states: calculated Hugoniot, 87–11,963 K. Barker liquid comparison: 100 K. Grimsditch: room T.\nGrimsditch density and elastic quantities are adopted/derived; they are not independent compression measurements.",
        fontsize=9,
        color="#536270",
    )
    fig.subplots_adjust(
        left=0.075, right=0.96, top=0.88, bottom=0.105, wspace=0.23, hspace=0.33
    )
    for ext in ("png", "pdf", "svg"):
        fig.savefig(out / f"argon-supporting-data.{ext}", dpi=190)
    plt.close(fig)

    def trace(rs, label, mode="markers", symbol="circle", visible=True):
        hover = []
        for r in rs:
            t = (
                f"{r['temperature_k']:g} K"
                if r["temperature_k"] is not None
                else "Room temperature; numeric T unspecified"
            )
            hover.append(
                f"{r['source']} · {r['phase']}<br>{r['role']}<br>{t}<br>{r['dataset']} : line {r['csv_line']}<br>{html.escape(r['notes'])}"
            )
        return go.Scatter(
            x=[r["pressure_gpa"] for r in rs],
            y=[r["volume_atomic_a3"] for r in rs],
            name=label,
            mode=mode,
            visible=visible,
            text=hover,
            hovertemplate="%{text}<br>P = %{x:.6g} GPa<br>V = %{y:.6g} Å³/atom<extra></extra>",
            marker={
                "size": 7,
                "symbol": symbol,
                "color": COLORS[rs[0]["source"]],
                "opacity": 0.75,
            },
            line={"color": COLORS[rs[0]["source"]], "dash": "dash", "width": 2},
        )

    graphs = []

    def graph(name, caption, traces, xlabel, ylabel, **layout):
        f = go.Figure(traces)
        f.update_layout(
            template="plotly_white",
            height=680,
            title=name,
            xaxis_title=xlabel,
            yaxis_title=ylabel,
            hovermode="closest",
            legend={"orientation": "h", "y": -0.18, "font": {"size": 11}},
            margin={"l": 75, "r": 35, "t": 60, "b": 170},
            **layout,
        )
        graphs.append((name, caption, f))

    traces = []
    for source in dict.fromkeys(r["source"] for r in measured):
        rs = [r for r in measured if r["source"] == source]
        if source == "Dewaele 2021":
            for label, group, marker in [
                ("fit subset", [r for r in rs if r["fit_selected"]], "diamond"),
                (
                    "other observations",
                    [r for r in rs if not r["fit_selected"]],
                    "circle-open",
                ),
            ]:
                traces.append(
                    trace(group, f"{source} · {label} ({len(group)})", symbol=marker)
                )
        else:
            traces.append(
                trace(
                    rs,
                    f"{source} · {'hcp' if source == 'Wittlinger 1997' else 'fcc'} ({len(rs)})",
                    symbol="x" if source == "Wittlinger 1997" else "circle",
                )
            )
    for source, role in dict.fromkeys(
        (r["source"], r["role"])
        for r in pv
        if r not in measured and r["phase"] != "liquid"
    ):
        rs = sorted(
            [
                r
                for r in pv
                if r["source"] == source
                and r["role"] == role
                and r["phase"] != "liquid"
            ],
            key=lambda r: r["pressure_gpa"],
        )
        curve = "curve" in role or "isotherm" in role or "fit" in role
        traces.append(
            trace(
                rs,
                f"{source} · {role} ({len(rs)})",
                "lines+markers" if curve else "markers",
                visible="legendonly",
            )
        )
    graph(
        "All solid P–V data",
        "All finite measured P–V coordinates are visible initially. Click legend entries to add published theory, adopted densities, or overlapping historical digitizations. These optional layers are not additional independent measurements. Hover gives source row and actual temperature.",
        traces,
        "Pressure (GPa)",
        "Volume (Å³ / atom)",
        xaxis={"range": [0, 145]},
    )
    graph(
        "Room-temperature compression",
        "Reported room-temperature data (numeric temperatures 290–300 K where available); hcp Wittlinger distinguished from fcc. No isothermal correction or recalibration applied. Zoom for low-pressure detail.",
        [
            trace(
                [r for r in room if r["source"] == source],
                f"{source} · {'hcp' if source == 'Wittlinger 1997' else 'fcc'}",
                symbol="x" if source == "Wittlinger 1997" else "circle",
            )
            for source in dict.fromkeys(r["source"] for r in room)
        ],
        "Pressure (GPa)",
        "Volume (Å³ / atom)",
    )
    ct = go.Scatter(
        x=[r["pressure_gpa"] for r in cold],
        y=[r["volume_atomic_a3"] for r in cold],
        mode="markers",
        text=[
            f"{r['source']} · line {r['csv_line']} · {r['temperature_k']} K<br>{html.escape(r['notes'])}"
            for r in cold
        ],
        marker={
            "size": 7,
            "color": [r["temperature_k"] for r in cold],
            "colorscale": "Viridis",
            "colorbar": {"title": "K"},
        },
        hovertemplate="%{text}<br>P=%{x:.6g} GPa<br>V=%{y:.6g} Å³/atom<extra></extra>",
        name="Cryogenic observations",
    )
    graph(
        "Cryogenic compression",
        "Anderson thesis raw data and Dewaele low-temperature measurements. Different temperatures are retained; no correction to a common isotherm. Log-pressure axis separates the 0.025–2 GPa thesis data from higher-pressure diffraction.",
        [ct],
        "Pressure (GPa, logarithmic)",
        "Volume (Å³ / atom)",
        xaxis={"type": "log"},
    )
    graph(
        "Near-sublimation volumes",
        "Xiao Table 5: pressure was not tabulated and is not invented. Error bars are source k=1 standard uncertainties.",
        [
            go.Scatter(
                x=[r["temperature_k"] for r in neutron],
                y=[r["volume_atomic_a3"] for r in neutron],
                mode="markers",
                name="Xiao 2025 · 22 observations",
                error_x={"array": [0.02] * 22},
                error_y={"array": [r["volume_error_atomic_a3"] for r in neutron]},
            )
        ],
        "Temperature (K)",
        "Volume (Å³ / atom)",
    )
    graph(
        "Liquid comparisons",
        "Separate from solid compression: calculated Ross Hugoniot (87–11,963 K), Barker 100 K experimental/theory comparison, and Grimsditch room-temperature adopted densities. Hover for provenance and temperature.",
        [
            trace(
                [r for r in liquid if r["source"] == source and r["role"] == role],
                f"{source} · {role}",
            )
            for source, role in dict.fromkeys((r["source"], r["role"]) for r in liquid)
        ],
        "Pressure (GPa)",
        "Volume (Å³ / atom)",
    )
    acoustic = tables["argon-grimsditch-1986-table1.csv"]
    graph(
        "Acoustic observations",
        "All 98 Grimsditch Table I rows; repeated pressures and solid anisotropic branches are retained. The source labels liquid/solid explicitly.",
        [
            go.Scatter(
                x=[float(r["pressure_gpa"]) for r in acoustic if r["phase"] == phase],
                y=[
                    float(r["brillouin_shift_cm_inverse"])
                    for r in acoustic
                    if r["phase"] == phase
                ],
                mode="markers",
                name=phase,
            )
            for phase in ["liquid", "fcc"]
        ],
        "Pressure (GPa)",
        "Brillouin shift (cm⁻¹)",
    )
    assert sum(len(t.x) for t in graphs[-1][2].data) == 98
    graph(
        "Elastic constraints",
        "Grimsditch Table II: adopted bulk modulus, acoustic envelope, derived C44, and upper/lower bounds. Printed errors have no invented confidence level.",
        [
            go.Scatter(
                x=[float(r["pressure_gpa"]) for r in elast],
                y=[float(r[name]) for r in elast],
                error_y={"array": [float(r[err]) for r in elast]},
                name=label,
                mode="markers",
            )
            for name, err, label in [
                ("bulk_modulus_gpa", "bulk_modulus_error_gpa", "B adopted"),
                ("cstar_gpa", "cstar_error_gpa", "C* envelope"),
                ("c11_upper_bound_gpa", "c11_error_gpa", "C11 upper bound"),
                ("c44_gpa", "c44_error_gpa", "C44 derived"),
                ("c12_lower_bound_gpa", "c12_error_gpa", "C12 lower bound"),
            ]
        ],
        "Pressure (GPa)",
        "Modulus or bound (GPa)",
    )

    # All original columns, including fit parameters and non-PV properties, remain inspectable.
    details = []
    for filename, rows in tables.items():
        headings = list(rows[0])
        header = "".join(f"<th>{html.escape(k)}</th>" for k in headings)
        body = "".join(
            "<tr>"
            + "".join(f"<td>{html.escape(r[k])}</td>" for k in headings)
            + "</tr>"
            for r in rows
        )
        details.append(
            f"<details><summary>{filename} · {len(rows)} rows</summary><div class='table'><table><thead><tr>{header}</tr></thead><tbody>{body}</tbody></table></div></details>"
        )
    payload = "[" + ",".join(f.to_json() for _, _, f in graphs) + "]"
    buttons = "".join(
        f"<button onclick='show({i})'>{html.escape(name)}</button>"
        for i, (name, _, _) in enumerate(graphs)
    )
    citations = " · ".join(
        f"<a href='https://doi.org/{doi}'>{name}</a>" for name, doi in DOIS.items()
    )
    page = (
        """<!doctype html><html><head><meta charset='utf-8'><title>Argon data atlas</title>
<style>body{font:16px system-ui;color:#253746;background:#f3f6f8;margin:0}main{max-width:1350px;margin:auto;padding:32px}h1{font-size:34px;margin-bottom:8px}p{line-height:1.6}.muted{color:#526574}nav{display:flex;flex-wrap:wrap;gap:8px}button{background:white;border:1px solid #b7c9d2;border-radius:6px;padding:10px 15px;cursor:pointer}button.active{background:#006e78;color:white}#plot{background:white;margin-top:12px;border-radius:8px}.table{overflow:auto;max-height:500px}td,th{padding:7px;white-space:nowrap;border-bottom:1px solid #ddd;font:12px monospace;text-align:left}th{position:sticky;top:0;background:#e6edf1}details{margin:10px 0;background:white;padding:12px}summary{cursor:pointer}a{color:#006e78}</style>
<script>"""
        + get_plotlyjs()
        + """</script></head><body><main><h1>Argon data atlas</h1>
<p class='muted'>Collected source tables • all volumes in Å³ per atom • snapshot 26 September 2026</p>
<p>Click a source in the legend to hide/show it; double-click to isolate it. Drag to zoom and double-click the plot to reset. All plotting code and data are embedded, so this file works offline.</p>
<nav>"""
        + buttons
        + "</nav><p id='caption'></p><div id='plot'></div>"
        + """
<p class='muted'>The plot retains reported pressures, temperatures, repeated rows and known source anomalies. Different pressure calibrations are not harmonized. Zero-point, thermal and sample-holder corrections are not silently applied. Error bars are shown in the dedicated Xiao/elastic views; source error and digitization columns remain in the raw tables.</p>
<h2>Complete source-table archive</h2><p>These tables also include fit coefficients, derived quantities, phase observations and repeated graphical representations. Row counts are not counts of independent experiments. Eight Dewaele rows lack pressure; 15 Anderson zero-volume sentinels are retained here but excluded from compression plots. Maltby contributes a calculated checkpoint, not a new experimental dataset.</p>
"""
        + "".join(details)
        + "<details><summary>Maltby 2024 · published calculated checkpoint</summary><pre>"
        + html.escape(json.dumps(records[-1]["original"], indent=2))
        + "</pre></details><p>"
        + citations
        + "</p>"
        + """
<script>const figures="""
        + payload
        + ";const captions="
        + json.dumps([c for _, c, _ in graphs])
        + ";function show(i){document.getElementById('caption').textContent=captions[i];document.querySelectorAll('nav button').forEach((b,j)=>b.classList.toggle('active',j===i));Plotly.react('plot',figures[i].data,figures[i].layout,{responsive:true,displaylogo:false,toImageButtonOptions:{format:'svg',filename:'argon-data'}})}show(0);</script></main></body></html>"
    )
    (out / "argon-data-atlas.html").write_text(page)
    fields = [k for k in records[0] if k != "original"]
    with (out / "argon-normalized-pv.csv").open("w") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows({k: r[k] for k in fields} for r in records)
    report = dict(
        git_revision=subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT, text=True
        ).strip(),
        source_inventory=inventory,
        normalization="angstrom^3 per atom; fcc cell /4, hcp cell /2; Vm*1e24/NA",
        normalized_rows=len(records),
        primary_compression_rows=len(measured),
        near_sublimation_rows=len(neutron),
        room_temperature_rows=len(room),
        cryogenic_rows=len(cold),
        counts_by_source_and_role=dict(
            Counter(r["source"] + ": " + r["role"] for r in records)
        ),
        exclusions={
            "Dewaele pressure missing": 8,
            "Anderson unreliable zero sentinels": 15,
        },
        note="Rows from earlier precursor/replotted Ross and adopted Grimsditch densities are non-independent, hidden initially. No source fitting or recalibration.",
    )
    (out / "argon-data-inventory.json").write_text(json.dumps(report, indent=2) + "\n")
    (out / "README.md").write_text(
        "# Argon data atlas\n\nOffline interactive atlas: `argon-data-atlas.html`. Static figures are supplied as PNG, PDF and SVG.\n\n"
        + f"{len(measured)} source compression rows and {len(neutron)} Xiao near-sublimation volume rows appear in the experimental overview. These counts include repeated measurements; they are not independent-sample counts.\n\n"
        + f"All {len(tables)} CSV source tables are embedded in the HTML. The inventory lists source hashes and counts. Normalized P-V export also includes calculated, adopted and overlapping graphical data, with explicit role labels. Missing pressures remain blank. Volume is per atom. Original CSV tables are unchanged.\n\n"
        + "Maltby has no newly recovered independent experimental table; its published Table 8 calculated checkpoint is available as an optional layer. Grimsditch acoustic shifts and elastic bounds have separate views. Phase observations and parameter summaries remain in the table archive.\n"
    )
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
