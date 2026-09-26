"""Test a combined-data interpretation without changing the published EOS.

Run: python -m scripts.fit_argon_errandonea_combined
Ross's exact Table I replaces indistinguishable open circles in Fig. 5.
The cryogenic diamond enters only explicitly mixed-temperature diagnostics.
"""

import csv
import json

import numpy as np
from scipy.optimize import brentq, least_squares

from scripts.reproduce_argon_errandonea_2006 import ROOT, bm3

INPUT = ROOT / "docs/data/argon-errandonea-2006-combined-inputs.json"
REPORT = ROOT / "docs/data/argon-errandonea-2006-combined-fit.json"
FIGURE = ROOT / "docs/data/argon-errandonea-2006-combined-fit.png"
PUBLISHED = np.array([143.0, 6.5, 5.1])
ERRORS = np.array([11.0, 0.5, 0.3])
PREFERRED = "room_temperature_without_ross_row17"


def observations():
    source = json.loads(INPUT.read_text())
    own_path = (
        ROOT / "peritheos/data/datasets/argon-fcc-errandonea-2006-figure5-digitized.csv"
    )
    with own_path.open() as stream:
        own = list(csv.DictReader(stream))
    result = [
        dict(
            source="errandonea_2006",
            source_row=i + 1,
            pressure_gpa=float(r["pressure_gpa"]),
            volume_a3=float(r["volume_a3"]),
            temperature_k=300,
        )
        for i, r in enumerate(own)
    ]
    result += [
        dict(
            source="ross_1986",
            source_row=int(r["source_row"]),
            pressure_gpa=float(r["pressure_gpa"]),
            volume_a3=float(r["volume_a3"]),
            temperature_k=298,
        )
        for r in source["ross"]["rows"]
    ]
    marker = source["anderson_swenson_marker"]
    result.append(
        dict(
            source="anderson_swenson_1975",
            source_row=1,
            pressure_gpa=marker["pressure_gpa"],
            volume_a3=marker["volume_a3"],
            temperature_k=None,
        )
    )
    for row in result:
        suspect = row["source"] == "ross_1986" and row["source_row"] == 17
        row["included_in_preferred_fit"] = (
            not suspect and row["temperature_k"] is not None
        )
        row["exclusion_reason"] = (
            "Suspected printed pressure error (247 kbar); adjacent near-identical volume is at 347 kbar. User requested exclusion; no replacement value asserted."
            if suspect
            else "Cryogenic marker, not room-temperature data."
            if row["temperature_k"] is None
            else None
        )
    return result


def inverse(pressure, pars):
    return np.array(
        [
            brentq(lambda v: bm3(v, *pars) - p, pars[0] * 0.05, pars[0], xtol=1e-11)
            if p > 0
            else pars[0]
            for p in pressure
        ]
    )


def fit(rows, objective):
    p = np.array([r["pressure_gpa"] for r in rows])
    v = np.array([r["volume_a3"] for r in rows])
    residual = (
        (lambda a: bm3(v, *a) - p)
        if objective == "pressure"
        else (lambda a: inverse(p, a) - v)
    )
    # Same bounds as the earlier own-data audit; expose any active bound.
    solutions = [
        least_squares(
            residual,
            start,
            bounds=([80, 0.01, 4], [300, 50, 15]),
            x_scale="jac",
            max_nfev=5000,
            ftol=1e-11,
            xtol=1e-11,
            gtol=1e-11,
        )
        for start in ([143, 6.5, 5.1], [180, 1.5, 7])
    ]
    best = min(solutions, key=lambda f: np.sum(f.fun**2))
    pars = best.x
    return dict(
        objective=f"equal_{objective}_weights",
        count=len(rows),
        parameters=dict(zip(("V0", "K0", "K0_prime"), pars.tolist())),
        pressure_rmse_gpa=float(np.sqrt(np.mean((bm3(v, *pars) - p) ** 2))),
        volume_rmse_a3=float(np.sqrt(np.mean((inverse(p, pars) - v) ** 2))),
        published_pressure_rmse_gpa=float(np.sqrt(np.mean((bm3(v) - p) ** 2))),
        solver_success=bool(best.success),
        active_bound_mask=best.active_mask.tolist(),
        parameters_inside_published_error_intervals=bool(
            np.all(abs(pars - PUBLISHED) <= ERRORS)
        ),
        note="Published error intervals have unspecified confidence; overlap is not statistical proof of fit reproduction.",
    )


def calculate():
    rows = observations()
    subsets = {
        "errandonea_only": [r for r in rows if r["source"] == "errandonea_2006"],
        "room_temperature_union": [r for r in rows if r["temperature_k"] is not None],
        "all_including_cryogenic_marker": rows,
        PREFERRED: [r for r in rows if r["included_in_preferred_fit"]],
    }
    return dict(
        status="diagnostic_combined_fit_not_original_fit_reproduction",
        preferred_fit=dict(
            subset=PREFERRED, objective="pressure", excluded_ross_rows=[17]
        ),
        source_counts={
            s: sum(r["source"] == s for r in rows)
            for s in ("errandonea_2006", "ross_1986", "anderson_swenson_1975")
        },
        published_parameters=dict(zip(("V0", "K0", "K0_prime"), PUBLISHED.tolist())),
        fit_bounds=dict(lower=[80, 0.01, 4], upper=[300, 50, 15]),
        qualifications=[
            "Eight recoverable Errandonea markers, not its full unavailable raw dataset.",
            "All 42 Ross Table I rows archived unchanged. Preferred fit excludes row17 (suspected printed pressure error); unfiltered fits retained as sensitivity checks.",
            "298 K Ross and 300 K Errandonea are pooled without a 2 K thermal correction.",
            "The Anderson–Swenson diamond is cryogenic; including it is a mixed-temperature sensitivity, not a valid 300 K isotherm.",
            "No experimental weighting or covariance reconstructed; equal-pressure and equal-volume objectives compared.",
            "No catalog coefficients or reproduction status changed.",
        ],
        fits={
            name: {
                objective: fit(data, objective) for objective in ("pressure", "volume")
            }
            for name, data in subsets.items()
        },
        observations=rows,
    )


def plot(report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
    rows = report["observations"]
    rt = list(report["fits"][PREFERRED]["pressure"]["parameters"].values())
    unfiltered = list(
        report["fits"]["room_temperature_union"]["pressure"]["parameters"].values()
    )
    for source, label, color, marker in [
        ("errandonea_2006", "Errandonea: 8 digitized, 300 K", "#d65a31", "o"),
        ("ross_1986", "Ross: 41 included, 298 K", "#2374ab", "o"),
        ("anderson_swenson_1975", "Anderson–Swenson: cryogenic", "#923c91", "D"),
    ]:
        selected = [
            r
            for r in rows
            if r["source"] == source
            and not (source == "ross_1986" and r["source_row"] == 17)
        ]
        p = np.array([r["pressure_gpa"] for r in selected])
        v = np.array([r["volume_a3"] for r in selected])
        axes[0].scatter(p, v, s=24, color=color, marker=marker, label=label, zorder=3)
        axes[1].scatter(p, bm3(v, *rt) - p, s=24, color=color, marker=marker)
    excluded = next(
        r for r in rows if r["source"] == "ross_1986" and r["source_row"] == 17
    )
    axes[0].scatter(
        excluded["pressure_gpa"],
        excluded["volume_a3"],
        s=50,
        color="0.45",
        marker="x",
        label="Ross row 17: excluded",
        zorder=4,
    )
    axes[1].scatter(
        excluded["pressure_gpa"],
        bm3(excluded["volume_a3"], *rt) - excluded["pressure_gpa"],
        s=50,
        color="0.45",
        marker="x",
    )
    p = np.linspace(0, 120, 301)
    for pars, label, style, color in [
        (PUBLISHED, "Published BM3", "--", "black"),
        (rt, "Preferred fit: row 17 excluded", "-", "#24824f"),
        (unfiltered, "Unfiltered RT fit (comparison)", ":", "#923c91"),
    ]:
        axes[0].plot(p, inverse(p, pars), style, color=color, lw=1.6, label=label)
    axes[0].set(xlabel="Pressure (GPa)", ylabel="Fcc cell volume (Å³)", xlim=(-2, 120))
    axes[0].legend(fontsize=7.8)
    axes[1].axhline(0, color="0.5", lw=1)
    axes[1].set(
        xlabel="Reported pressure (GPa)",
        ylabel="298/300 K fit − reported pressure (GPa)",
    )
    axes[1].annotate(
        "Ross row 17: excluded, source retained",
        xy=(24.7, float(bm3(66.28871956, *rt) - 24.7)),
        xytext=(36, 7),
        arrowprops=dict(arrowstyle="->"),
        fontsize=8,
    )
    fig.suptitle("Combined-data diagnostic • equal pressure weights", fontsize=13)
    fig.savefig(FIGURE, dpi=180)
    plt.close(fig)


if __name__ == "__main__":
    report = calculate()
    REPORT.write_text(json.dumps(report, indent=2) + "\n")
    plot(report)
    print(REPORT)
    for name, fits in report["fits"].items():
        for objective, result in fits.items():
            print(name, objective, result["parameters"], result["pressure_rmse_gpa"])
