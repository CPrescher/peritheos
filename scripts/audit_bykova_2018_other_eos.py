"""Read-only numerical audit of Bykova (2018) EOSs outside the IV/V BM2 work.

Repository maintenance tool. Default prints JSON; --write writes only this
audit's results/plot. No material, dataset or shared ledger is modified.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.optimize import brentq, least_squares

from peritheos.eos.rt.bm import BM3

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
READINGS = ROOT / "docs/data/bykova-2018-other-eos-readings.json"
RESULT = ROOT / "docs/data/bykova-2018-other-eos-audit.json"
PLOT = ROOT / "docs/literature-reproductions/assets/bykova-2018-other-eos-checks.svg"
EXPERIMENTAL = (547.2, 103.0, 3.02)
V_DIAGNOSTIC = (427.3921208, 185.26, 3.10)
SOURCES = {
    "output/bykova-2018-other-eos-audit-2026-10-24/AnaCernok_thesis.pdf": "093882fd5953cc0eab8517e6ed4ddc4e6003628097956b0b5029d5a4822fdd2c",
    "output/coesite-iv-bykova-2018-audit-2026-10-23/article.pdf": "776ca4ca61872b7540c1ebaa70c29046f84a41c2d0700725d009eba203bd3cb5",
    "output/coesite-iv-bykova-2018-audit-2026-10-23/supplement.pdf": "9d4bcb36ae6a89b19f5bbcc3fc671f3c1bfc6530a76915382c6132d3d907e459",
    "peritheos/data/datasets/coesite-i-ii-cernok-2014-table1-pv.csv": "733b94709bd150635b1fe5e7975a652055cc745e604dad723f2b06146b3ce575",
    "peritheos/data/datasets/coesite-ii-iii-bykova-2018-table2-pv.csv": "8aca860f04909fc827fc09ea6e24039d0364dcc21f3c774b07cb03cdcfe945e2",
    "peritheos/data/datasets/coesite-v-bykova-2018-table10-calc-pv.csv": "fd4ffe396fe224622cadaec3b4688be816af407c3996f950bfb439fcb6d68b31",
}


def pressure(volume: Any, parameters: Any) -> Any:
    """Standard BM3, independently written in Eulerian finite strain."""
    v0, k0, kp = parameters
    strain = 0.5 * ((v0 / np.asarray(volume)) ** (2.0 / 3.0) - 1.0)
    return 3 * k0 * strain * (1 + 2 * strain) ** 2.5 * (1 + 1.5 * (kp - 4) * strain)


def inverse_pressure(p: float, parameters: Any) -> float:
    return brentq(
        lambda v: pressure(v, parameters) - p, 0.6 * parameters[0], parameters[0]
    )


def summary(residual: Any) -> dict[str, float]:
    r = np.asarray(residual)
    return {
        "rms_gpa": float(np.sqrt(np.mean(r**2))),
        "maximum_absolute_gpa": float(np.max(np.abs(r))),
    }


def graphical(
    section: dict[str, Any], key: str, parameters: Any, bound_factor: float = 1
) -> dict[str, Any]:
    pt = np.array(section["pressure_ticks_gpa_and_x_pixel"])
    vt = np.array(section["volume_ticks_a3_per_formula_unit_and_y_pixel"])
    pm = np.polyfit(pt[:, 1], pt[:, 0], 1)
    vm = np.polyfit(vt[:, 1], vt[:, 0], 1)
    bound = section["coordinate_bound_pixels"] * bound_factor
    dp, dv = abs(pm[0]) * bound, abs(vm[0]) * bound
    rows = []
    for i, (x, y) in enumerate(section[key], 1):
        p, v = float(np.polyval(pm, x)), float(np.polyval(vm, y))
        r = float(pressure(v * 16, parameters) - p)
        lower = float(pressure((v + dv) * 16, parameters) - (p + dp))
        upper = float(pressure((v - dv) * 16, parameters) - (p - dp))
        rows.append(
            {
                "reading": i,
                "x_pixel": x,
                "y_pixel": y,
                "pressure_gpa": p,
                "volume_a3_per_formula_unit": v,
                "residual_gpa": r,
                "residual_bound_gpa": [lower, upper],
                "curve_crosses_reading_box": lower <= 0 <= upper,
            }
        )
    return {
        "rows": rows,
        "pressure_bound_gpa": float(dp),
        "volume_bound_a3_per_formula_unit": float(dv),
        **summary([r["residual_gpa"] for r in rows]),
        "incompatible_readings": [
            r["reading"] for r in rows if not r["curve_crosses_reading_box"]
        ],
    }


def reproduce() -> dict[str, Any]:
    # Source PDFs are a local audit archive; published CSV/readings suffice to
    # replay elsewhere. When present, verify that the archive is unchanged.
    for name, digest in SOURCES.items():
        path = ROOT / name
        if path.exists():
            assert hashlib.sha256(path.read_bytes()).hexdigest() == digest, name
        elif not name.startswith("output/"):
            raise FileNotFoundError(path)
    readings = json.loads(READINGS.read_text())
    material = json.loads(
        (ROOT / "peritheos/data/materials/coesite.eosmat").read_text()
    )
    record = next(
        r
        for r in material["eos_records"]
        if r["identifier"] == "coesite_i_iii_bykova_2018_300k_bm3"
    )
    assert [record["parameters"][k]["value"] for k in ("V0", "K0", "K0_prime")] == list(
        EXPERIMENTAL
    )
    held = json.loads((ROOT / "peritheos/data/materials/coesite_v.eosmat").read_text())
    assert held["eos_records"] == []
    assert (
        held["source"]["held_eos_candidate"]["diagnostic_reconstruction"]["V0_a3"]
        == V_DIAGNOSTIC[0]
    )
    rows = []
    for filename in (
        "coesite-i-ii-cernok-2014-table1-pv.csv",
        "coesite-ii-iii-bykova-2018-table2-pv.csv",
    ):
        with (DATA / filename).open(newline="") as stream:
            rows.extend(
                dict(row, source_file=filename) for row in csv.DictReader(stream)
            )
    assert len(rows) == 24
    p = np.array([float(r["pressure_gpa"]) for r in rows])
    v = np.array([float(r["volume_a3_z16_equivalent_cell"]) for r in rows])
    for row in rows:
        assert np.isclose(
            float(row["volume_a3_z16_equivalent_cell"]),
            float(row["volume_a3_source_cell"])
            * 16
            / float(row["formula_units_per_source_cell"]),
            atol=1e-8,
        )
    original = pressure(v, EXPERIMENTAL) - p
    # Replay the one existing diagnostic exactly; no alternative masks,
    # weighting or replacement model is explored.
    fit = least_squares(
        lambda x: pressure(v, x) - p,
        EXPERIMENTAL,
        bounds=((273.6, 10.3, -10), (820.8, 1030, 20)),
        max_nfev=5000,
    )
    assert fit.success
    fitted = pressure(v, fit.x) - p
    model = BM3(V0=EXPERIMENTAL[0], K0=EXPERIMENTAL[1], K0_prime=EXPERIMENTAL[2])
    np.testing.assert_allclose(model.pressure(v), original + p, atol=1e-10, rtol=1e-12)
    for i, row in enumerate(rows):
        row["published_curve_residual_gpa"] = float(original[i])
        row["existing_diagnostic_residual_gpa"] = float(fitted[i])
    phase_results = {}
    for phase in ("coesite-I", "coesite-II", "coesite-III"):
        mask = np.array([r["phase"] == phase for r in rows])
        phase_results[phase] = {"rows": int(sum(mask)), **summary(original[mask])}
    fig1 = graphical(readings["figure1"], "curve_pixels", EXPERIMENTAL)
    fig6 = readings["figure6"]
    markers = graphical(fig6, "marker_pixels", V_DIAGNOSTIC)
    curve = graphical(fig6, "curve_pixels", V_DIAGNOSTIC)
    extra = [
        r["residual_gpa"]
        for r in markers["rows"]
        if r["reading"] != fig6["anchor_marker"]
    ]
    rounding_v0 = [
        brentq(lambda v0: pressure(342.5716, (v0, 185.26, 3.1)) - a, 400, 460)
        for a in (56.5, 57.5)
    ]
    transformed_pressure = float(pressure(352.3491, V_DIAGNOSTIC))
    # Reproduce the previously audited IV result without writing its files.
    from scripts.reproduce_bykova_2018_coesite_iv import reproduce as replay_iv

    iv = replay_iv()
    return {
        "source_doi": "10.1038/s41467-018-07265-z",
        "audit_date": "2026-10-24",
        "source_sha256": SOURCES,
        "residual_definition": "P(BM3 at printed/digitized volume) minus source P, GPa",
        "experimental_i_ii_iii": {
            "published_parameters_z16": list(EXPERIMENTAL),
            "source_rows": rows,
            "phase_results": phase_results,
            "published_curve": summary(original),
            "existing_partial_unweighted_diagnostic": {
                "parameters_z16": fit.x.tolist(),
                **summary(fitted),
                "rss_per_degree_of_freedom_gpa2": float(sum(fitted**2) / 21),
                "qualification": "Descriptive residual scale, not statistical reduced chi-square; author objective/selection/weights unavailable.",
            },
            "figure1_fitted_curve": fig1,
            "python_bm3_maximum_difference_gpa": float(
                np.max(abs(model.pressure(v) - pressure(v, EXPERIMENTAL)))
            ),
        },
        "theoretical_iv_existing_audit": iv,
        "theoretical_v": {
            "diagnostic_parameters_z16": list(V_DIAGNOSTIC),
            "anchor_residual_gpa": float(pressure(342.5716, V_DIAGNOSTIC) - 57),
            "nearest_integer_rounding_assumption_v0_interval_a3": rounding_v0,
            "figure6_calculated_markers": markers,
            "additional_markers_excluding_anchor": {
                "rows": len(extra),
                **summary(extra),
            },
            "figure6_smooth_curve": curve,
            "eight_pixel_marker_sensitivity": graphical(
                fig6, "marker_pixels", V_DIAGNOSTIC, 2
            ),
            "transformed_table10_state": {
                "source_label": "coesite-IV starting structure; text says relaxation transforms to V",
                "printed_pressure_gpa": 55,
                "printed_volume_a3": 352.3491,
                "diagnostic_v_pressure_gpa": transformed_pressure,
                "residual_gpa": transformed_pressure - 55,
                "qualification": "Not a second unambiguously assigned equilibrium coesite-V EOS datum. Pre/post-relaxation pressure interpretation unresolved.",
            },
        },
        "other_source_owned_am05": [
            {
                "source_label": "stishovite/CaCl2",
                "K0_gpa": 298.90,
                "K0_prime": 3.66,
                "V0": None,
                "parameter_uncertainty": None,
                "finding": "New Table 11 result; zero-pressure volume and numerical grid absent; combined branch/transition mapping unresolved.",
            },
            {
                "source_label": "seifertite",
                "K0_gpa": 310.14,
                "K0_prime": 4.09,
                "V0": None,
                "parameter_uncertainty": None,
                "finding": "New Table 11 result; zero-pressure volume and numerical grid absent.",
            },
        ],
    }


def plot(result: dict[str, Any]) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(
        2, 2, figsize=(12, 8), gridspec_kw={"height_ratios": [2, 1]}
    )
    exp = result["experimental_i_ii_iii"]
    ax, resax = axes[:, 0]
    rows = exp["source_rows"]
    for phase, marker, color in zip(
        ("coesite-I", "coesite-II", "coesite-III"),
        ("o", "s", "^"),
        ("#38855c", "#75539b", "#cf722a"),
    ):
        subset = [r for r in rows if r["phase"] == phase]
        p = [float(r["pressure_gpa"]) for r in subset]
        v = [float(r["volume_a3_z16_equivalent_cell"]) / 16 for r in subset]
        ax.errorbar(
            p,
            v,
            xerr=[float(r["pressure_sigma_gpa"]) for r in subset],
            yerr=[float(r["volume_sigma_a3_z16_equivalent_cell"]) / 16 for r in subset],
            fmt=marker,
            ms=4,
            color=color,
            label=phase,
        )
        resax.scatter(
            p,
            [r["published_curve_residual_gpa"] for r in subset],
            marker=marker,
            s=17,
            color=color,
        )
    ps = np.linspace(0, 37, 160)
    ax.plot(
        ps,
        [inverse_pressure(p, EXPERIMENTAL) / 16 for p in ps],
        color="black",
        label="Published BM3",
    )
    diagnostic = exp["existing_partial_unweighted_diagnostic"]["parameters_z16"]
    ax.plot(
        ps,
        [inverse_pressure(p, diagnostic) / 16 for p in ps],
        "--",
        color="#777777",
        label="Existing subset diagnostic",
    )
    resax.scatter(
        [float(r["pressure_gpa"]) for r in rows],
        [r["existing_diagnostic_residual_gpa"] for r in rows],
        marker="x",
        color="#777777",
        s=14,
        label="Subset diagnostic",
    )
    ax.set_title("Experimental I/II/III: partial 24-row replay")
    ax.legend(fontsize=8)
    ax, resax = axes[:, 1]
    theory = result["theoretical_v"]
    ps = np.linspace(25, 66, 160)
    ax.plot(
        ps,
        [inverse_pressure(p, V_DIAGNOSTIC) / 16 for p in ps],
        color="#633b7b",
        label="Held anchor-based diagnostic",
    )
    for key, label, marker, color in (
        (
            "figure6_calculated_markers",
            "Calculated markers (digitized)",
            "^",
            "#8c609c",
        ),
        ("figure6_smooth_curve", "Smooth curve (digitized)", ".", "#333333"),
    ):
        data = theory[key]
        x = [r["pressure_gpa"] for r in data["rows"]]
        y = [r["volume_a3_per_formula_unit"] for r in data["rows"]]
        ax.errorbar(
            x,
            y,
            xerr=data["pressure_bound_gpa"],
            yerr=data["volume_bound_a3_per_formula_unit"],
            fmt=marker,
            color=color,
            ms=4,
            label=label,
        )
        residual = np.array([r["residual_gpa"] for r in data["rows"]])
        limits = np.array([r["residual_bound_gpa"] for r in data["rows"]])
        resax.errorbar(
            x,
            residual,
            yerr=np.array([residual - limits[:, 0], limits[:, 1] - residual]),
            fmt=marker,
            color=color,
            ms=4,
        )
    ax.scatter(
        [57],
        [342.5716 / 16],
        marker="*",
        s=85,
        color="#b56e17",
        label="Integer-pressure anchor",
    )
    ax.set_title("Theoretical V: graphical checks, no new fit")
    ax.legend(fontsize=8)
    for ax in axes[0]:
        ax.set_ylabel("Volume per SiO2 (Å³)")
        ax.grid(alpha=0.15)
    for ax in axes[1]:
        ax.axhline(0, color="#777777", lw=0.8)
        ax.set_ylabel("Pressure residual (GPa)")
        ax.set_xlabel("Pressure (GPa)")
        ax.grid(alpha=0.15)
    fig.text(
        0.5,
        0.01,
        "Bykova et al. (2018), doi:10.1038/s41467-018-07265-z. Graphical bounds are analyst reading widths; no confidence band.",
        ha="center",
        fontsize=8,
    )
    fig.tight_layout(rect=(0, 0.035, 1, 1))
    PLOT.parent.mkdir(parents=True, exist_ok=True)
    preview = ROOT / "output/bykova-2018-other-eos-audit-2026-10-24/checks.png"
    preview.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(PLOT)
    fig.savefig(preview, dpi=150)
    plt.close(fig)


def compare(actual: Any, expected: Any) -> None:
    if isinstance(actual, dict):
        assert actual.keys() == expected.keys()
        for key in actual:
            compare(actual[key], expected[key])
    elif isinstance(actual, list):
        assert len(actual) == len(expected)
        for a, b in zip(actual, expected):
            compare(a, b)
    elif isinstance(actual, (float, int)) and not isinstance(actual, bool):
        np.testing.assert_allclose(actual, expected, rtol=1e-6, atol=1e-7)
    else:
        assert actual == expected, (actual, expected)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--check", action="store_true")
    parser.add_argument("--plot", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    if args.write:
        RESULT.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    if args.check:
        compare(result, json.loads(RESULT.read_text()))
        print("Bykova other-EOS audit replay matches saved numerical evidence.")
    if args.plot:
        plot(result)
    if not any((args.write, args.check, args.plot)):
        print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
