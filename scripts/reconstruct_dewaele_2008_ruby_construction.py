"""Repository-only, partial source replay of Dewaele (2008) ruby construction.

This is an evidence audit, not a pressure standard or a recovered author fit.
Default execution needs NumPy/SciPy/Matplotlib. --recover-pdf additionally needs
pdfplumber and the exact supplied PDF. All products stay in the scoped output
directory; existing audit reports and material records are read-only.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import itertools
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/dewaele-2008-ruby-construction"
DATA = ROOT / "peritheos/data/datasets"
PDF_SHA = "a06fd5cb9433b99f35b20582be5214b1e805014c150a2907fdf8933776cf6daf"
METALS = ("Co", "Ni", "Mo", "Ag", "Al", "Cu", "Pt", "Ta", "W", "Au")
NAMES = dict(
    Co="cobalt_hcp",
    Ni="nickel",
    Mo="molybdenum",
    Ag="silver",
    Al="aluminum",
    Cu="copper",
    Pt="platinum",
    Ta="tantalum",
    W="tungsten",
    Au="gold",
)
# Direct transcription from 2008 Table IV. Columns: V0 atomic A3, fixed K0
# GPa, DAC K0', RSW K0', DAC/RSW printed 95% widths on K0'. No covariance.
TABLE4 = {
    "Co184": (11.096, 184.0, 4.73, 4.93, 0.06, 0.04),
    "Co190": (11.096, 190.0, 4.37, 4.67, 0.04, 0.03),
    "Ni": (10.942, 183.0, 4.99, 5.13, 0.02, 0.01),
    "Mo": (15.558, 261.0, 4.06, 4.18, 0.05, 0.02),
    "Ag": (17.021, 101.0, 5.97, 6.19, 0.02, 0.01),
}
# Preserve the damaged indexed transcription; do not silently decode it into
# an authoritative six-metal coefficient table. Direct PDF access needs login.
INDEXED_2004_RSW_TOKENS = {
    "Au": "167.5, 8.2",
    "Pt": "277.4, 6.6",
    "Ta": "194.3, 7.6",
    "W": "296.4, 4.5",
    "Cu": "133.5, 3.3",
    "Al": "73.4, 5.0",
}


def ruby(log_ratio, a=1920.0, b=9.61):
    return a / b * np.expm1(b * np.asarray(log_ratio))


def log_ratio(pressure, a=1904.0, b=9.5):
    return np.log1p(b * np.asarray(pressure) / a) / b


def vinet(volume_ratio, k0, kp):
    x = np.cbrt(np.asarray(volume_ratio))
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def compression(pressure, k0, kp):
    return np.array(
        [
            brentq(lambda r: float(vinet(r, k0, kp)) - p, 0.25, 1.0)
            for p in np.atleast_1d(pressure)
        ]
    )


def recover_figure3(pdf_path):
    import pdfplumber

    if hashlib.sha256(pdf_path.read_bytes()).hexdigest() != PDF_SHA:
        raise ValueError("Figure mapping requires the inspected publisher PDF")
    with pdfplumber.open(pdf_path) as pdf:
        page = pdf.pages[5]
        curve = page.curves[0]
        assert curve["stroking_color"] == (0.99608, 0.66275, 0.0)
        xy = np.array(curve["pts"])
        # Major tick centers, verified against printed axes and Figure caption.
        p = (xy[:, 0] - 340.008) / (516.262 - 340.008) * 150
        dp = (253.658988 - xy[:, 1]) / (253.658988 - 76.119988) * 20
        result = {
            "source_pdf_sha256": PDF_SHA,
            "source": "Dewaele et al. 2008, Figure 3, journal page 104102-6",
            "kind": "vector centerline of published ruby fit; calculated graph checkpoints",
            "pdf_curve_index": 0,
            "stroke_width_pt": curve["linewidth"],
            "mapping": {
                "x0_pt": 340.008,
                "x150_pt": 516.262,
                "y0_pt": 253.658988,
                "y20_pt": 76.119988,
            },
            "reading_halfwidth_gpa": 0.03,
            "points_mao_pressure_and_pressure_difference_gpa": np.c_[p, dp].tolist(),
            "qualification": "Centerline tolerance, not statistical uncertainty; not original calibration targets.",
        }
    (OUT / "figure3-source.json").write_text(json.dumps(result, indent=2) + "\n")


def source_inventory():
    rows, inputs = [], []
    for symbol in ("Al", "Cu", "Pt", "Ta", "W", "Au"):
        path = DATA / f"{NAMES[symbol]}-dewaele-2004-table1-compression.csv"
        items = list(csv.DictReader(path.open()))
        inputs.append(
            {
                "path": str(path.relative_to(ROOT)),
                "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                "rows": len(items),
            }
        )
        for index, row in enumerate(items, 1):
            rows.append(
                {
                    "source": "2004 Table I",
                    "source_material_row": index,
                    "metal": symbol,
                    "run": "not encoded in existing CSV",
                    "mao_token": row["ruby_pressure_classical_gpa"],
                    "revised_token": row["ruby_pressure_revised_gpa"],
                    "volume_token_a3_per_atom": row["atomic_volume_a3"],
                    "helium_gauge": False,
                    "source_eos_fit_included": True,
                    "ruby2008_calibration_selection": "unreported"
                    if symbol != "Au"
                    else "historical comparison only; 2008 uses revised Au source",
                }
            )
    path = DATA / "dewaele-2008-table2-complete.csv"
    items = list(csv.DictReader(path.open()))
    inputs.append(
        {
            "path": str(path.relative_to(ROOT)),
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "rows": len(items),
        }
    )
    symbols = {v: k for k, v in NAMES.items()}
    for index, row in enumerate(items, 1):
        if row["material"] not in symbols:
            continue  # Zn is not a member of the final ten-metal calibration.
        rows.append(
            {
                "source": "2008 Table II",
                "source_material_row": index,
                "metal": symbols[row["material"]],
                "run": row["run"],
                "mao_token": row["pressure_mao_token"],
                "revised_token": row["pressure_revised_token"],
                "volume_token_a3_per_atom": row["volume_token"],
                "helium_gauge": row["pressure_from_helium_eos"] == "True",
                "source_eos_fit_included": row["fit_included"] == "True",
                "ruby2008_calibration_selection": "unreported",
            }
        )
    for row in rows:
        # These are derived coordinates from rounded pressures, never recovered
        # spectrometer measurements. He-derived pressures have no ruby reading.
        is_ruby = not row["helium_gauge"]
        row["derived_ln_lambda_ratio_from_mao"] = (
            float(log_ratio(float(row["mao_token"]), b=7.665)) if is_ruby else None
        )
        row["derived_ln_lambda_ratio_from_2004"] = (
            float(log_ratio(float(row["revised_token"]))) if is_ruby else None
        )
    with (OUT / "pressure-volume-input-inventory.csv").open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    return rows, inputs


def table4_replay(graph):
    result = {}
    for name, (v0, k0, kd, ks, ed, es) in TABLE4.items():
        metal = "Co" if name.startswith("Co") else name
        xy = np.array(graph[metal]["pressure_difference_points_gpa"])
        p, dp = xy.T
        r = compression(p, k0, kd)
        model = vinet(r, k0, ks) - p
        # Shared fixed K0 and V0 are kept identical in both arms. This box is
        # only an interval compatibility diagnostic; it is not a 95% band.
        corners = []
        for bd, bs in itertools.product((kd - ed, kd + ed), (ks - es, ks + es)):
            cr = compression(p, k0, bd)
            corners.append(vinet(cr, k0, bs) - p)
        low, high = np.min(corners, axis=0), np.max(corners, axis=0)
        result[name] = {
            "kind": "published Vinet coefficient replay; not original shock reduction",
            "source": "2008 Table IV, columns 4 and 7; Section IV",
            "V0_atomic_a3_fixed_shared": v0,
            "K0_gpa_fixed_shared": k0,
            "DAC_K0_prime": kd,
            "RSW_K0_prime": ks,
            "printed_95pct_widths_DAC_RSW": [ed, es],
            "graph_checkpoint_count": len(p),
            "graph_domain_gpa": [float(p.min()), float(p.max())],
            "nominal_rms_graph_difference_gpa": float(
                np.sqrt(np.mean((model - dp) ** 2))
            ),
            "nominal_max_graph_difference_gpa": float(max(abs(model - dp))),
            "all_graph_points_within_reading_bound_0_03_gpa": bool(
                np.all(abs(model - dp) <= 0.03)
            ),
            "all_graph_points_in_coefficient_box_plus_reading_bound": bool(
                np.all((dp >= low - 0.03) & (dp <= high + 0.03))
            ),
            "coefficient_box_is_statistical_confidence_band": False,
            "original_shock_observations_recovered": False,
        }
    return result


def figure3_check():
    path = OUT / "figure3-source.json"
    if not path.exists():
        return {"status": "not recovered; run --recover-pdf with supplied PDF"}
    source = json.loads(path.read_text())
    p, dp = np.array(source["points_mao_pressure_and_pressure_difference_gpa"]).T
    logr = log_ratio(p, b=7.665)
    residual = ruby(logr) - p - dp
    fit = least_squares(
        lambda q: ruby(logr, *q) - p - dp,
        [1920.0, 9.61],
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    return {
        "kind": "independent source-figure equation check; not author optimization",
        "checkpoints": len(p),
        "mao_pressure_range_gpa": [float(min(p)), float(max(p))],
        "rms_difference_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_difference_gpa": float(max(abs(residual))),
        "passes_0_03_gpa_centerline_reading_bound": bool(np.all(abs(residual) <= 0.03)),
        "diagnostic_graph_fit_A_gpa_B": fit.x.tolist(),
        "author_workflow_reproduced": False,
    }


def paired_pressure_exception(rows):
    row = next(
        r
        for r in rows
        if r["metal"] == "Ni"
        and r["mao_token"] == "49.0"
        and r["revised_token"] == "51.0"
    )
    # Use only the source's adopted Ni branch, not a newly fitted replacement.
    p_ni = float(vinet(float(row["volume_token_a3_per_atom"]) / 10.942, 183, 4.99))
    return {
        "source": "2008 Table II, run 6, shared Ni/Zn state",
        "printed_mao_gpa": 49.0,
        "printed_revised_gpa": 51.0,
        "formula_revised_from_49_gpa": float(
            ruby(log_ratio(49.0, b=7.665), 1904.0, 9.5)
        ),
        "formula_mao_from_51_gpa": float(ruby(log_ratio(51.0), 1904.0, 7.665)),
        "Ni_table4_pressure_at_source_volume_gpa": p_ni,
        "finding": "Source pair is internally inconsistent. Ni EOS is fitted to these data, so proximity to one pressure cannot establish a corrected observation. No correction applied.",
    }


def dependency_map():
    """Carry shared ancestors into each pressure difference, without a covariance guess."""
    nodes = {
        "ruby2004": {"parents": ["comparison2004"], "kind": "shared input calibration"},
        "comparison2004": {
            "parents": ["RSW2004", "DAC2004_measured_volumes"],
            "kind": "upstream six-metal calibration",
        },
        "RSW2004": {
            "parents": [],
            "kind": "multiple older shock reductions; row lineage incomplete",
        },
        "DAC2004_measured_volumes": {"parents": [], "kind": "measured diffraction"},
        "helium1998": {
            "parents": [],
            "kind": "adopted secondary EOS; upstream calibration unresolved here",
        },
        "Au_stress2007": {
            "parents": [],
            "kind": "corrected diffraction input, numerical correction unrecovered",
        },
        "KennedyKeeler1972": {
            "parents": [],
            "kind": "shared compilation/reduction source",
        },
        "HixsonFritz1992": {
            "parents": [],
            "kind": "shared experimental/reduction source",
        },
        "Wang2002": {"parents": [], "kind": "shared MFP reduction model/source"},
        "Holmes1989": {"parents": [], "kind": "Pt shock reduction"},
        "Al_preliminary2007": {
            "parents": ["SWDB", "Al_zero_pressure_thermochemistry"],
            "kind": "initial fit; final joint Table I not substituted",
        },
        "SWDB": {
            "parents": [],
            "kind": "shock database; version and selected rows unavailable",
        },
        "Al_zero_pressure_thermochemistry": {
            "parents": [],
            "kind": "heat capacity, expansion and elasticity",
        },
    }
    shock_source = {
        "Co": "KennedyKeeler1972",
        "Ni": "KennedyKeeler1972",
        "Ag": "KennedyKeeler1972",
        "Mo": "HixsonFritz1992",
        "W": "HixsonFritz1992",
        "Cu": "Wang2002",
        "Ta": "Wang2002",
        "Au": "Wang2002",
        "Pt": "Holmes1989",
        "Al": "Al_preliminary2007",
    }
    for metal in METALS:
        elastic = f"{metal}_elastic_constraint"
        nodes[elastic] = {
            "parents": [],
            "kind": "shared fixed DAC/RSW K0"
            if metal in ("Co", "Ni", "Mo", "Ag")
            else "DAC elastic/reference constraint; older source branch unresolved",
        }
        parents = ["ruby2004", elastic]
        if metal == "Ni":
            parents.append("helium1998")
        if metal == "Au":
            parents.append("Au_stress2007")
        nodes[f"{metal}_DAC"] = {"parents": parents, "kind": "fitted DAC Vinet curve"}
        nodes[f"{metal}_RSW"] = {
            "parents": [shock_source[metal]]
            + ([elastic] if metal in ("Co", "Ni", "Mo", "Ag") else []),
            "kind": "fitted reduced-shock representation; full inputs incomplete",
        }
        nodes[f"{metal}_difference"] = {
            "parents": [f"{metal}_DAC", f"{metal}_RSW"],
            "kind": "common-volume pressure difference",
        }
    nodes["ruby2008"] = {
        "parents": [f"{m}_difference" for m in METALS] + ["ruby2004"],
        "kind": "output calibration; averaging and residual norm unresolved",
    }

    def ancestors(key, active=()):
        if key in active:
            raise ValueError("Calibration history must not contain an undated cycle")
        result = set(nodes[key]["parents"])
        for parent in nodes[key]["parents"]:
            result.update(ancestors(parent, (*active, key)))
        return result

    inherited = {m: sorted(ancestors(f"{m}_difference")) for m in METALS}
    shared = {
        f"{a}/{b}": sorted(set(inherited[a]) & set(inherited[b]))
        for a, b in itertools.combinations(METALS, 2)
    }
    result = {
        "nodes": nodes,
        "difference_ancestors": inherited,
        "shared_ancestors_by_metal_pair": shared,
        "numeric_covariance": None,
        "qualification": "Common ancestors propagate dependencies, not a numerical correlation coefficient. Incomplete shock genealogy may add further dependencies. The verified Table IV shared elastic constraints occur once per metal, not independently in their two arms. Identical constraints are not asserted for unresolved older branches.",
    }
    (OUT / "dependency-map.json").write_text(json.dumps(result, indent=2) + "\n")
    return {
        "node_count": len(nodes),
        "metal_pairs_with_shared_ancestors": sum(bool(v) for v in shared.values()),
        "numeric_covariance_recovered": False,
    }


def write_plot(report, graph):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(1, 2, figsize=(11.2, 4.4), layout="constrained")
    p = np.linspace(0, 160, 321)
    prior = report["existing_diagnostics"]
    for label, key in [
        ("Ruby fit to common-domain mean", "common_range_mean"),
        ("Ruby fit to pooled individual domains", "individual_ranges_pooled"),
        ("Ruby fit to graph continuations", "full_range_graph_continuation"),
    ]:
        q = prior[key]
        grid = p[p <= q["pressure_range_gpa"][1]] if key == "common_range_mean" else p
        axes[0].plot(
            grid, ruby(log_ratio(grid), q["A_gpa"], q["B"]) - grid, label=label
        )
    axes[0].plot(
        p,
        ruby(log_ratio(p)) - p,
        color="black",
        lw=2.2,
        label="Published A=1920, B=9.61",
    )
    axes[0].set(
        xlabel="2004 ruby pressure (GPa)", ylabel="Correction to 2004 pressure (GPa)"
    )
    axes[0].legend(fontsize=7)
    path = OUT / "figure3-source.json"
    if path.exists():
        p3, d3 = np.array(
            json.loads(path.read_text())[
                "points_mao_pressure_and_pressure_difference_gpa"
            ]
        ).T
        axes[1].plot(p3, d3, color="orange", lw=4, label="Figure 3 vector centerline")
        axes[1].plot(
            p3,
            ruby(log_ratio(p3, b=7.665)) - p3,
            color="black",
            ls="--",
            label="Printed ruby formula",
        )
    axes[1].set(
        xlabel="Mao 1986 ruby pressure (GPa)", ylabel="Correction to Mao pressure (GPa)"
    )
    axes[1].legend(fontsize=8)
    for ax in axes:
        ax.grid(alpha=0.2)
    fig.savefig(OUT / "calibration-comparison.png", dpi=180)
    fig.savefig(OUT / "calibration-comparison.svg")
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recover-pdf", type=Path)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.recover_pdf:
        recover_figure3(args.recover_pdf)
    graph = json.loads(
        (ROOT / "docs/data/dewaele-2008-followup-source.json").read_text()
    )["figure2b"]["metals"]
    rows, inputs = source_inventory()
    coverage = {}
    for metal in METALS:
        selected = [
            float(r["revised_token"].strip("()"))
            for r in rows
            if r["metal"] == metal and r["source_eos_fit_included"]
        ]
        curve = np.array(graph[metal]["pressure_difference_points_gpa"])
        coverage[metal] = {
            "selected_DAC_pressure_range_gpa": [min(selected), max(selected)]
            if metal != "Au"
            else None,
            "historical_Au_range_not_substituted_gpa": [min(selected), max(selected)]
            if metal == "Au"
            else None,
            "visible_graph_pressure_range_gpa": [
                float(curve[:, 0].min()),
                float(curve[:, 0].max()),
            ],
            "graph_endpoint_may_be_clipped_by_plot_boundary": True,
        }
    prior = json.loads((ROOT / "docs/data/dewaele-2008-followup.json").read_text())[
        "ruby_crosscalibration"
    ]
    report = {
        "format": "peritheos.dewaele-2008-ruby-construction.v1",
        "status": "partial_source_reconstruction; original optimization unresolved",
        "author_workflow_reproduced": False,
        "recommended_pressure_standard": False,
        "metals": list(METALS),
        "Zn_excluded": True,
        "supplied_2008_pdf_sha256": PDF_SHA,
        "input_files": inputs,
        "inventory_note": "Printed PV rows plus derived wavelength ratios; no measured wavelengths. Historical Au is retained but cannot replace the 2007 corrected Au input. He rows have no inferred ruby ratio.",
        "inventory_row_count": len(rows),
        "pressure_support": coverage,
        "source_table4_comparisons": table4_replay(graph),
        "aluminum": {
            "source": "2008 Section IV p.6; Dorogokupets-Oganov 2007 Section III p.4",
            "adopted_RSW_K0_gpa_K0_prime": [73.46, 4.52],
            "2007_final_Table_I_K0_gpa_K0_prime": [72.67, 4.62],
            "stage": "Initial shock plus zero-pressure thermochemical/elastic fit, before static DAC data are added",
            "final_table_is_not_the_adopted_input": True,
            "raw_SWDB_row_selection_weights_and_preliminary_full_thermal_coefficients": None,
        },
        "2004_indexed_table": {
            "access_level": "Publisher-indexed pages; direct PDF inaccessible without login",
            "raw_RSW_tokens": INDEXED_2004_RSW_TOKENS,
            "tokens_accepted_as_authoritative_numeric_coefficients": False,
            "reason": "Digit ordering damaged by extraction across bold K0 and K0-prime typography; no guessed numerical table is used for a full ten-metal replay.",
        },
        "figure3_published_formula_check": figure3_check(),
        "published_correction_at_2004_160_gpa": float(ruby(log_ratio(160.0)) - 160.0),
        "paired_pressure_exception": paired_pressure_exception(rows),
        "dependency_map": dependency_map(),
        "existing_diagnostics": prior,
        "uncertainty": {
            "shared_dependencies": [
                "2004 ruby A and B",
                "co-compressed sample ruby states",
                "same ultrasonic K0 in DAC and RSW arms",
                "Ni run-6 helium EOS",
                "shared shock/reduction compilations",
            ],
            "numeric_covariance_recovered": False,
            "propagation": "For r_i=P_SW_i(V_i)-P_ruby, C_r=J_shared C_shared J_shared^T+C_other. For an average, Var(mean)=w^T C_r w. Do not substitute independent diagonal errors for the shared term.",
        },
        "missing_for_author_reproduction": [
            "Undamaged 2004 Table II and actual shock-isotherm numerical inputs/fit selections",
            "Kennedy-Keeler 1972 Co/Ni/Ag isotherm rows and exact source shock/reduction lineage",
            "Hixson-Fritz 1992 Mo/W isotherm tables and thermodynamic reduction inputs",
            "2007 Dewaele-Loubeyre stress-corrected Au data, correction and 2004-scale fit",
            "Dorogokupets-Oganov preliminary Al full model, SWDB version/rows and weights",
            "Original ten-metal averaging support, grid/point choice, weights, residual norm, constraints/order and covariance",
            "Raw ruby wavelength readings and simultaneous run-6 He volumes",
        ],
    }
    (OUT / "reconstruction.json").write_text(json.dumps(report, indent=2) + "\n")
    with (OUT / "source-table4-vinet-pairs.csv").open("w", newline="") as stream:
        writer = csv.writer(stream)
        writer.writerow(
            [
                "branch",
                "V0_atomic_a3_fixed_shared",
                "K0_gpa_fixed_shared",
                "DAC_K0_prime",
                "RSW_K0_prime",
                "DAC_95pct_width",
                "RSW_95pct_width",
            ]
        )
        writer.writerows((name, *values) for name, values in TABLE4.items())
    write_plot(report, graph)
    print(
        json.dumps(
            {
                "status": report["status"],
                "inventory_rows": len(rows),
                "figure3": report["figure3_published_formula_check"],
                "table4": report["source_table4_comparisons"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
