"""Render the conditional Sakai network with source/refit curves and residuals."""

from __future__ import annotations

import argparse
import json

import numpy as np

from scripts.reconstruction_diagnostics import DiagnosticComparison, panel
from scripts.reconstruction_figures import current_curve, curve_svg
from scripts.report_numbers import number, value_error
from scripts.report_plot_themes import dark_path, dark_svg, themed_images
from scripts.reproduce_sakai_2025_gold import pressure
from scripts.reproduce_sakai_2025_network import OUTPUT, ROOT

PAGE = ROOT / "docs/literature-reproductions/sakai-2025-network.md"


def comparison(stage, mode):
    rows = stage["reductions"]
    v = np.array([r["target_volume_a3"] for r in rows])
    p = np.array([r["calculated_pressure_gpa"] for r in rows])
    if stage["anchor"]:
        v = np.append(v, stage["anchor"]["volume_a3"])
        p = np.append(p, stage["anchor"]["pressure_gpa"])
    grid = np.linspace(v.min(), v.max(), 160)
    source, refit = [
        pressure(v, stage[k]) for k in ("published_parameters", "fitted_parameters")
    ]
    return DiagnosticComparison(
        [
            panel(
                stage["target"] + ": " + mode.replace("_", " "),
                v,
                p,
                source,
                refit,
                "Volume (Å³/cell)",
                "Pressure (GPa)",
                "GPa",
                swap_axes=True,
                observed_label="Recovered pressure reductions / reference constraint"
                if mode == "sequential_with_recovered_history"
                else "Paired-volume pressure reductions / named constraint",
                curves=[
                    (
                        "300 K",
                        grid,
                        pressure(grid, stage["published_parameters"]),
                        pressure(grid, stage["fitted_parameters"]),
                    )
                ],
            )
        ],
        "Conditional subset replay, with reduced pressures dependent on upstream calibrants; not the complete original fit.",
        source_title="Published coefficients on subset",
    )


def render(report, check=False):
    lines = [
        "# Sakai (2025): measured-pair sequential calibration",
        "",
        "Audit date: 2026-10-19. **Eight conditional subset calibration stages are executable.** "
        "Published catalog coefficients and classifications are preserved. A coefficient mismatch "
        "on an incomplete subset does not establish failure of the full published fit.",
        "",
        "[Source article](https://doi.org/10.1038/s43246-025-00792-5) · "
        "[Nine-material audit](sakai-2025-nine-materials.md) · [Full five-family Au check](sakai-2025-gold.md) · "
        "[Numerical results, hashes and every residual](../data/sakai-2025-network-reconstruction.json)",
        "",
        "## Inputs and method",
        "",
        "The checksum-pinned Supplementary Data 1 cell archive supplies **393 measured volume pairs** "
        "on S1–S7, including all 48 Mo–W and 19 Fe–W pairs. S6 stores W volume in column J; "
        "D/F are Fe lattice a/c. Original numeric tokens, source cells, run markers and missing "
        "errors are retained in `peritheos/data/datasets/sakai-2025-network-paired-volumes.csv`. "
        "These rows overlap the earlier Au ancestor archive and represent the same experiments.",
        "",
        "The NaCl control adds 103 existing Sakai2011/Dorfman2012 measured pairs and the "
        "32.3 GPa / 27.25 Å³ Matsui reference point identified in Table S8. This is one equal-residual "
        "reference constraint, not a new paired observation or hard equality. Sata2002 peak-combination "
        "choices and the other Table S8 earlier-study constraints are omitted from these subset controls. "
        "These omissions describe this audit's scope, not exclusions attributed to the original authors.",
        "",
        "Cu stays at the published Table 2 anchor. Published-upstream controls hold the other "
        "Table 2 calibrants fixed. Sequential fits carry fitted Re into Pt, Pt into MgO, W into "
        "Au/Mo/Fe, and fitted Pt/MgO/Au into NaCl. V0/K0/K∞′ remain fixed except that Fe and "
        "NaCl refine V0/K0/K0′, as in Table 2. **No Table S9 calculated output is fitted.** "
        "Equal pressure residuals are an explicit analyst assumption. Numerical bounds keep positive "
        "K0 and reference volume above observed volumes; the implemented finite K′/V0/K0 bounds are "
        "diagnostic safeguards, not source constraints.",
        "",
        "**Errors:** parameter errors below are conditional residual-scaled regression standard errors. "
        "They omit upstream uncertainty, shared calibrant covariance and fixed-input errors. Source "
        "`err.` fields do not establish this statistical convention. RMS is an unweighted descriptive "
        "pressure error. Residual = Pmodel(Vtarget) − Pcalibrant(Vcalibrant). Both curves use identical "
        "inputs. Pressure is horizontal and volume vertical. The NaCl reference constraint enters its RMS.",
        "",
        "## Results",
        "",
        "Re from S1 gives K0′=4.410969, rounding to the published 4.411 under the stated Cu anchor "
        "and objective. This is conditional central-value agreement. Other subsets show the effect "
        "of restricted observations and sequential calibration; the earlier full five-family Au "
        "central-value parity result remains separate. Original weights, covariance, complete historical "
        "constraints and full-network workflow parity remain unresolved.",
        "",
        "### Recovered historical constraints",
        "",
        "The third mode adds traceable earlier inputs: Dewaele2004 Cu–W rows paired by "
        "their shared Table I coordinates, Anzellini2014 ruby run 1 and measured W run 3, "
        "Dewaele2004 Pt/ruby, Dewaele2000 300 K MgO/Pt rows with a reported Pt lattice "
        "parameter, Dorfman2012 MgO/Pt and Mo/MgO pairs, Dewaele2008 Mo/ruby, the "
        "previously audited five-family 210-state Au selection, Dewaele2006 hcp-Fe/W, "
        "and Sakai2014 300 K Fe/NaCl-B2. Blank calibrant coordinates and Fe/NaCl-B1 "
        "states are excluded; they are not replaced by back-calculated marker volumes. "
        "Ruby pressures are rounding-limited coordinate reconstructions using the "
        "source Mao1986 or Dorogokupets2007 equations and Shen2020, not recovered optical spectra.",
        "",
        "This historical control fits W before Re because Re's older input uses W. "
        "It then carries Re → Pt → MgO, W/Re/MgO → Au, W/MgO → Mo, Pt/MgO/Au → NaCl "
        "and W/NaCl → Fe. Its W, Re and Au coefficients move close to their source values. "
        "The upstream-fitted Au K0′ is approximately 5.78075; the earlier 5.78010 central-value "
        "parity result conditions on published calibrants and remains distinct. Some original "
        "constraints, row masks and weights remain unavailable or ambiguous, so this is still "
        "a conditional reconstruction rather than full author-workflow reproduction. "
        "The source hashes and every contributing source/row identity are in the numerical report.",
        "",
    ]
    for mode, stages in report["modes"].items():
        lines += [
            "### " + mode.replace("_", " ").capitalize(),
            "",
            "| Material | Inputs | Fitted parameters ± conditional errors | Published RMS (GPa) | Refit RMS (GPa) |",
            "|---|---:|---|---:|---:|",
        ]
        for s in stages:
            indices = [0, 1, 2] if len(s["free_parameters"]) == 3 else [2]
            coefs = []
            for name, i, error in zip(
                s["free_parameters"], indices, s["conditional_parameter_errors"]
            ):
                value, width = value_error(s["fitted_parameters"][i], error)
                unit = " Å³/cell" if name == "V0" else " GPa" if name == "K0" else ""
                coefs.append(f"{name} = {value} ± {width}{unit}")
            lines.append(
                f"| {s['target']} | {s['observations']} | {'; '.join(coefs)} | {number(s['published_rms_gpa'], 3)} | {number(s['refit_rms_gpa'], 3)} |"
            )
        lines += [""]
        for s in stages:
            c = comparison(s, mode)
            rel = f"../data/sakai-2025-network/{mode}-{s['target']}.svg"
            path = PAGE.parent / rel
            if check:
                if not path.exists() or not current_curve(path.read_bytes(), c):
                    raise ValueError(f"Stale Sakai figure: {path.name}")
                if (PAGE.parent / dark_path(rel)).read_bytes() != dark_svg(
                    path.read_bytes()
                ):
                    raise ValueError(f"Stale dark figure: {path.name}")
            else:
                path.parent.mkdir(parents=True, exist_ok=True)
                if not path.exists() or not current_curve(path.read_bytes(), c):
                    path.write_bytes(curve_svg(c))
                (PAGE.parent / dark_path(rel)).write_bytes(dark_svg(path.read_bytes()))
            lines += [
                '<details markdown="1">',
                f"<summary>{s['target']}: curves and residuals</summary>",
                "",
                *themed_images(
                    s["target"] + " conditional subset curves and residuals", rel
                ),
                "",
                "</details>",
                "",
            ]
    lines += [
        "## Reproduction",
        "",
        "Repository-only maintenance commands:",
        "",
        "```sh",
        "uv run python -m scripts.reproduce_sakai_2025_network",
        "uv run python -m scripts.report_sakai_2025_network",
        "# Both commands also support --check",
        "```",
        "",
        "Source rights and the contributor-only CC0 scope remain those of "
        "[the original recovery](sakai-2025-gold.md#rights-and-integration-scope). "
        "No author outreach or replacement EOS is part of this audit.",
        "",
        "## Public dataset example",
        "",
        "The paired observations are accessible through the installed package API:",
        "",
        "```python",
        "from peritheos import Material, get_material_document",
        "",
        "mgo = Material.from_eosmat(",
        '    get_material_document("mgo"), record_identifiers=["mgo_sakai_2025_rydberg_stacey_1"]',
        ")",
        'pairs = mgo.get_dataset("mgo_sakai_2025_table_s7_measured_pairs")',
        "print(len(pairs))  # 81 measured MgO–Pt volume pairs",
        'print(pairs["target_volume_a3"], pairs["calibrant_volume_a3"])',
        "```",
        "",
    ]
    content = "\n".join(lines)
    if check:
        if PAGE.read_text(encoding="utf-8") != content:
            raise ValueError("Stale Sakai network page")
    else:
        PAGE.write_text(content, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    # Plot saved results; rendering must never rerun a regression.
    render(
        json.loads(OUTPUT.read_text(encoding="utf-8")), check=parser.parse_args().check
    )
