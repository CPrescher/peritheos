"""Verify saved actual GUI runs; this command does not automate or emulate fitting."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import numpy as np

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from scripts.audit_fes_author_inputs import read_data, summary
from scripts.reconstruction_diagnostics import DiagnosticComparison, panel
from scripts.reconstruction_figures import current_curve, curve_svg
from scripts.report_numbers import number, value_error
from scripts.report_plot_themes import dark_path, dark_svg, themed_images
from scripts.run_fes_eosfit_console import saved_covariance

ROOT = Path(__file__).resolve().parents[1]
FOLDER = ROOT / "docs/data/fes-gui-2026-10-19"
AUTHOR = ROOT / "docs/data/fes-author-eosfit/sources"
PAGE = ROOT / "docs/literature-reproductions/fes-gui-replay.md"


def model(path):
    text = path.read_text(encoding="utf-8")
    parameters = {
        int(k): float(v) for k, v in re.findall(r"^Param =\s*(\d+)\s+(\S+)", text, re.M)
    }
    tref = float(re.search(r"^Tref =\s*(\S+)", text, re.M)[1])
    cold = BM3(parameters[1] / 10, parameters[2], parameters[3])
    thermal = int(re.search(r"^Thermal =\s*(\d+)", text, re.M)[1])
    return (
        (
            MieGruneisenDebye(
                cold,
                tref,
                parameters[11],
                parameters[18],
                parameters[19],
                parameters[13],
            )
            if thermal
            else cold
        ),
        parameters,
        tref,
    )


def residual_groups(data, residual):
    masks = {
        "all": np.ones(len(data), dtype=bool),
        "thermal": data[:, 2] > 300,
        "cold": data[:, 2] == 300,
    }
    return {
        key: summary(residual[mask], np.flatnonzero(mask))
        for key, mask in masks.items()
        if np.any(mask)
    }


def reproduce():
    source = AUTHOR / "FeS6EOSfittot-SataOhfuji-300K.dat"
    selected = read_data(source)
    cold, hot = selected[selected[:, 2] == 300], selected[selected[:, 2] > 300]
    np.testing.assert_array_equal(read_data(FOLDER / "cold.dat"), cold)
    np.testing.assert_array_equal(read_data(FOLDER / "hot.dat"), hot)
    cases = {}
    for name, data, initial, observed_gui_max in (
        ("gamma-unweighted", selected, AUTHOR / "FeS6.eos", 5.18),
        ("gamma-weighted", selected, AUTHOR / "FeS6.eos", 5.28),
        ("cold-weighted", cold, FOLDER / "cold-initial.eos", 4.19),
        ("staged-weighted", hot, FOLDER / "thermal-initial.eos", 5.29),
    ):
        fitted, params, tref = model(FOLDER / (name + ".eos"))
        original, _, _ = model(initial)

        def evaluate(eos):
            return (
                eos.pressure(data[:, 4] / 10)
                if name == "cold-weighted"
                else eos.pressure(data[:, 4] / 10, data[:, 2])
            )

        source_residual = evaluate(original) - data[:, 0]
        residual = evaluate(fitted) - data[:, 0]
        covariance = saved_covariance(
            FOLDER / (name + ".eos"), "cold" if name == "cold-weighted" else "morard"
        )
        keys = [1, 2, 3] if name == "cold-weighted" else [18]
        esds = np.sqrt(np.diag(covariance))
        maximum = float(np.max(np.abs(residual)))
        if abs(maximum - observed_gui_max) > 0.02:
            raise ValueError(f"Saved/native and GUI-displayed residual differ: {name}")
        cases[name] = {
            "rows": len(data),
            "source_initial_eos": str(initial.relative_to(ROOT)),
            "fitted_eos": str((FOLDER / (name + ".eos")).relative_to(ROOT)),
            "refined_parameters": {
                str(k): {"value": params[k], "error": float(e)}
                for k, e in zip(keys, esds)
            },
            "conditional_covariance": covariance.tolist(),
            "reference_temperature_k": tref,
            "source": residual_groups(data, source_residual)
            if name != "cold-weighted"
            else {"cold": summary(source_residual, np.arange(len(data)))},
            "refit": residual_groups(data, residual)
            if name != "cold-weighted"
            else {"cold": summary(residual, np.arange(len(data)))},
            "gui_displayed_max_abs_gpa": observed_gui_max,
            "gui_saved_native_max_difference_from_display_gpa": abs(
                maximum - observed_gui_max
            ),
            "residuals": [
                {
                    "row": i + 1,
                    "pressure_gpa": float(row[0]),
                    "temperature_k": float(row[2]),
                    "volume_cm3_mol": float(row[4]),
                    "initial_residual_gpa": float(s),
                    "refit_residual_gpa": float(r),
                }
                for i, (row, s, r) in enumerate(zip(data, source_residual, residual))
            ],
        }
    paths = [
        source,
        AUTHOR / "FeS6.eos",
        *sorted(FOLDER.glob("*.eos")),
        *sorted(FOLDER.glob("*.dat")),
        *sorted(FOLDER.glob("*.png")),
    ]
    return {
        "audit_date": "2026-10-19",
        "gui_version": "20211014",
        "executable_sha256": "e7d19f1f79cf544e0828b83c25a9a9e981e97c2fce925c1d895e9943fea5968e",
        "executable_source": "https://www.rossangel.com/Download/EoSFit7installer_MacOS.dmg",
        "execution": "Unmodified official Mac GUI in an authenticated temporary Xvfb display; extracted X11 dependencies relocated locally, no system installation.",
        "selection": "Complete confirmed selected 167-row author input, excluding the 11 rejected quenched rows. Cold stage uses 21 literature rows; thermal stage uses unchanged 146 hot rows, with duplicate retained.",
        "protocol": "Gamma-only on 167 rows with saved author cold model at Tr=298 K, equal weights and supplied P/V/T errors; separate three-parameter weighted BM3 cold fit at 300 K, followed by full-MGD gamma-only weighted fit on hot rows at Tr=300 K. Theta=417 K, n=2, q=1; q-compromise off.",
        "staged_handoff": "GUI data loading clears the model. The cold solution is saved and its three printed parameter tokens inserted into an otherwise unchanged author EOS at Tr=300 K before loading the hot data. This checks the GUI workflow at saved-file precision, not an unrounded in-memory handoff. Initial files and outputs are archived; no estimation step or new model family is used.",
        "residual_sign": "Pcalculated - Pobserved",
        "statistics": "RMS and pointwise residuals are native Peritheos evaluations of GUI-saved models, not RMS read from GUI output. Screenshots independently retain displayed fit convergence, coefficients, esds and maximum residual. Covariance is the actual GUI-saved conditional covariance; cold covariance is not propagated into thermal gamma errors.",
        "limitation": "The author's exact 20210609 executable is unavailable. This 20211014 GUI comparison supports the tested controls only; original exact GUI selections/controls and saved-vs-printed source-output differences remain unresolved.",
        "cases": cases,
        "artifact_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths
        },
        "published_residual_envelope_reproduced": False,
    }


def comparison(name):
    data = (
        read_data(FOLDER / "hot.dat")
        if name == "staged-weighted"
        else read_data(AUTHOR / "FeS6EOSfittot-SataOhfuji-300K.dat")
    )
    data = data[data[:, 2] > 300]
    initial = (
        FOLDER / "thermal-initial.eos"
        if name == "staged-weighted"
        else AUTHOR / "FeS6.eos"
    )
    source, _, _ = model(initial)
    refit, _, _ = model(FOLDER / (name + ".eos"))
    v, p, t = data[:, 4], data[:, 0], data[:, 2]
    grid = np.linspace(v.min(), v.max(), 160)
    temps = np.unique(t)
    curves = [
        (
            f"{int(temp)} K",
            grid,
            source.pressure(grid / 10, temp),
            refit.pressure(grid / 10, temp),
        )
        for temp in temps[[0, len(temps) // 2, -1]]
    ]
    return DiagnosticComparison(
        [
            panel(
                name + ": unchanged 146 hot rows",
                v,
                p,
                source.pressure(v / 10, t),
                refit.pressure(v / 10, t),
                "Volume (cm³/mol FeS)",
                "Pressure (GPa)",
                "GPa",
                swap_axes=True,
                sigma=data[:, 1],
                x_sigma=data[:, 5],
                curves=curves,
                observed_label="Author-selected hot observations",
            )
        ],
        "Native evaluation of GUI-saved coefficients; three labelled isotherms and every hot-row residual. No uncertainty band or propagation of cold/calibration covariance.",
        source_title="GUI cold curve / initial gamma"
        if name == "staged-weighted"
        else "Saved author EOS",
    )


def write_report(report, check=False):
    lines = [
        "# Morard FeS: actual GUI fitting comparison",
        "",
        "Audit date: 2026-10-19. "
        "**The distributed 20211014 GUI reproduces the tested console fits closely and still exceeds "
        "the ±3 GPa envelope.** The published coefficients are accepted unchanged on "
        "parameter agreement; the residual discrepancy is informational.",
        "",
        "[Original audit and author-input confirmation](fes-eosfit-direct.md) · "
        "[Archived model files, covariance, hashes and residuals](../data/fes-gui-2026-10-19/manifest.json)",
        "",
        "The unmodified official Mac GUI was run on the confirmed author-used inputs. "
        "The accepted 167-row selection excludes all 11 rejected quenched states. "
        "Full MGD is selected, q-compromise off, θ₀=417 K, n=2 and q=1 fixed. "
        "The exact author 20210609 executable remains unavailable; this later 2021 build was "
        "already pinned in the version audit. These results concern the tested controls, "
        "and do not identify every original GUI setting.",
        "",
        "## Results",
        "",
        "| GUI fit | γ₀ ± conditional esd | Hot-row RMS (GPa) | Hot max | Hot rows outside ±3 |",
        "|---|---:|---:|---:|---:|",
    ]
    for name in ("gamma-unweighted", "gamma-weighted", "staged-weighted"):
        c = report["cases"][name]
        p = c["refined_parameters"]["18"]
        v, e = value_error(p["value"], p["error"])
        s = c["refit"]["thermal"]
        lines.append(
            f"| {name} | {v} ± {e} | {number(s['rmse_gpa'], 3)} | {number(s['max_abs_gpa'], 3)} | {len(s['rows_exceeding_3_gpa'])}/146 |"
        )
    lines += [
        "",
        "The weighted 167-row GUI gamma differs from the console by 0.00001; "
        "the staged gamma differs by 0.00002. These differences are much smaller "
        "than the conditional parameter errors. Full output precision remains in the saved models and manifest.",
        "",
        "The GUI cold fit gives V₀ = 15.4 ± 0.5 cm³/mol, K₀ = 120 ± 30 GPa and K′ = 5.0 ± 0.7. "
        "The next GUI data load clears the active model, so the three saved cold tokens were "
        "inserted into a copy of the author EOS at Tr=300 K for the hot stage. This is a "
        "saved-precision handoff; the archived initial files show exactly what was loaded. "
        "No estimation or q-compromise switch was used. Gamma-only controls use the saved "
        "author cold model at Tr=298 K and retain all 167 rows; the staged hot fit uses 146 rows.",
        "",
        "RMS and rowwise residuals below are native evaluations of GUI-saved models. "
        "The GUI screenshots independently preserve convergence, selected weights, refined "
        "parameters, esds and maximum residuals. GUI and native maxima agree within 0.02 GPa "
        "at displayed precision. Parameter covariance is the actual saved GUI matrix. "
        "Gamma errors condition on fixed cold coefficients and omit cold/calibrant uncertainty. "
        "Curves show labelled isotherms, while residuals include each row's measured temperature.",
        "",
    ]
    for name in ("gamma-unweighted", "gamma-weighted", "staged-weighted"):
        rel = f"../data/fes-gui-2026-10-19/{name}-curves.svg"
        path = PAGE.parent / rel
        c = comparison(name)
        if check:
            if not current_curve(path.read_bytes(), c) or (
                PAGE.parent / dark_path(rel)
            ).read_bytes() != dark_svg(path.read_bytes()):
                raise ValueError("Stale GUI comparison plot")
        else:
            if not path.exists() or not current_curve(path.read_bytes(), c):
                path.write_bytes(curve_svg(c))
            (PAGE.parent / dark_path(rel)).write_bytes(dark_svg(path.read_bytes()))
        lines += [
            "### " + name,
            "",
            *themed_images(name + " curves and residuals", rel),
            "",
            f"[GUI fit output screenshot](../data/fes-gui-2026-10-19/{name}.png)",
            "",
        ]
    lines += [
        "[Cold-fit GUI output](../data/fes-gui-2026-10-19/cold-weighted.png)",
        "",
        "## Reproduction",
        "",
        "Repository-only `uv run python -m scripts.audit_fes_gui_replay` "
        "recomputes the saved-model pressure comparisons and verifies the archived inputs. "
        "Add `--check` to check the manifest, page and plots. It does not rerun the GUI. "
        "To repeat the actual fits in EosFit: open the selected data and initial EOS, enable "
        "only the listed refinement flags and weights, click Fit, accept the converged solution "
        "and save it. For the staged case, first fit the cold dataset with thermal model None "
        "and V0/K0/K′ free, then load the hot dataset and archived thermal initial EOS, fix all "
        "but gamma0, enable all P/V/T weights and fit full MGD.",
        "",
        "The GUI comparison is now completed for these controls. Remaining limits are the "
        "original exact build/workflow and the source-output discrepancy, not missing author-used "
        "input files. No further input-file request or replacement EOS is implied.",
        "",
    ]
    content = "\n".join(lines)
    if check:
        if PAGE.read_text(encoding="utf-8") != content:
            raise ValueError("Stale GUI audit page")
    else:
        PAGE.write_text(content, encoding="utf-8")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    path = FOLDER / "manifest.json"
    if args.check:
        from scripts.report_plot_cache import same_inputs

        if not same_inputs(json.loads(path.read_text(encoding="utf-8")), report):
            raise ValueError("Stale GUI replay manifest")
    else:
        path.write_text(
            json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    write_report(report, check=args.check)
