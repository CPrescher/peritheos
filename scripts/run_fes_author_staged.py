"""Compare published and reproduced cold curves with gamma-only EosFit fits.

Only the author-confirmed 167-row input is read. Its 21 literature cold rows
define stage one; the unchanged 146 hot rows define stage two. The cold
solution is carried between stages in the same official EosFit process.
"""

from __future__ import annotations

import argparse
import json
import os
import re
from pathlib import Path

import numpy as np

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from scripts.audit_fes_author_inputs import digest, read_data, run_macro, summary
from scripts.run_fes_eosfit_console import (
    enrich_saved_model,
    parse_console,
    saved_covariance,
)

ROOT = Path(__file__).resolve().parents[1]
AUTHOR = ROOT / "docs/data/fes-author-eosfit/sources"
OUTPUT = ROOT / "docs/data/fes-author-staged"
PUBLISHED = {"V0": 15.40, "K0": 115.5, "Kp": 4.99}


def macro(folder, cold_source, weighted, tref):
    commands = [
        "input",
        f"load {AUTHOR / 'FeS6.eos'}",
        "th",
        "0",
    ]
    if cold_source == "published":
        commands += ["pr", "2", "3", "15.4", "115.5", "4.99"]
    commands += ["x", f"read {folder / 'cold.dat'}"]
    if cold_source == "reproduced":
        commands += ["fit", "n", "y", "y", "y"]
        commands += ["y" if weighted else "n"] * 3
        commands += ["", "y", "n"]
    commands += [
        f"save {folder / 'cold-fitted.eos'}",
        "y",
        "clear",
        "y",
        "input",
        "th",
        "7",
        str(tref),
        "n",
        "",
        "417",
        "2",
        "2.42",
        "1",
        "x",
        f"read {folder / 'thermal.dat'}",
        f"save {folder / 'thermal-initial.eos'}",
        "y",
        "fit",
        "n",
        "n",
        "n",
        "n",
        "n",
        "y",
        "n",
    ]
    commands += ["y" if weighted else "n"] * 3
    commands += ["", "y", "n", f"save {folder / 'fitted.eos'}", "y", "exit"]
    return commands


def verify_native(report):
    """Check pressure replay against the independent native implementation."""
    data = read_data(AUTHOR / report["source_file"])
    hot = data[data[:, 2] > 300]
    for case in report["cases"].values():
        fixed = case["fixed_cold_parameters"]
        gamma = case["thermal_stage"]["refined_parameters_final_cycle"]["Gamm0"][
            "value"
        ]
        model = MieGruneisenDebye(
            BM3(fixed["V0"] / 10, fixed["K0"], fixed["Kp"]),
            case["reference_temperature_k"],
            417,
            gamma,
            1,
            2,
        )
        native = model.pressure(hot[:, 4] / 10, hot[:, 2]) - hot[:, 0]
        external = np.array(
            [
                row["residual_gpa"]
                for row in case["thermal_stage"]["fitted_replay"]["rows"]
            ]
        )
        difference = float(np.max(np.abs(native - external)))
        assert difference < 0.0012
        case["native_crosscheck"] = {
            "implementation": "Peritheos native full BM3-MGD",
            "max_residual_difference_gpa": difference,
            "thermal_residual_summary": summary(native, np.arange(146)),
            "qualification": "Uses EosFit's final printed coefficients; differences include console parameter/data rounding. This is numerical agreement, not independent physical validation.",
        }
    return report


def run(executable, output):
    output.mkdir(parents=True, exist_ok=False)
    source = AUTHOR / "FeS6EOSfittot-SataOhfuji-300K.dat"
    data = read_data(source)
    assert data.shape == (167, 6)
    assert not np.any((data[:, 2] == 300) & (data[:, 1] == 0.0056))
    cold, hot = data[data[:, 2] == 300], data[data[:, 2] > 300]
    assert len(cold) == 21 and len(hot) == 146
    lines = source.read_text(encoding="utf-8").splitlines()
    cold_text = (
        lines[0]
        + "\n"
        + "\n".join(
            line for line in lines[1:] if line.strip() and float(line.split()[2]) == 300
        )
        + "\n"
    )
    hot_text = (
        lines[0]
        + "\n"
        + "\n".join(
            line for line in lines[1:] if line.strip() and float(line.split()[2]) > 300
        )
        + "\n"
    )
    report = {
        "audit_date": "2026-10-03",
        "source_file": source.name,
        "source_sha256": digest(source),
        "executable_sha256": digest(executable),
        "selection": {"literature_cold": 21, "thermal": 146, "quenched": 0},
        "selection_provenance": "User reports Guillaume Morard's personal communication confirming exclusion of the 11 unsuitable quenched observations.",
        "fixed_thermal_parameters": {
            "theta0_k": 417,
            "n": 2,
            "q": 1,
            "q_compromise": False,
        },
        "published_thermal_parameter": {"gamma0": 2.42, "esd": 0.03},
        "qualification": "Full BM3-MGD. Thermal objective uses only the 146 unchanged hot observations. Gamma errors are conditional on each fixed cold curve, without propagation of its covariance. Tr=300 K is consistent with the cold isotherm; Tr=298 K is a saved-author-file control. Both are reported, not silently interchanged. Audit results do not replace published database coefficients.",
        "cases": {},
    }
    for tref in (300, 298):
        for weighted in (False, True):
            if tref == 298 and not weighted:
                continue
            for cold_source in ("published", "reproduced"):
                name = f"{cold_source}-{'errors' if weighted else 'unit'}-tr{tref}"
                folder = output / name
                # run_macro creates the case directory; inputs are prepared in
                # an adjacent directory so its no-overwrite guard remains useful.
                input_folder = output / (name + "-inputs")
                input_folder.mkdir()
                (input_folder / "cold.dat").write_text(cold_text, encoding="utf-8")
                (input_folder / "thermal.dat").write_text(hot_text, encoding="utf-8")
                commands = macro(folder, cold_source, weighted, tref)
                commands = [
                    c.replace(
                        str(folder / "cold.dat"), str(input_folder / "cold.dat")
                    ).replace(
                        str(folder / "thermal.dat"), str(input_folder / "thermal.dat")
                    )
                    for c in commands
                ]
                stdout = run_macro(executable, folder, commands)
                first, thermal_text = stdout.split("EOSFIT-7.6>clear", 1)
                cold_case = None
                if cold_source == "reproduced":
                    cold_case = parse_console(
                        first,
                        "combined",
                        n_observations=21,
                        parameter_names=["V0", "K0", "Kp"],
                    )
                    assert cold_case["converged"]
                    cold_case["covariance_parameter_order"] = ["V0", "K0", "Kp"]
                    cold_case["saved_refined_parameter_covariance"] = saved_covariance(
                        folder / "cold-fitted.eos", "cold"
                    ).tolist()
                    fixed = {
                        k: v["value"]
                        for k, v in cold_case["refined_parameters_final_cycle"].items()
                    }
                else:
                    fixed = PUBLISHED
                thermal = parse_console(thermal_text, "morard", n_observations=146)
                assert thermal["converged"]
                enrich_saved_model(thermal, folder, "morard")
                final = thermal_text.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
                for key, value in fixed.items():
                    found = re.search(rf"^\s*{key}\s+0\s+([-\d.]+)", final, re.M)
                    assert found and abs(float(found[1]) - value) < 1e-5
                residual = np.array(
                    [r["residual_gpa"] for r in thermal["fitted_replay"]["rows"]]
                )
                report["cases"][name] = {
                    "cold_source": cold_source,
                    "reference_temperature_k": tref,
                    "weighting": "supplied P/V/T uncertainties"
                    if weighted
                    else "equal pressure weights",
                    "fixed_cold_parameters": fixed,
                    "cold_stage": cold_case,
                    "thermal_stage": thermal,
                    "thermal_residual_summary": summary(residual, np.arange(146)),
                    "artifact_sha256": {
                        p.name: digest(p) for p in sorted(folder.iterdir())
                    },
                    "input_sha256": {
                        p.name: digest(p) for p in sorted(input_folder.iterdir())
                    },
                }
                print(
                    name,
                    fixed,
                    thermal["refined_parameters_final_cycle"],
                    report["cases"][name]["thermal_residual_summary"],
                    flush=True,
                )
    verify_native(report)
    (output / "manifest.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    parser.add_argument(
        "--library-path",
        type=Path,
        help="Optional macOS X11 library directory for a portable EosFit installation",
    )
    args = parser.parse_args()
    if args.library_path:
        os.environ["DYLD_LIBRARY_PATH"] = str(args.library_path.resolve())
    run(args.executable.resolve(), args.output.resolve())
