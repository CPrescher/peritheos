#!/usr/bin/env python3
"""Fit cold FeS first in EosFit; then fix that curve and refine gamma on all rows.

Both stages run in the same unmodified EosFit7c process, so the cold solution
is carried forward in memory without rounding it through a saved EOS file.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

from scipy.constants import Avogadro

from scripts.run_fes_eosfit_console import (
    COLD_DATA,
    THERMAL_DATA,
    data_lines,
    enrich_saved_model,
    parse_console,
    saved_covariance,
)


def cold_molar_lines():
    lines = data_lines("cold")
    lines[2] = "Vscale cm^3/mol"
    for i in range(5, len(lines)):
        values = list(map(float, lines[i].split()))
        values[2] *= Avogadro / 4e24
        values[3] *= Avogadro / 4e24
        lines[i] = " ".join(format(x, ".17g") for x in values)
    return lines


def staged_macro(folder, weighted, *, q_compromise=False):
    commands = [
        f"log {folder / 'run.log'}",
        f"read {folder / 'cold.dat'}",
        "input",
        "pr",
        "2",
        "3",
        "14.89877624024",
        "148",
        "4.53",
        "pscale",
        "GPa",
        "vscale",
        "cm^3/mol",
        "x",
        f"save {folder / 'cold-initial.eos'}",
        "y",
        "list",
        "fit",
        "n",
        "y",
        "y",
        "y",
    ]
    commands += ["y" if weighted else "n"] * 2
    commands += [
        "",
        "y",
        "n",
        f"save {folder / 'cold-fitted.eos'}",
        "y",
        "clear",
        "y",  # Only the in-memory data are cleared; EOS retained.
        f"read {folder / 'all.dat'}",
        "input",
        "th",
        "7",
        "300",
        "y" if q_compromise else "n",
        "",
        "417",
        "2",
        "2.42",
    ]
    if not q_compromise:
        commands += ["1"]
    commands += [
        "x",
        f"save {folder / 'thermal-initial.eos'}",
        "y",
        "list",
        "fit",
        "n",
        "n",
        "n",
        "n",
        "n",
        "y",
    ]
    if not q_compromise:
        commands += ["n"]  # q fixed only in full MGD; gamma alone refined.
    commands += ["y" if weighted else "n"] * 3
    commands += ["", "y", "n", f"save {folder / 'fitted.eos'}", "y"]
    return "\n".join(commands) + "\n"


def run(executable, output, *, q_compromise=False):
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "cold_source_sha256": hashlib.sha256(COLD_DATA.read_bytes()).hexdigest(),
        "thermal_source_sha256": hashlib.sha256(THERMAL_DATA.read_bytes()).hexdigest(),
        "procedure": "Same-process cold three-parameter BM3 fit, then gamma-only full MGD fit on all 159 rows with the cold solution fixed.",
        "qualification": "Recovered 13 cold VI observations are the Sata 2010 re-report, not proven complete original author inputs. Theta=417 K, n=2, q=1 and Tr=300 K fixed. Gamma errors are conditional on the fitted cold curve, not propagated cold uncertainty. Cold temperature errors are unavailable; all-data weighted run retains the isotherm assumption and EosFit warning.",
        "cases": {},
    }
    if q_compromise:
        report["procedure"] = (
            "Same-process cold three-parameter BM3 fit, then gamma-only MGD q-compromise fit on all 159 rows with the cold solution fixed."
        )
        report["thermal_model"] = (
            "MGD q-compromise: theta constant, gamma/V constant; q undefined"
        )
        report["qualification"] = report["qualification"].replace(
            "Theta=417 K, n=2, q=1 and Tr=300 K fixed.",
            "Theta=417 K, n=2 and Tr=300 K fixed; q undefined in q-compromise.",
        )
    for weighted in (False, True):
        name = "errors" if weighted else "unit"
        folder = output / name
        folder.mkdir()
        (folder / "cold.dat").write_text("\n".join(cold_molar_lines()) + "\n")
        (folder / "all.dat").write_text("\n".join(data_lines("combined")) + "\n")
        (folder / "run.mcr").write_text(
            staged_macro(folder, weighted, q_compromise=q_compromise)
        )
        result = subprocess.run(
            [str(executable)],
            input=f"macro {folder / 'run.mcr'}\nexit\n",
            capture_output=True,
            text=True,
            cwd=folder,
            timeout=120,
        )
        (folder / "stdout.txt").write_text(result.stdout)
        (folder / "stderr.txt").write_text(result.stderr)
        if result.returncode or not (folder / "fitted.eos").is_file():
            raise RuntimeError(f"Staged EosFit failed; inspect {folder}")
        cold_text = result.stdout.split("EOSFIT-7.6>clear")[0]
        cold = parse_console(cold_text, "cold")
        cold["covariance_parameter_order"] = ["V0_cm3_per_mol_fes", "K0_gpa", "Kprime"]
        cold["saved_refined_parameter_covariance"] = saved_covariance(
            folder / "cold-fitted.eos", "cold"
        ).tolist()
        thermal_text = result.stdout.split(f"save {folder / 'thermal-initial.eos'}", 1)[
            1
        ]
        thermal = parse_console(thermal_text, "staged")
        enrich_saved_model(thermal, folder, "staged")
        if not cold["converged"] or not thermal["converged"]:
            raise RuntimeError(f"Stage did not converge: {folder}")
        thermal["fixed_cold_parameters"] = cold["refined_parameters_final_cycle"]
        for field in ("initial_source_replay", "fitted_replay"):
            rows = thermal[field]["rows"][:146]
            thermal[field]["thermal_rows_exceeding_3_gpa"] = [
                row["row"] for row in rows if abs(row["residual_gpa"]) > 3
            ]
            thermal[field]["thermal_rows_rmse_gpa"] = (
                sum(row["residual_gpa"] ** 2 for row in rows) / 146
            ) ** 0.5
        # The final-cycle fixed values must exactly preserve the cold solution.
        final = thermal_text.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
        for parameter, value in cold["refined_parameters_final_cycle"].items():
            match = re.search(
                rf"\b{parameter}\s+0\s+([-\d.]+)\s+\[NOT REFINED\]", final
            )
            if not match or float(match.group(1)) != value["value"]:
                raise RuntimeError(f"Cold parameter changed: {parameter}")
        case = {
            "cold_stage": cold,
            "thermal_stage": thermal,
            "files_sha256": {
                p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                for p in sorted(folder.iterdir())
            },
        }
        report["cases"][name] = case
        print(
            name,
            cold["refined_parameters_final_cycle"],
            thermal["refined_parameters_final_cycle"],
            flush=True,
        )
    (output / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--q-compromise", action="store_true")
    args = parser.parse_args()
    run(
        args.executable.resolve(), args.output.resolve(), q_compromise=args.q_compromise
    )


if __name__ == "__main__":
    main()
