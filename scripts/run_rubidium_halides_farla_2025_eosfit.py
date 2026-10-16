#!/usr/bin/env python3
"""Refit the printed Farla 2025 observations with an official EosFit7c binary.

Requires an external runnable executable and its runtime/display dependencies.
Each case preserves the generated input, row mapping, macro, raw console output,
native log and saved models/covariance. No published coefficients are replaced.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro

from scripts.investigate_rubidium_halides_farla_2025 import AUTHOR_THETA, OLD_NA
from scripts.reproduce_rubidium_halides_farla_2025 import (
    DATA,
    ERRORS,
    PREFIX,
    SOURCE,
)
from scripts.reproduce_rubidium_halides_farla_2025 import read_rows as source_rows
from scripts.run_fes_eosfit_console import parse_console

PARAMETERS = ["K0", "Kp", "Gamm0", "q"]
CASES = {
    "unweighted_printed_theta": {"weighted": False},
    "weighted_printed_theta": {"weighted": True},
    "weighted_precise_theta": {"weighted": True, "precise_theta": True},
    "weighted_all_lvp_precise_theta": {
        "weighted": True,
        "precise_theta": True,
        "include_cold_lvp": True,
    },
    "weighted_floor_0p005_precise_theta": {
        "weighted": True,
        "precise_theta": True,
        "pressure_floor": 0.005,
    },
    "weighted_legacy_na_precise_theta": {
        "weighted": True,
        "precise_theta": True,
        "legacy_na": True,
    },
}


def prepare_input(
    compound,
    *,
    include_cold_lvp=False,
    pressure_floor=0.0,
    legacy_na=False,
    include_dac=True,
    **_,
):
    """Convert Z=1 cell volumes to molar units, retaining source row identity."""
    rows = [
        r
        for r in source_rows(compound, "observations")
        if r["volume_a3"]
        and (include_dac or r["source_table"] != "S1")
        and (
            include_cold_lvp
            or r["source_table"] == "S1"
            or float(r["temperature_k"]) > 300
        )
    ]
    na = OLD_NA if legacy_na else Avogadro
    lines = [
        f"Title Farla2025 {compound} printed observations",
        "Pscale GPa",
        "Vscale cm^3/mol",
        "Tscale K",
        "Format 1 P SIGP V SIGV T SIGT",
    ]
    mapping = []
    for index, row in enumerate(rows, 1):
        values = [
            float(row["pressure_gpa"]),
            max(float(row["pressure_uncertainty_gpa"] or 0), pressure_floor),
            float(row["volume_a3"]) * na * 1e-24,
            float(row["volume_uncertainty_a3"]) * na * 1e-24,
            float(row["temperature_k"]),
            float(row["temperature_uncertainty_k"] or 0),
        ]
        lines.append(" ".join(format(value, ".17g") for value in values))
        mapping.append(
            {
                "eosfit_row": index,
                "source_table": row["source_table"],
                "source_row": int(row["source_row"]),
                "run": row["run"],
                "temperature_k": values[4],
            }
        )
    return "\n".join(lines) + "\n", mapping


def macro(
    compound,
    folder,
    *,
    weighted,
    precise_theta=False,
    fixed_static=None,
    thermal_start=None,
    refine_static=True,
    **_,
):
    molar, theta, coefficients = SOURCE[compound]
    coefficients = list(coefficients)
    if fixed_static is not None:
        molar, coefficients[0], coefficients[1] = fixed_static
    if thermal_start is not None:
        coefficients[2:] = thermal_start
    if precise_theta:
        theta = AUTHOR_THETA[compound]
    commands = [
        f"log {folder / 'run.log'}",
        f"read {folder / 'input.dat'}",
        "input",
        "pr",
        "2",
        "3",
        str(molar),
        str(coefficients[0]),
        str(coefficients[1]),
        "th",
        "7",
        "300",
        "n",
        str(molar),
        str(theta),
        "2",
        str(coefficients[2]),
        str(coefficients[3]),
        "pscale",
        "GPa",
        "vscale",
        "cm^3/mol",
        "x",
        f"save {folder / 'initial.eos'}",
        "y",
        "list",
    ]
    # Repeat after updating once to check/polish the default stopping tolerance.
    for _ in range(2):
        flag = "y" if refine_static else "n"
        commands += ["fit", "n", "n", flag, flag, "n", "y", "y"]
        commands += ["y" if weighted else "n"] * 3
        commands += ["", "y", "n"]
    commands += [f"save {folder / 'fitted.eos'}", "y"]
    return "\n".join(commands) + "\n"


def covariance(path, indices=(1, 2, 17, 18)):
    text = path.read_text().split("Variance-Covariance matrix=")[1]
    full = []
    for line in text.splitlines():
        numbers = re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)
        if len(numbers) == 59:
            full.append([float(x) for x in numbers])
        if len(full) == 59:
            break
    if len(full) != 59:
        raise ValueError("Missing native covariance matrix")
    return np.array(full)[np.ix_(indices, indices)]


def compare(compound, fit_parameters):
    source = np.array(SOURCE[compound][2])
    fit = np.array([fit_parameters[p]["value"] for p in PARAMETERS])
    delta = fit - source
    return {
        "published_parameters": source.tolist(),
        "eosfit_parameters": fit.tolist(),
        "eosfit_minus_published": delta.tolist(),
        "relative_difference_percent": (100 * delta / source).tolist(),
        "difference_over_published_esd": (delta / ERRORS[compound]).tolist(),
        "maximum_difference_over_published_esd": float(
            max(abs(delta / ERRORS[compound]))
        ),
        "all_coefficients_round_to_publication": bool(np.all(abs(delta) < 0.005)),
    }


def run(executable, output):
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "audit_date": "2026-10-16",
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "source": "https://www.rossangel.com/Download/EoSFit7installer_MacOS.dmg",
        "version": None,
        "parameter_order": PARAMETERS,
        "source_observation_hashes": {
            c: hashlib.sha256(
                (DATA / f"{PREFIX}-{c}-observations.csv").read_bytes()
            ).hexdigest()
            for c in SOURCE
        },
        "author_contact": "none",
        "original_author_input_files": False,
        "original_author_fit_parity": False,
        "qualification": "Direct official executable refits of generated inputs from printed observations; conditional on row selection, fixed anchors and explicit uncertainty choices. Not recovered author fit files, calibration verification or physical validation.",
        "cases": {},
    }
    independent = json.loads(
        (DATA.parents[2] / "docs/data" / f"{PREFIX}-follow-up.json").read_text()
    )
    for compound in SOURCE:
        report["cases"][compound] = {}
        for name, options in CASES.items():
            folder = output / compound / name
            folder.mkdir(parents=True)
            data, mapping = prepare_input(compound, **options)
            (folder / "input.dat").write_text(data)
            (folder / "row-map.json").write_text(json.dumps(mapping, indent=2) + "\n")
            (folder / "run.mcr").write_text(macro(compound, folder, **options))
            process = subprocess.run(
                [str(executable)],
                input=f"macro {folder / 'run.mcr'}\nexit\n".encode(),
                capture_output=True,
                cwd=folder,
                timeout=60,
            )
            (folder / "stdout.txt").write_bytes(process.stdout)
            (folder / "stderr.txt").write_bytes(process.stderr)
            if process.returncode != 0 or not (folder / "fitted.eos").is_file():
                raise RuntimeError(
                    f"EosFit failed for {compound}/{name}; inspect {folder}"
                )
            stdout = process.stdout.decode("utf-8", errors="replace")
            parsed = parse_console(
                stdout,
                "rubidium",
                n_observations=len(mapping),
                parameter_names=PARAMETERS,
            )
            if not parsed["converged"]:
                raise RuntimeError(f"No convergence for {compound}/{name}")
            version = re.search(r"Version:\s*([\d.]+)", stdout)
            report["version"] = version.group(1) if version else "unknown"
            parsed["input_options"] = options
            parsed["observations"] = len(mapping)
            parsed["refinement_passes"] = 2
            parsed["convergence_messages"] = stdout.count("Least Squares converged")
            parsed["saved_parameter_covariance"] = covariance(
                folder / "fitted.eos"
            ).tolist()
            parsed["covariance_parameter_order"] = PARAMETERS
            parsed["covariance_qualification"] = (
                "Native covariance with EosFit's own chi-squared esd scaling; conditional on supplied errors, selection and fixed parameters."
            )
            parsed["published_comparison"] = compare(
                compound, parsed["refined_parameters_final_cycle"]
            )
            if name == "weighted_precise_theta":
                ref = independent["compounds"][compound]["iterated_weight_refits"][
                    "author_theta_and_eosfit_R"
                ]["fit_parameters"]
                parsed["independent_irls_parameters"] = ref
                parsed["eosfit_minus_independent_irls"] = (
                    np.array(parsed["published_comparison"]["eosfit_parameters"]) - ref
                ).tolist()
            parsed["files_sha256"] = {
                f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                for f in sorted(folder.iterdir())
            }
            report["cases"][compound][name] = parsed
            print(
                compound,
                name,
                parsed["published_comparison"]["eosfit_parameters"],
                flush=True,
            )
            (output / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.executable.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
