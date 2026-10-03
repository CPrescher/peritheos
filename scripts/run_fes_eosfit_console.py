#!/usr/bin/env python3
"""Run an external official EosFit7c executable on unchanged FeS observations.

Requires a runnable console (on macOS including its X11 dependencies/display).
This does not install software or emulate the EosFit refinement driver. Each
case runs the real executable in a fresh directory and retains inputs, macro,
stdout, EosFit log, saved EOS and its covariance. No coefficients are promoted.
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

from scripts.audit_morard_2026_fes import (
    CELL_PER_MOLAR,
)
from scripts.audit_morard_2026_fes import (
    DATA as THERMAL_DATA,
)
from scripts.audit_morard_2026_fes import (
    observations as thermal_observations,
)
from scripts.audit_sata_2010_fes_vi import DATA as COLD_DATA
from scripts.audit_sata_2010_fes_vi import VR, elastic_properties, pressure
from scripts.audit_sata_2010_fes_vi import observations as cold_observations


def data_lines(dataset):
    """All volume conversions are derived; no source observation is altered."""
    thermal = dataset != "cold"
    lines = [
        f"Title FeS audit {dataset} unchanged observations",
        "Pscale GPa",
        "Vscale cm^3/mol" if thermal else "Vscale A^3/cell",
        "Tscale K",
        "Format 1 P SIGP V SIGV T SIGT" if thermal else "Format 1 P SIGP V SIGV T",
    ]
    if thermal:
        d = thermal_observations()
        for i in range(len(d["source_row"])):
            values = [
                d["pressure_gpa"][i],
                d["sigma_pressure_gpa"][i],
                d["volume_angstrom3"][i] / CELL_PER_MOLAR,
                d["sigma_volume_angstrom3"][i] / CELL_PER_MOLAR,
                d["temperature_k"][i],
                d["sigma_temperature_k"][i],
            ]
            lines.append(" ".join(format(x, ".17g") for x in values))
    if dataset != "morard":
        for row in cold_observations():
            if not row["included"]:
                continue
            factor = Avogadro / 4e24 if thermal else 1
            values = [
                row["pressure_gpa"],
                row["error_pressure_gpa"],
                row["volume_angstrom3"] * factor,
                row["error_volume_angstrom3"] * factor,
                300,
            ]
            if thermal:
                values.append(0)  # 300 K cold isotherm; no quoted T error.
            lines.append(" ".join(format(x, ".17g") for x in values))
    return lines


def macro(dataset, weighted, folder, *, q_compromise=False):
    """Explicit MGD choice; theta=417, n=2, Tr=300 and full-model q=1."""
    cold = dataset == "cold"
    commands = [
        f"log {folder / 'run.log'}",
        f"read {folder / 'input.dat'}",
        "input",
        "pr",
        "2",
        "3",
    ]
    commands += ["98.96", "148", "4.53"] if cold else ["15.4", "115.5", "4.99"]
    if not cold:
        commands += [
            "th",
            "7",
            "300",
            "y" if q_compromise else "n",
            "15.4",
            "417",
            "2",
            "2.42",
        ]
        if not q_compromise:
            commands += ["1"]
    commands += ["pscale", "GPa", "vscale", "A^3/cell" if cold else "cm^3/mol", "x"]
    commands += [f"save {folder / 'initial.eos'}", "y", "list", "fit", "n"]
    if cold:
        commands += ["y", "y", "y"]
    else:
        commands += ["y", "y", "y"] if dataset == "combined" else ["n"] * 3
        commands += ["n", "y"]  # theta fixed, gamma free.
        if not q_compromise:
            commands += ["n"]  # q fixed only in the full model.
    commands += ["y" if weighted else "n"] * (2 if cold else 3)
    commands += ["", "y", "n", f"save {folder / 'fitted.eos'}", "y"]
    return "\n".join(commands) + "\n"


def parse_console(stdout, dataset, *, n_observations=None, parameter_names=None):
    """Extract final-cycle precision, not rounded Param lines in .eos files."""
    block = stdout.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
    block = block.split("INPUT <CR> TO CONTINUE")[0]
    fitted = {}
    for line in block.splitlines():
        match = re.match(
            r"\s*(V0|K0|Kp|ThMGD|Natom|Gamm0|q)\s+1\s+"
            r"([-\d.]+)\s+([-\d.]+)\s+([-\d.]+)",
            line,
        )
        if match:
            name, value, shift, error = match.groups()
            fitted[name] = {"value": float(value), "esd": float(error)}
    expected = {"V0", "K0", "Kp"} if dataset == "cold" else {"Gamm0"}
    if dataset == "combined":
        expected |= {"V0", "K0", "Kp"}
    if parameter_names is not None:
        expected = set(parameter_names)
    if set(fitted) != expected:
        raise ValueError(f"Unexpected refined parameters: {fitted}")
    statistics = re.findall(r"CHI\^2 \(WEIGHTED\)=\s*([\d.]+)", stdout)
    maxima = re.findall(r"MAXIMUM DELTA-PRESSURE=\s*([-\d.]+)", stdout)
    n = (
        n_observations
        if n_observations is not None
        else {"cold": 13, "morard": 146, "combined": 159, "staged": 159}[dataset]
    )

    def table(text, fitted):
        result = []
        count = (14 if dataset == "cold" else 15) - (0 if fitted else 1)
        for line in text.splitlines():
            fields = line.split()
            if len(fields) == count and fields[0].isdigit() and fields[1] == "1":
                index = 8 if dataset == "cold" else 9
                result.append(
                    {"row": int(fields[0]), "residual_gpa": -float(fields[index])}
                )
        if len(result) != n:
            raise ValueError(f"Expected {n} console residuals, found {len(result)}")
        residual = np.array([row["residual_gpa"] for row in result])
        return {
            "rows": result,
            "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "min_gpa": float(residual.min()),
            "max_gpa": float(residual.max()),
            "rows_exceeding_3_gpa": [
                row["row"] for row in result if abs(row["residual_gpa"]) > 3
            ],
            "sign": "Pcalculated-Pobserved",
            "precision": "As printed by EosFit, not full internal precision",
        }

    warnings = sorted(
        {line.strip() for line in stdout.splitlines() if "WARNING" in line}
    )
    return {
        "converged": "Least Squares converged" in stdout,
        "refined_parameters_final_cycle": fitted,
        "eosfit_weighted_chi2_per_dof": float(statistics[-1]) if statistics else None,
        "eosfit_max_delta_pressure_gpa": float(maxima[-1]),
        "eosfit_residual_sign": "Pobserved-Pcalculated",
        "warnings": warnings,
        "initial_source_replay": table(
            stdout.split("RESULTS FROM CYCLE")[0], fitted=False
        ),
        "fitted_replay": table(
            stdout.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1], fitted=True
        ),
    }


def saved_covariance(path, dataset):
    """Extract the program's fixed-width covariance, preserving its scaling."""
    text = path.read_text().split("Variance-Covariance matrix=")[1]
    matrix = []
    for line in text.splitlines():
        numbers = re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)
        if len(numbers) == 59:
            matrix.append([float(x) for x in numbers])
        if len(matrix) == 59:
            break
    if len(matrix) != 59:
        raise ValueError("Missing saved EosFit covariance")
    indices = [0, 1, 2] if dataset == "cold" else [17]
    if dataset == "combined":
        indices = [0, 1, 2, 17]
    return np.asarray(matrix)[np.ix_(indices, indices)]


def enrich_saved_model(case, folder, dataset):
    covariance = saved_covariance(folder / "fitted.eos", dataset)
    case["saved_refined_parameter_covariance"] = covariance.tolist()
    case["covariance_parameter_order"] = (
        ["V0", "K0", "Kp"] if dataset == "cold" else ["Gamm0"]
    )
    if dataset == "combined":
        case["covariance_parameter_order"] = ["V0", "K0", "Kp", "Gamm0"]
    case["covariance_qualification"] = (
        "The actual saved EosFit covariance, with the program's esd scaling. "
        "Conditional on the input errors, model and selected rows; no "
        "inter-study calibration covariance or physical validation."
    )
    if dataset == "cold":
        values = case["refined_parameters_final_cycle"]
        theta = np.array([values[k]["value"] for k in ("V0", "K0", "Kp")])

        def convert(x):
            v0, k0, kp = x
            kr, krp = elastic_properties(VR, [0, k0, kp], vr=v0)
            return np.array([pressure(VR, [0, k0, kp], vr=v0), kr, krp])

        jac = np.empty((3, 3))
        for j in range(3):
            h = abs(theta[j]) * 1e-5
            step = np.eye(3)[j] * h
            jac[:, j] = (convert(theta + step) - convert(theta - step)) / (2 * h)
        reference_covariance = jac @ covariance @ jac.T
        case["sata_fixed_reference_representation"] = {
            "Vr_angstrom3": VR,
            "parameter_order": ["Pr_gpa", "Kr_gpa", "Kr_prime"],
            "parameters": convert(theta).tolist(),
            "covariance": reference_covariance.tolist(),
            "esds": np.sqrt(np.diag(reference_covariance)).tolist(),
            "qualification": "Curve re-expression and local covariance propagation, not a second EosFit fit.",
        }


def run(executable, output, *, q_compromise=False):
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "source": "https://www.rossangel.com/Download/EoSFit7installer_MacOS.dmg",
        "version": None,
        "cold_source_sha256": hashlib.sha256(COLD_DATA.read_bytes()).hexdigest(),
        "thermal_source_sha256": hashlib.sha256(THERMAL_DATA.read_bytes()).hexdigest(),
        "qualification": (
            "Direct official executable refits with explicit diagnostic objectives. "
            "These are not the authors' recovered input files. The combined fit "
            "uses all 146 thermal rows and all 13 re-reported cold VI rows, with "
            "no independent recalibration or inter-study scale correction. "
            "EosFit's reported esds/covariance are retained in its outputs."
        ),
        "cases": {},
    }
    if q_compromise:
        report["thermal_model"] = (
            "MGD q-compromise: theta constant, gamma/V constant; q undefined"
        )
    for dataset in ("cold", "morard", "combined"):
        for weighted in (False, True):
            name = dataset + ("-errors" if weighted else "-unit")
            folder = output / name
            folder.mkdir()
            (folder / "input.dat").write_text("\n".join(data_lines(dataset)) + "\n")
            (folder / "run.mcr").write_text(
                macro(dataset, weighted, folder, q_compromise=q_compromise)
            )
            result = subprocess.run(
                [str(executable)],
                # The Mac build restores its last DISLIN working directory;
                # subprocess cwd alone does not determine its data directory.
                # Absolute filenames preserve the user's existing WDIR preference.
                input=f"macro {folder / 'run.mcr'}\nexit\n",
                text=True,
                capture_output=True,
                cwd=folder,
                timeout=120,
            )
            (folder / "stdout.txt").write_text(result.stdout)
            (folder / "stderr.txt").write_text(result.stderr)
            if result.returncode != 0 or not (folder / "fitted.eos").is_file():
                raise RuntimeError(f"EosFit failed: {name}; inspect {folder}")
            version = re.search(r"Version:\s*([\d.]+)", result.stdout)
            report["version"] = version.group(1) if version else "unknown"
            case = parse_console(result.stdout, dataset)
            enrich_saved_model(case, folder, dataset)
            if not case["converged"]:
                raise RuntimeError(f"EosFit did not converge: {name}")
            case["dataset_rows"] = {"cold": 13, "morard": 146, "combined": 159}[dataset]
            case["weights"] = (
                "all supplied P/V/T errors" if weighted else "unit weights"
            )
            if dataset == "cold" and weighted:
                case["weights"] = "supplied P/V errors"
            case["files_sha256"] = {
                file.name: hashlib.sha256(file.read_bytes()).hexdigest()
                for file in sorted(folder.iterdir())
            }
            report["cases"][name] = case
            print(name, case["refined_parameters_final_cycle"], flush=True)
    (output / "manifest.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


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
