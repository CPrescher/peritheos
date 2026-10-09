"""Replay the author-described three stages with the official EosFit7c console.

The user relayed an in-person explanation from Miozzi on 2026-10-03:
cold-only fit, thermal fit with cold coefficients fixed, then release the
coefficients starting from stage two. The user clarified that V0 is released
and theta0 stays fixed at 420 K. Exact weights remain unavailable.
The optional q-compromise tests either retain that approximation or switch
to full MGD after stage two, as allowed by the EosFit GUI manual.
No library coefficients are changed.
Requires an external executable and its working display/runtime environment.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from scripts.reconstruct_miozzi_2020_iron import (
    NA,
    bm3,
    constant_q_gamma_theta,
    debye_energy,
    pressure_column,
    read_data,
)
from scripts.refit_miozzi_2020_tange import target_pressures

FACTOR = NA / 2e24
PARAMETERS = {"V0": 0, "K0": 1, "Kp": 2, "ThMGD": 10, "Gamm0": 17, "q": 18}


def input_data(data, pressures, cold, molar_volume_multiplier=1):
    select = data["temperature_k"] <= 300 if cold else np.ones(len(pressures), bool)
    lines = [
        "Title Miozzi author-described staged reconstruction",
        "Pscale GPa",
        "Vscale cm^3/mol",
        "Tscale K",
        "Format 1 P SIGP V SIGV T" if cold else "Format 1 P SIGP V SIGV T SIGT",
    ]
    for i in np.flatnonzero(select):
        coordinates = [
            pressures[i],
            data["pressure_error_gpa"][i],
            data["volume_a3"][i] * FACTOR * molar_volume_multiplier,
            data["volume_error_a3"][i] * FACTOR * molar_volume_multiplier,
            data["temperature_k"][i],
        ]
        if not cold:
            coordinates.append(data["temperature_error_k"][i])
        # EosFit's input format requires numbers. Undefined errors are marked
        # zero only in this console input, never in the canonical source CSVs.
        lines.append(
            " ".join(f"{x if np.isfinite(x) else 0:.17g}" for x in coordinates)
        )
    return "\n".join(lines) + "\n"


def fit_commands(free, weighted, thermal, q_compromise=False):
    names = ["V0", "K0", "Kp"]
    if thermal:
        names += ["ThMGD", "Gamm0"] + ([] if q_compromise else ["q"])
    return [
        "fit",
        "n",
        *["y" if n in free else "n" for n in names],
        *["y" if weighted else "n"] * (3 if thermal else 2),
        "",
        "y",
        "n",
    ]


def macro(
    folder,
    weighted,
    variant,
    q_compromise=False,
    switch_to_full=False,
    atoms=1,
    molar_volume_multiplier=1,
):
    thermal_free = ["Gamm0"] + ([] if q_compromise else ["q"])
    final_free = ["V0", "K0", "Kp", *thermal_free]
    if variant == "theta_free":
        thermal_free += ["ThMGD"]
        final_free += ["ThMGD"]
    elif variant == "V0_theta_fixed":
        final_free.remove("V0")
    commands = [
        f"log {folder / 'run.log'}",
        f"read {folder / 'cold.dat'}",
        "input",
        "pr",
        "2",
        "3",
        f"{22.81 * FACTOR * molar_volume_multiplier:.17g}",
        "129",
        "6.24",
        "pscale",
        "GPa",
        "vscale",
        "cm^3/mol",
        "x",
    ]
    commands += fit_commands(["V0", "K0", "Kp"], weighted, False)
    commands += [
        f"save {folder / 'stage1.eos'}",
        "y",
        "clear",
        "y",
        f"read {folder / 'all.dat'}",
        "input",
        "th",
        "7",
        "300",
        "y" if q_compromise else "n",
        "",
        "420",
        str(atoms),
        "1.11",
    ]
    if not q_compromise:
        commands += ["0.3"]
    commands += ["x"]
    commands += fit_commands(thermal_free, weighted, True, q_compromise)
    commands += [f"save {folder / 'stage2.eos'}", "y"]
    if switch_to_full:
        # Retain the in-memory cold curve and gamma from stage two. Changing
        # the thermal model prompts for these again; blank means retain.
        commands += [
            "input",
            "th",
            "7",
            "300",
            "n",
            "",
            "420",
            str(atoms),
            "",
            "0.3",
            "x",
        ]
        q_compromise = False
        final_free += ["q"]
    commands += fit_commands(final_free, weighted, True, q_compromise)
    commands += [f"save {folder / 'stage3.eos'}", "y"]
    return "\n".join(commands) + "\n"


def covariance(path, names, molar_volume_multiplier=1):
    if not path.exists():
        return None
    text = path.read_text(encoding="utf-8").split("Variance-Covariance matrix=")[-1]
    matrix = []
    for line in text.splitlines():
        values = re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)
        if len(values) == 59:
            matrix.append([float(x) for x in values])
        if len(matrix) == 59:
            break
    if len(matrix) != 59:
        return None
    indices = [PARAMETERS[name] for name in names]
    block = np.array(matrix)[np.ix_(indices, indices)]
    conversion = np.array(
        [
            1 / (FACTOR * molar_volume_multiplier) if name == "V0" else 1
            for name in names
        ]
    )
    return (block * np.outer(conversion, conversion)).tolist()


def model_pressure(
    v, t, values, q_compromise=False, atoms=1, molar_volume_multiplier=1
):
    """Independent pressure check; q-compromise has no adjustable q."""
    p = bm3(v, values["V0"], values["K0"], values["Kp"])
    if "Gamm0" not in values:
        return p
    if q_compromise:
        return (
            p
            + values["Gamm0"]
            * (
                debye_energy(values["ThMGD"], t, atoms)
                - debye_energy(values["ThMGD"], 300, atoms)
            )
            / (values["V0"] * FACTOR * molar_volume_multiplier)
            / 1000
        )
    gamma, theta = constant_q_gamma_theta(
        v / values["V0"], values["Gamm0"], values["q"], values["ThMGD"]
    )
    return (
        p
        + gamma
        * (debye_energy(theta, t, atoms) - debye_energy(theta, 300, atoms))
        / (v * FACTOR * molar_volume_multiplier)
        / 1000
    )


def parse_stages(
    text,
    folder,
    data,
    pressures,
    q_compromise=False,
    switch_to_full=False,
    atoms=1,
    molar_volume_multiplier=1,
):
    blocks = text.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[1:]
    stages = []
    for i, block in enumerate(blocks[:3]):
        table = block.split("INPUT <CR> TO CONTINUE")[0]
        values, errors = {}, {}
        for line in table.splitlines():
            match = re.match(
                r"\s*(V0|K0|Kp|ThMGD|Gamm0|q)\s+([01])\s+"
                r"([-\d.]+)(?:\s+([-\d.]+)\s+([-\d.]+))?",
                line,
            )
            if not match:
                continue
            name, flag, value, _, error = match.groups()
            conversion = 1 / (FACTOR * molar_volume_multiplier) if name == "V0" else 1
            values[name] = float(value) * conversion
            if flag == "1" and error is not None:
                errors[name] = float(error) * conversion
        select = (
            data["temperature_k"] <= 300 if i == 0 else np.ones(len(pressures), bool)
        )
        p = model_pressure(
            data["volume_a3"][select],
            data["temperature_k"][select],
            values,
            q_compromise and not (switch_to_full and i == 2),
            atoms,
            molar_volume_multiplier,
        )
        residual = p - pressures[select]
        stage = {
            "parameters": values,
            "conditional_eosfit_standard_errors": errors,
            "rows": int(np.sum(select)),
            "free_parameters": list(errors),
            "covariance": covariance(
                folder / f"stage{i + 1}.eos", list(errors), molar_volume_multiplier
            ),
            "rmse_gpa_from_printed_coefficients": float(np.sqrt(np.mean(residual**2))),
            "residual_range_gpa": [float(residual.min()), float(residual.max())],
        }
        stages.append(stage)
    return stages


def run(
    executable,
    output,
    q_compromise=False,
    switch_to_full=False,
    atoms=1,
    molar_volume_multiplier=1,
):
    output.mkdir(parents=True, exist_ok=False)
    _, data, hashes = read_data()
    report = {
        "protocol_provenance": "User's 2026-10-03 account of an in-person conversation with Miozzi: RT-only; all data with RT coefficients fixed; release V0,K0,Kprime,gamma0,q starting from stage two. User explicitly clarified that theta0 stays fixed at 420 K and V0 is free in the final stage.",
        "qualification": f"Diagnostic replay, not recovered author inputs. Room temperature means all 36 rows at 298/300 K. n={atoms}, Tr=300 and theta0=420 K fixed. Unit weights in every stage. Original weights/handling of incomplete errors remain unspecified. Blanks are written as zero solely for unused console-input error columns; canonical observations are untouched. No calibration covariance or new physical validation.",
        "source_csv_sha256": hashes,
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "normalization": {
            "eosfit_atoms_per_formula_unit": atoms,
            "molar_volume_multiplier_relative_to_moles_of_Fe": molar_volume_multiplier,
            "volume_input_conversion": f"V_cell * N_A / (2e24) * {molar_volume_multiplier}",
            "reported_V0_units": "angstrom^3 per two-atom hcp Fe cell; conversion undone for comparison",
        },
        "thermal_model": (
            "q_compromise_then_full_mgd"
            if switch_to_full
            else "q_compromise_mgd"
            if q_compromise
            else "full_mgd"
        ),
        "cases": {},
    }
    if atoms == 2:
        report["protocol_provenance"] += (
            " User subsequently relayed Miozzi's confirmation that the author"
            " model was full MGD with n=2 and V0 approximately 6.87 cm3/mol."
            " The additional doubled-volume run is a normalization control,"
            " not the author-reported setup."
        )
        report["qualification"] += (
            " Doubling n without doubling molar volumes changes the physical"
            " thermal-pressure normalization. Doubling both is equivalent to"
            " n=1 with volumes per mole of Fe. The user confirmed that the"
            " author used n=2 with V0 approximately 6.87 cm3/mol; this"
            " provenance is an in-person report, not an original input file."
        )
    if q_compromise:
        report["qualification"] += (
            " Q-compromise keeps theta(V)=420 K and gamma(V)/V=gamma0/V0;"
            " q is absent, not a free or fitted parameter. Stage two releases"
            " gamma0 only."
        )
    if switch_to_full:
        report["qualification"] += (
            " In this variant only stage two uses q-compromise; before stage"
            " three switch to full MGD, retaining stage-two cold coefficients"
            " and gamma0 and setting the initial q to 0.3. Stage three releases"
            " V0,K0,Kprime,gamma0,q, with theta0 still fixed at 420 K."
        )
    elif q_compromise:
        report["qualification"] += (
            " Stage three retains q-compromise and releases V0,K0,Kprime,gamma0."
            " This tests the GUI approximation, not the paper's reported"
            " variable-q model."
        )
    for calibration in ["printed", "tange_vinet", "speziale_variable_q_debye"]:
        pressures = data["pressure_gpa"].copy()
        if calibration == "tange_vinet":
            pressures = target_pressures(data, "vinet")
        elif calibration != "printed":
            pressures = pressure_column(data, calibration)[0]
        for weighted in [False]:
            for variant in ["theta_fixed"]:
                name = f"{calibration}-{'errors' if weighted else 'unit'}-{variant}"
                folder = output / name
                folder.mkdir()
                (folder / "cold.dat").write_text(
                    input_data(data, pressures, True, molar_volume_multiplier),
                    encoding="utf-8",
                )
                (folder / "all.dat").write_text(
                    input_data(data, pressures, False, molar_volume_multiplier),
                    encoding="utf-8",
                )
                (folder / "run.mcr").write_text(
                    macro(
                        folder,
                        weighted,
                        variant,
                        q_compromise,
                        switch_to_full,
                        atoms,
                        molar_volume_multiplier,
                    ),
                    encoding="utf-8",
                )
                try:
                    result = subprocess.run(
                        [str(executable)],
                        input=f"macro {folder / 'run.mcr'}\nexit\n",
                        text=True,
                        capture_output=True,
                        cwd=folder,
                        timeout=50,
                    )
                    (folder / "stdout.txt").write_text(result.stdout, encoding="utf-8")
                    (folder / "stderr.txt").write_text(result.stderr, encoding="utf-8")
                    stages = parse_stages(
                        result.stdout,
                        folder,
                        data,
                        pressures,
                        q_compromise,
                        switch_to_full,
                        atoms,
                        molar_volume_multiplier,
                    )
                    case = {
                        "returncode": result.returncode,
                        "completed_stages": len(stages),
                        "converged_stages": result.stdout.count(
                            "Least Squares converged"
                        ),
                        "stages": stages,
                        "warnings": sorted(
                            {
                                line.strip()
                                for line in result.stdout.splitlines()
                                if "WARNING" in line
                            }
                        ),
                    }
                    case["stages_saved"] = all(
                        (folder / f"stage{i}.eos").exists() for i in range(1, 4)
                    )
                    if len(stages) != 3 or not case["stages_saved"]:
                        case["failure"] = (
                            "Three-stage protocol did not complete; inspect console output"
                        )
                except (
                    subprocess.TimeoutExpired,
                    ValueError,
                    KeyError,
                    FloatingPointError,
                ) as error:
                    case = {"failure": str(error)}
                case["files_sha256"] = {
                    p.name: hashlib.sha256(p.read_bytes()).hexdigest()
                    for p in sorted(folder.iterdir())
                }
                report["cases"][name] = case
                (output / "report.json").write_text(
                    json.dumps(report, indent=2, allow_nan=False) + "\n",
                    encoding="utf-8",
                )
                print(
                    name,
                    json.dumps(
                        {
                            k: v
                            for k, v in case.items()
                            if k not in ["files_sha256", "warnings", "stages"]
                        }
                    ),
                    flush=True,
                )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    model = parser.add_mutually_exclusive_group()
    model.add_argument("--q-compromise", action="store_true")
    model.add_argument("--q-compromise-start", action="store_true")
    parser.add_argument("--atoms", type=int, choices=[1, 2], default=1)
    parser.add_argument(
        "--molar-volume-multiplier", type=int, choices=[1, 2], default=1
    )
    args = parser.parse_args()
    run(
        args.executable.resolve(),
        args.output.resolve(),
        args.q_compromise or args.q_compromise_start,
        args.q_compromise_start,
        args.atoms,
        args.molar_volume_multiplier,
    )


if __name__ == "__main__":
    main()
