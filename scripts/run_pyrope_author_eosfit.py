"""Replay the supplied pyrope fit with the external official EosFit console.

Repository-only; requires an already runnable console and display/runtime.
Original attachments are copied unchanged into a private output directory.
No upstream input or catalog file is changed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from decimal import Decimal
from pathlib import Path

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_pyrope_angel_2022_observations import DATA, properties
from scripts.audit_pyrope_author_inputs import NAMES, PREFIX, read_log
from scripts.refit_pyrope_angel_2022_subset import GAS_CONSTANT


def final_refinement(text):
    block = text.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
    parameters, errors = {}, {}
    for label, key in zip(["K0", "Kp", "ThMGD", "Gamm0", "Sca 3"], NAMES):
        match = re.search(
            r"^\s*" + label + r"\s+1\s+(\S+)\s+(\S+)\s+(\S+)", block, re.M
        )
        parameters[key], errors[key] = float(match[1]), float(match[3])
    covariance_text = block.split("VARIANCE-COVARIANCE MATRIX RESCALED")[1].split(
        "UPDATE EOS"
    )[0]
    covariance = [
        [float(v) for v in line.split(":")[1].split()]
        for line in covariance_text.splitlines()
        if ":" in line
    ]
    chi = float(re.search(r"CHI\^2 \(WEIGHTED\)=\s*(\S+)", block)[1])
    return dict(
        parameters=parameters,
        residual_scaled_esds=errors,
        reduced_chi_square=chi,
        printed_rescaled_covariance=covariance,
        converged="Least Squares converged" in text,
    )


def prepare(directory, output):
    output.mkdir(parents=True, exist_ok=True)
    for name, target in [
        ("pyrope_qcompMGD_March2023.eos", "initial.eos"),
        ("Pyrope_PVT_used.dat", "pvt.dat"),
        ("Pyrope_PKT_used.dat", "pkt.dat"),
    ]:
        shutil.copy2(directory / name, output / target)
    commands = [
        f"log {output}/run.log",
        f"load {output}/initial.eos",
        f"read {output}/pvt.dat",
        "v",
        f"read {output}/pkt.dat",
        "a",
        "v",
        "group",
        "3",
        "118,-132",
        "",
        "",
        "fit",
        "n",
        "n",
        "y",
        "y",
        "y",
        "y",
        "n",
        "n",
        "y",
        "y",
        "y",
        "y",
        "",
        "y",
        "n",
        f"save {output}/fitted.eos",
        "y",
        "exit",
    ]
    macro = "\n".join(commands) + "\n"
    (output / "commands.txt").write_text(macro)
    return macro


def execute(executable, output, macro):
    result = subprocess.run(
        [str(executable)],
        input=macro,
        text=True,
        capture_output=True,
        cwd=output,
        timeout=45,
    )
    (output / "stdout.txt").write_text(result.stdout)
    (output / "stderr.txt").write_text(result.stderr)
    if result.returncode:
        raise ValueError(f"EosFit exited with {result.returncode}")
    return result.stdout


def compare(directory, output, executable):
    stdout = (output / "stdout.txt").read_text()
    native = final_refinement(stdout)
    eos_text = (output / "fitted.eos").read_text()
    covariance_rows = [
        re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)
        for line in eos_text.split("Variance-Covariance matrix=")[1].splitlines()
    ]
    covariance_tokens = [row for row in covariance_rows if len(row) == 59]
    full_covariance = np.array([[float(x) for x in row] for row in covariance_tokens])
    indices = [1, 2, 10, 17, 52]
    native_covariance = full_covariance[np.ix_(indices, indices)]
    native["saved_rescaled_covariance"] = native_covariance.tolist()
    native["saved_covariance_widths"] = dict(
        zip(NAMES, np.sqrt(np.diag(native_covariance)).tolist())
    )
    rounding = {}
    for i, name in zip(indices, NAMES):
        variance = full_covariance[i, i]
        half_step = 0.5 * 10.0 ** Decimal(covariance_tokens[i][i]).as_tuple().exponent
        lower, upper = np.sqrt(variance - half_step), np.sqrt(variance + half_step)
        width = native["residual_scaled_esds"][name]
        rounding[name] = dict(
            saved_variance_half_step=half_step,
            width_lower=float(lower),
            width_upper=float(upper),
        )
        if not lower - 5.1e-6 <= width <= upper + 5.1e-6:
            raise ValueError(
                "Saved covariance width disagrees beyond printed precision"
            )
    native["saved_covariance_rounding_checks"] = rounding
    author_text, author_rows = read_log(directory / "pyrope_25Oct2021.log")
    author = final_refinement(author_text)
    # Only the final list is checked, after the definitive group assignments.
    final_list = stdout.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[-1]
    actual_rows = [
        line.split()
        for line in final_list.splitlines()
        if re.match(r"^\s*\d+\s+[01]\s+\d+\s+(V|KS)\s", line)
    ]
    if len(actual_rows) != 194 or any(
        a[:10] != b[:10] for a, b in zip(actual_rows, author_rows)
    ):
        raise ValueError("Final native observed rows/groups differ from author's log")
    audit = json.loads((DATA / (PREFIX + "-audit.json")).read_text())
    ours = audit["independent_fits"][0]
    rows = []
    for name in NAMES:
        reference = author["parameters"][name]
        esd = author["residual_scaled_esds"][name]
        rows.append(
            dict(
                parameter=name,
                author_value=reference,
                author_rescaled_esd=esd,
                eosfit_value=native["parameters"][name],
                eosfit_rescaled_esd=native["residual_scaled_esds"][name],
                independent_value=ours["parameters"][name],
                independent_rescaled_esd=ours["conditional_residual_scaled_widths"][
                    name
                ],
                eosfit_minus_author=native["parameters"][name] - reference,
                independent_minus_author=ours["parameters"][name] - reference,
                eosfit_difference_over_author_esd=abs(
                    native["parameters"][name] - reference
                )
                / esd,
                independent_difference_over_author_esd=abs(
                    ours["parameters"][name] - reference
                )
                / esd,
                independent_error_relative_to_eosfit=ours[
                    "conditional_residual_scaled_widths"
                ][name]
                / native["residual_scaled_esds"][name]
                - 1,
            )
        )
    widths = np.sqrt(np.diag(native_covariance))
    native_correlation = native_covariance / np.outer(widths, widths)
    independent_correlation = np.array(
        audit["independent_covariance_check"]["correlation"]
    )
    return dict(
        scope="Replay of supplied processed inputs, not recovery of upstream reduction decisions",
        author_log=author,
        eosfit_replay=native,
        independent_reduced_chi_square=ours["conditional_reduced_chi_square"],
        comparison=rows,
        independent_error_comparison=dict(
            maximum_relative_width_difference=max(
                abs(row["independent_error_relative_to_eosfit"]) for row in rows
            ),
            maximum_absolute_correlation_difference=float(
                np.max(np.abs(native_correlation - independent_correlation))
            ),
            native_theta_gamma_correlation=float(native_correlation[2, 3]),
            independent_theta_gamma_correlation=float(independent_correlation[2, 3]),
        ),
        final194_observation_error_and_group_tokens_match_author_log=True,
        eosfit_version=re.search(r"Version:\s*(\S+)", stdout)[1],
        executable_sha256=hashlib.sha256(executable.read_bytes()).hexdigest(),
        native_exact_final_log_value_match=all(
            r["eosfit_minus_author"] == 0 for r in rows
        ),
        error_qualification="All compared widths use residual-rescaled local covariance; shared data mean differences are not independent-sample statistical tests",
    )


def fixed_properties(directory, output, executable):
    target = output / "fixed-properties"
    target.mkdir(exist_ok=True)
    text = (directory / "pyrope_qcompMGD_March2023.eos").read_text()
    p, _ = source("pyrope")
    # These coefficients are exactly representable in binary; this is a local
    # proposal for Ross's suggested comparison, not a set already agreed with him.
    proposed = {1: 112.0, 2: 170.0, 3: 4.5, 11: 800.0, 13: 20.0, 18: 1.25}
    for index, value in proposed.items():
        text = re.sub(
            r"(?m)^(Param =\s*" + str(index) + r"\s+).*",
            lambda m: m[1] + str(value),
            text,
        )
        p[index] = value
    eos = target / "proposed.eos"
    eos.write_text(text)
    states = [
        (pressure, temperature)
        for pressure in [0.0, 10.0, 30.0]
        for temperature in [100.0, 298.0, 1000.0, 1700.0]
    ]
    macro = "\n".join(
        [f"load {eos}", "props", "n"]
        + [f"{a} {b}" for a, b in states]
        + ["", "exit", ""]
    )
    (target / "commands.txt").write_text(macro)
    stdout = execute(executable, target, macro)
    chunks = re.split(r"At\s+P\s*=.*?The parameter values are:", stdout)[1:]
    if len(chunks) != len(states):
        raise ValueError("Unexpected property output")
    rows = []
    tolerances = {"volume_cm3_mol": 1e-4, "kt_gpa": 0.001, "ks_gpa": 0.001}
    for (pressure, temperature), chunk in zip(states, chunks):
        expected = properties(p, pressure, temperature, GAS_CONSTANT)
        actual = {
            key: float(re.search(r"(?m)^\s*" + label + r"\s+(\S+)", chunk)[1])
            for key, label in [
                ("volume_cm3_mol", "V"),
                ("kt_gpa", "K"),
                ("ks_gpa", "Ks"),
            ]
        }
        delta = {key: actual[key] - expected[key] for key in actual}
        rows.append(
            dict(
                pressure_gpa=pressure,
                temperature_k=temperature,
                independent=expected,
                eosfit=actual,
                differences=delta,
            )
        )
    return dict(
        scope="Proposed exactly representable coefficient set, not yet agreed with author; local same-parameter property check only",
        parameters=proposed,
        reference_temperature_k=298.0,
        gas_constant=GAS_CONSTANT,
        tolerances=tolerances,
        tolerance_basis="Previously declared native printed-property/single-precision tolerances, unchanged",
        all_within_tolerance=all(
            abs(r["differences"][key]) <= tolerances[key]
            for r in rows
            for key in tolerances
        ),
        states=rows,
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument(
        "--reuse-fit",
        action="store_true",
        help="Parse existing private fit output rather than running it again",
    )
    args = parser.parse_args()
    directory = args.source_dir.resolve()
    output = args.output.resolve()
    executable = args.executable.resolve()
    if not args.reuse_fit:
        macro = prepare(directory, output)
        (output / "fitted.eos").unlink(missing_ok=True)
        execute(executable, output, macro)
        if not (output / "fitted.eos").exists():
            raise ValueError("No fresh fit saved")
    result = compare(directory, output, executable)
    result["fixed_parameter_property_check"] = fixed_properties(
        directory, output, executable
    )
    (DATA / (PREFIX + "-eosfit.json")).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["comparison"], indent=2))
    print(
        "Fixed property check",
        result["fixed_parameter_property_check"]["all_within_tolerance"],
    )
