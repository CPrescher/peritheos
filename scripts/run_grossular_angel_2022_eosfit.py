"""Official same-coefficient grossular property checks (repository-only).

Run with --executable /path/to/eosfit7c in a working EosFit/X11 environment.
This performs no least-squares refinement and changes no source parameters.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

from scripts.audit_garnet_angel_2022 import FOLDER, source
from scripts.audit_grossular_angel_2022_observations import DATA, GAS_CONSTANT
from scripts.audit_pyrope_angel_2022_observations import properties

STATES = [(0.0001, t) for t in (20, 100, 298, 700, 1000, 1350)] + [
    (7.674, 298),
    (9.17, 300),
    (16.91, 1650),
]
TOLERANCES = {
    "volume_cm3_mol": 0.0001,
    "kt_gpa": 0.001,
    "ks_gpa": 0.001,
    "alpha_per_k": 2e-9,
    "gamma": 0.0001,
    "theta_k": 0.001,
}
PATTERNS = dict(
    volume_cm3_mol="V",
    kt_gpa="K",
    ks_gpa="Ks",
    alpha_per_k="alpha",
    gamma="Gamma",
    theta_k="DebT",
)


def check_output(stdout):
    chunks = re.split(r"At\s+P\s*=.*?The parameter values are:", stdout)[1:]
    if len(chunks) != len(STATES):
        raise ValueError("Unexpected official grossular PROPS output")
    p, _ = source("grossular")
    checked = []
    for (pressure, temperature), chunk in zip(STATES, chunks):
        expected = properties(p, pressure, temperature, GAS_CONSTANT) | {
            "theta_k": p[11]
        }
        actual = {
            key: float(re.search(r"(?m)^\s*" + label + r"\s+(\S+)", chunk)[1])
            for key, label in PATTERNS.items()
        }
        actual["alpha_per_k"] /= 1e5
        delta = {key: actual[key] - expected[key] for key in actual}
        checked.append(
            dict(
                pressure_gpa=pressure,
                temperature_k=temperature,
                independent=expected,
                official=actual,
                official_minus_independent=delta,
                all_within_tolerances=all(abs(delta[k]) < TOLERANCES[k] for k in delta),
            )
        )
    return dict(
        scope="same_coefficient_property_evaluation_not_fitting_workflow",
        gas_constant=GAS_CONSTANT,
        absolute_tolerances=TOLERANCES,
        tolerance_basis="Previously declared EosFit printed property precision and source single-precision allowances, unchanged from the pyrope property audit.",
        states=checked,
        all_within_tolerances=all(s["all_within_tolerances"] for s in checked),
    )


def run(executable, output):
    output.mkdir(parents=True, exist_ok=True)
    path = (FOLDER / "sources/grossular-official.eos").resolve()
    macro = "\n".join(
        [f"load {path}", "props", "n"]
        + [f"{p} {t}" for p, t in STATES]
        + ["", "exit", ""]
    )
    (output / "commands.txt").write_text(macro)
    process = subprocess.run(
        [str(executable)],
        input=macro,
        text=True,
        capture_output=True,
        cwd=output,
        timeout=45,
    )
    (output / "stdout.txt").write_text(process.stdout)
    (output / "stderr.txt").write_text(process.stderr)
    if process.returncode:
        raise ValueError("Official EosFit returned failure")
    report = check_output(process.stdout) | {
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "source_eos_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
    }
    (output / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DATA / "eosfit-properties")
    args = parser.parse_args()
    print(
        run(args.executable.resolve(), args.output.resolve())["all_within_tolerances"]
    )
