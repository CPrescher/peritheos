"""Run the official external EosFit7c calculator on the two source EOS files.

Requires a runnable official console and its display/X11 dependencies. This
command installs nothing and performs no refinement. All inputs and stdout
are retained; rounded calculator results are explicitly calculated checks.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from scripts.audit_garnet_angel_2022 import FOLDER, NAMES, independent_pressure, source


def run(executable, output):
    output.mkdir(parents=True, exist_ok=True)
    report = {
        "scope": "official_calculator_parameterization_check_not_observations",
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "executable_source": "https://www.rossangel.com/Download/EoSFit7installer_MacOS.dmg",
        "pressure_output_decimal_places": 4,
        "absolute_tolerance_gpa": 0.0001,
        "eosfit_gas_constant": float(np.float32(8.314)),
        "tolerance_basis": "0.00005 GPa output rounding plus 0.00005 GPa allowance for the source calculator single-precision evaluation; gas constant separately harmonized. Default Peritheos retains its documented R.",
        "records": {},
    }
    for name in NAMES:
        p, _ = source(name)
        states = [
            (p[1] * ratio, t)
            for ratio in [0.85, 0.95, 1.0, 1.03]
            for t in [90.0, 298.0, 700.0, 1000.0, 1700.0]
        ]
        commands = [
            f"load {FOLDER / 'sources' / (name + '-official.eos')}",
            "pcal",
            "n",
            "n",
        ]
        commands += [f"{v:.12g} {t:.12g}" for v, t in states]
        commands += ["", "exit"]
        text = "\n".join(commands) + "\n"
        (output / f"{name}-commands.txt").write_text(text, encoding="utf-8")
        r = subprocess.run(
            [str(executable)],
            input=text,
            text=True,
            capture_output=True,
            cwd=output,
            timeout=60,
        )
        (output / f"{name}-stdout.txt").write_text(r.stdout, encoding="utf-8")
        (output / f"{name}-stderr.txt").write_text(r.stderr, encoding="utf-8")
        parsed = re.findall(
            r"PRESSURE AT volume=\s*(\S+) AND T =\s*(\S+)K IS\s*(\S+)", r.stdout
        )
        if r.returncode or len(parsed) != len(states):
            raise ValueError(
                f"{name}: calculator failed: exit {r.returncode}, states {len(parsed)}"
            )
        report["version"] = re.search(r"Version:\s*(\S+)", r.stdout)[1]
        rows = []
        for (v, t), (_, _, pressure) in zip(states, parsed):
            expected = independent_pressure(p, v, t)
            harmonized = independent_pressure(
                p, v, t, gas_constant=report["eosfit_gas_constant"]
            )
            actual = float(pressure)
            if abs(actual - harmonized) > report["absolute_tolerance_gpa"]:
                raise ValueError(
                    f"{name} non-parity at {v}, {t}: {actual} vs {expected}"
                )
            rows.append(
                {
                    "volume_cm3_mol": v,
                    "temperature_k": t,
                    "eosfit_pressure_gpa": actual,
                    "independent_pressure_gpa": expected,
                    "difference_gpa": actual - expected,
                    "source_gas_constant_pressure_gpa": harmonized,
                    "source_gas_constant_difference_gpa": actual - harmonized,
                }
            )
        report["records"][name] = rows
    (output / "checks.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=FOLDER / "eosfit")
    args = parser.parse_args()
    run(args.executable.resolve(), args.output.resolve())
