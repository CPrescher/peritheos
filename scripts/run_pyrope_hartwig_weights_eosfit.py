"""Official-console checks of two user-requested Hartwig sensitivity cases."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from decimal import Decimal
from pathlib import Path

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_pyrope_angel_2022_observations import DATA, properties
from scripts.audit_pyrope_hartwig_weights import weighted_rows
from scripts.refit_pyrope_angel_2022_expanded import EXPANDED_NAMES, expanded_rows
from scripts.refit_pyrope_angel_2022_subset import GAS_CONSTANT
from scripts.run_pyrope_expanded_eosfit import prepare


def run(executable, output):
    sensitivity = json.loads((DATA / "hartwig-weight-sensitivity.json").read_text())
    reports = {}
    for w in (500, 1000):
        target = output / str(w)
        rows = weighted_rows(expanded_rows(True), w)
        macro = prepare(target, True, rows=rows)
        initial = target / "initial.eos"
        initial.write_text(
            initial.read_text().replace(
                "DIAGNOSTIC expanded 191 channels",
                f"DIAGNOSTIC Hartwig weight {w}, 191 channels",
            )
        )
        for name in ("joint-fitted.eos", "run.log"):
            (target / name).unlink(missing_ok=True)
        process = subprocess.run(
            [str(executable)],
            input=macro,
            text=True,
            capture_output=True,
            cwd=target,
            timeout=60,
        )
        (target / "stdout.txt").write_text(process.stdout)
        (target / "stderr.txt").write_text(process.stderr)
        saved = target / "joint-fitted.eos"
        if process.returncode or not saved.exists():
            raise ValueError("Official sensitivity replay did not save a fresh result")
        tokens = {
            int(i): v
            for i, v in re.findall(r"^Param =\s*(\d+)\s+(\S+)", saved.read_text(), re.M)
        }
        indices = [2, 3, 11, 18, 51, 52, 53]
        actual = dict(zip(EXPANDED_NAMES, [float(tokens[i]) for i in indices]))
        fit = next(
            v for v in sensitivity["cases"]["191"] if v["weight_multiplier"] == w
        )["fit_from_published_start"]
        delta = {key: actual[key] - fit["parameters"][key] for key in actual}
        rounding = dict(
            zip(
                EXPANDED_NAMES,
                [0.5 * 10.0 ** Decimal(tokens[i]).as_tuple().exponent for i in indices],
            )
        )
        chi = float(re.findall(r"CHI\^2 \(WEIGHTED\)=\s*(\S+)", process.stdout)[-1])
        checked = {
            "scope": "User-requested diagnostic Hartwig sensitivity, not author-workflow reproduction",
            "n_channels": 191,
            "Hartwig_inverse_variance_multiplier": w,
            "parameters": actual,
            "numerical_convergence_reported": "Least Squares converged"
            in process.stdout,
            "conditional_reduced_chi_square": chi,
            "minus_independent_irls": delta,
            "saved_half_last_digits": rounding,
            "comparison_tolerances": {
                "parameter_fraction_of_conditional_unscaled_width": 0.01,
                "plus_saved_half_last_digit": True,
                "absolute_reduced_chi_square": 0.001,
            },
            "within_declared_tolerances": all(
                abs(delta[k])
                < 0.01 * fit["conditional_unscaled_widths"][k] + rounding[k]
                for k in actual
            )
            and abs(chi - fit["conditional_reduced_chi_square"]) < 0.001,
            "differences_over_conditional_unscaled_widths": {
                k: abs(delta[k]) / fit["conditional_unscaled_widths"][k] for k in actual
            },
            "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            "scope_limit": "The earlier strict solver-comparison threshold is retained. Any failure remains a failure; source coefficient and workflow parity are not established by convergent diagnostic fits.",
        }
        (target / "checks.json").write_text(json.dumps(checked, indent=2) + "\n")
        reports[str(w)] = checked
        print(
            w,
            actual,
            "strict comparison",
            checked["within_declared_tolerances"],
            flush=True,
        )
    return reports


def property_check(executable, output):
    """Direct same-coefficient equation checks across Hartwig and modulus T/P."""
    output.mkdir(parents=True, exist_ok=True)
    states = [(0, t) for t in (103, 150, 297.5, 600, 650, 1073)] + [
        (10, 300),
        (10, 1000),
        (19.86, 1700),
    ]
    files = {
        "published": DATA.parent / "sources/pyrope-official.eos",
        "weight_500": output.parent / "500/joint-fitted.eos",
        "weight_1000": output.parent / "1000/joint-fitted.eos",
    }
    tolerances = {
        "volume_cm3_mol": 0.0001,
        "kt_gpa": 0.001,
        "ks_gpa": 0.001,
        "alpha_per_k": 2e-9,
        "gamma": 0.0001,
        "theta_k": 0.001,
    }
    patterns = {
        "volume_cm3_mol": "V",
        "kt_gpa": "K",
        "ks_gpa": "Ks",
        "alpha_per_k": "alpha",
        "gamma": "Gamma",
        "theta_k": "DebT",
    }
    cases = {}
    for name, path in files.items():
        commands = (
            [f"load {path}", "props", "n"]
            + [f"{p} {t}" for p, t in states]
            + ["", "exit"]
        )
        macro = "\n".join(commands) + "\n"
        (output / f"{name}-commands.txt").write_text(macro)
        process = subprocess.run(
            [str(executable)],
            input=macro,
            text=True,
            capture_output=True,
            cwd=output,
            timeout=45,
        )
        (output / f"{name}-stdout.txt").write_text(process.stdout)
        (output / f"{name}-stderr.txt").write_text(process.stderr)
        chunks = re.split(r"At\s+P\s*=.*?The parameter values are:", process.stdout)[1:]
        if process.returncode or len(chunks) != len(states):
            raise ValueError("Unexpected official PROPS output")
        p, _ = source("pyrope")
        p.update(
            {
                int(i): float(v)
                for i, v in re.findall(
                    r"^Param =\s*(\d+)\s+(\S+)", path.read_text(), re.M
                )
            }
        )
        checked = []
        for (pressure, temperature), chunk in zip(states, chunks):
            expected = properties(p, pressure, temperature, GAS_CONSTANT) | {
                "theta_k": p[11]
            }
            actual = {
                key: float(re.search(r"(?m)^\s*" + label + r"\s+(\S+)", chunk)[1])
                for key, label in patterns.items()
            }
            actual["alpha_per_k"] /= 1e5
            delta = {key: actual[key] - expected[key] for key in actual}
            checked.append(
                {
                    "pressure_gpa": pressure,
                    "temperature_k": temperature,
                    "independent": expected,
                    "official": actual,
                    "official_minus_independent": delta,
                    "all_within_tolerances": all(
                        abs(delta[key]) < tolerances[key] for key in delta
                    ),
                }
            )
        cases[name] = {
            "eos_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "states": checked,
        }
    report = {
        "scope": "Direct equation/property evaluation with identical coefficients, not author-workflow or fitting-solver reproduction",
        "gas_constant": GAS_CONSTANT,
        "n_coefficient_sets": 3,
        "n_states_per_coefficient_set": len(states),
        "absolute_tolerances": tolerances,
        "tolerance_basis": "Same previously used four-decimal property output and source single-precision allowances; theta0 permits 0.001 K for single-precision output at hundreds of kelvin.",
        "all_within_tolerances": all(
            s["all_within_tolerances"] for c in cases.values() for s in c["states"]
        ),
        "cases": cases,
    }
    (output / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
    print(
        "27 same-coefficient property states:",
        report["all_within_tolerances"],
        flush=True,
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DATA / "hartwig-weight-eosfit")
    args = parser.parse_args()
    run(args.executable.resolve(), args.output.resolve())
    property_check(args.executable.resolve(), args.output.resolve() / "equation-checks")
