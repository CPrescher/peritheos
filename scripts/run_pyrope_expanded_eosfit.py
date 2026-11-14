"""Replay the two expanded conditional pyrope fits in official EosFit7c.

Requires the external console and its runtime/display, installs nothing.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from decimal import Decimal
from pathlib import Path

from scripts.audit_pyrope_angel_2022_observations import DATA
from scripts.refit_pyrope_angel_2022_expanded import (
    EXPANDED_NAMES,
    expanded_rows,
    quartz_pressure,
)


def prepare(output, include_temperature, *, rows=None):
    output.mkdir(parents=True, exist_ok=True)
    if rows is None:
        rows = expanded_rows(include_temperature)
    groups = [
        [r for r in rows if r["dataset"] == "milani_vt"],
        [r for r in rows if r["dataset"] == "hartwig_vt"],
        [r for r in rows if r["dataset"] == "sinogeikin_2000_ks"],
        [r for r in rows if r["kind"] == "V" and r["scale_index"] is None],
        [r for r in rows if r["kind"] == "KS" and r["scale_index"] is None],
    ]
    ordered = [r for group in groups for r in group]
    assert len(ordered) == len(rows)
    commands = [f"log {output / 'run.log'}", f"load {output / 'initial.eos'}"]
    for i, group in enumerate(groups, 1):
        lines = [
            f"Title DIAGNOSTIC expanded pyrope group {i}",
            "Pscale GPa",
            "Vscale cm3/mol",
            "Tscale K",
        ]
        if group[0]["kind"] == "KS":
            lines += ["Datatype moduli adiabatic"]
        lines += ["Format 1 P sigP V sigV T sigT"]
        lines += [
            " ".join(
                f"{value:.17g}"
                for value in (
                    r["P"],
                    r["sP"] or 0,
                    r["Y"],
                    r["sY"],
                    r["T"],
                    r["sT"] or 0,
                )
            )
            for r in group
        ]
        path = output / f"group-{i}.dat"
        path.write_text("\n".join(lines) + "\n")
        commands += [f"read {path}"]
        if i > 1:
            commands += ["a"]
    # Both the V and KS unscaled files belong to fixed-scale group 4.
    commands += ["group", "4", f"63,-{len(rows)}", "", "", "fit", "n"]
    # V0 fixed, four EOS coefficients and scales 1..3 free, scale 4 fixed.
    commands += ["n", "y", "y", "y", "y", "y", "y", "y", "n"]
    commands += [
        "y",
        "y",
        "y",
        "",
        "y",
        "n",
        f"save {output / 'joint-fitted.eos'}",
        "y",
        "exit",
    ]
    text = (DATA.parent / "sources/pyrope-official.eos").read_text()
    text = re.sub(
        r"(?m)^Title =.*",
        f"Title =DIAGNOSTIC expanded {len(rows)} channels; not author project",
        text,
    )
    text = re.sub(r"(?m)^Comment =.*\n", "", text)
    for index, value in ((51, 1.001), (52, 1.0003), (53, 0.97828), (54, 1)):
        text = re.sub(
            r"(?m)^(Param =\s*" + str(index) + r"\s+).*",
            lambda m: m[1] + f"{value:<16.7f}(Conditional group scale)",
            text,
        )
    (output / "initial.eos").write_text(text)
    (output / "row-map.json").write_text(json.dumps(ordered, indent=2) + "\n")
    macro = "\n".join(commands) + "\n"
    (output / "commands.txt").write_text(macro)
    return macro


def run(executable, output):
    reference = json.loads((DATA / "expanded-refits.json").read_text())
    reports = {}
    for include_t in (False, True):
        count = 191 if include_t else 182
        target = output / str(count)
        macro = prepare(target, include_t)
        # A stale saved EOS must never be mistaken for this run's result.
        # These files belong solely to this constructed replay.
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
        if process.returncode or not (target / "joint-fitted.eos").exists():
            raise ValueError("Expanded EosFit replay did not save a fresh result")
        fitted = (target / "joint-fitted.eos").read_text()
        tokens = {
            int(i): v for i, v in re.findall(r"^Param =\s*(\d+)\s+(\S+)", fitted, re.M)
        }
        indices = [2, 3, 11, 18, 51, 52, 53]
        actual = dict(zip(EXPANDED_NAMES, [float(tokens[i]) for i in indices]))
        fit = reference["cases"][str(count)]["refit"]
        diffs = {key: actual[key] - fit["parameters"][key] for key in actual}
        rounding = dict(
            zip(
                EXPANDED_NAMES,
                [0.5 * 10.0 ** Decimal(tokens[i]).as_tuple().exponent for i in indices],
            )
        )
        chi = float(re.findall(r"CHI\^2 \(WEIGHTED\)=\s*(\S+)", process.stdout)[-1])
        checked = {
            "scope": "Official console check of constructed recovered inputs, not full 194-channel author-workflow reproduction",
            "n_channels": count,
            "free_parameters": 7,
            "parameters": actual,
            "numerical_convergence_reported": "Least Squares converged"
            in process.stdout,
            "minus_independent_irls": diffs,
            "saved_half_last_digits": rounding,
            "conditional_reduced_chi_square": chi,
            "comparison_tolerances": {
                "parameter_fraction_of_conditional_width": 0.01,
                "plus_saved_half_last_digit": True,
                "absolute_reduced_chi_square": 0.001,
            },
            "within_declared_tolerances": all(
                abs(diffs[key])
                < 0.01 * fit["conditional_unscaled_widths"][key] + rounding[key]
                for key in actual
            )
            and abs(chi - fit["conditional_reduced_chi_square"]) < 0.001,
            "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
            "differences_over_conditional_unscaled_widths": {
                key: abs(value) / fit["conditional_unscaled_widths"][key]
                for key, value in diffs.items()
            },
            "numerical_limit": "Strict prior comparison threshold retained and not passed. Available CFML EoS source Get_Pressure_K uses modulus residual stopping tolerance 0.001*target modulus; native residual/weight comparisons show small differences. Attribution to this stopping rule is plausible, not isolated proof.",
        }
        (target / "checks.json").write_text(json.dumps(checked, indent=2) + "\n")
        reports[str(count)] = checked
        print(count, actual, checked["within_declared_tolerances"], flush=True)
    return reports


def calibration_check(executable, source_dir, output):
    """Check inverse-old/evaluate-new mapping with native PCAL, not PSHIFT UI."""
    output.mkdir(parents=True, exist_ok=True)
    milani = [r for r in expanded_rows() if r["dataset"] == "milani_pv"]
    commands = []
    for name in ("quartz-1997.eos", "quartz-2016.eos"):
        commands += [f"load {source_dir / name}", "pcal", "n", "n"]
        commands += [
            f"{r['calculated_equivalent_quartz_volume_a3']:.12g} 298" for r in milani
        ]
        commands += [""]
    commands += ["exit"]
    macro = "\n".join(commands) + "\n"
    (output / "commands.txt").write_text(macro)
    process = subprocess.run(
        [str(executable)],
        input=macro,
        text=True,
        capture_output=True,
        cwd=output,
        timeout=30,
    )
    (output / "stdout.txt").write_text(process.stdout)
    (output / "stderr.txt").write_text(process.stderr)
    parsed = re.findall(r"PRESSURE AT volume=\s*(\S+) IS\s*(\S+)", process.stdout)
    if len(parsed) != 22:
        raise ValueError(f"Unexpected quartz PCAL output: {len(parsed)}")
    checks = []
    for updated in (False, True):
        for i, r in enumerate(milani):
            actual = float(parsed[i + 11 * updated][1])
            expected = quartz_pressure(
                r["calculated_equivalent_quartz_volume_a3"], updated=updated
            )
            checks.append(
                {
                    "source_row": r["row"],
                    "updated": updated,
                    "native_pressure_gpa": actual,
                    "independent_pressure_gpa": expected,
                    "difference_gpa": actual - expected,
                }
            )
    report = {
        "scope": "Independent quartz BM3/BM4 old-pressure inversion/new-pressure evaluation checked against native PCAL. PSHIFT itself opens an unavailable dialog in this headless runtime and was not replayed.",
        "absolute_tolerance_gpa": 0.0006,
        "tolerance_basis": "Three-decimal console rounding (0.0005 GPa) plus 0.0001 GPa single-precision allowance",
        "all_within_tolerance": all(abs(r["difference_gpa"]) < 0.0006 for r in checks),
        "checks": checks,
    }
    (output / "checks.json").write_text(json.dumps(report, indent=2) + "\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DATA / "expanded-eosfit")
    parser.add_argument("--quartz-source-dir", type=Path)
    args = parser.parse_args()
    run(args.executable.resolve(), args.output.resolve())
    if args.quartz_source_dir:
        calibration_check(
            args.executable.resolve(),
            args.quartz_source_dir.resolve(),
            args.output.resolve() / "quartz-calibration",
        )
