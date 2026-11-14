"""Replay the constructed 92-channel pyrope subset in official EosFit7c.

Requires an external console and working display/runtime. No installations,
author project, missing observations or increased Hartwig weights supplied.
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
from scripts.refit_pyrope_angel_2022_subset import NAMES, prepared_rows


def prepare(output):
    output.mkdir(parents=True, exist_ok=True)
    rows = prepared_rows()
    for name in ("milani_vt", "hartwig_vt", "zou_pv", "zou_ks"):
        lines = [
            f"Title CONDITIONAL recovered pyrope subset {name}",
            "Pscale GPa",
            "Vscale cm3/mol",
            "Tscale K",
        ]
        if name == "zou_ks":
            lines += ["Datatype moduli adiabatic"]
        lines += ["Format 1 P sigP V sigV T sigT"]
        for row in rows:
            if row["dataset"] == name:
                # Zero is the console encoding for an omitted variance term;
                # raw row map retains None, not a claimed zero physical error.
                values = (
                    row["P"],
                    row["sP"] or 0,
                    row["Y"],
                    row["sY"],
                    row["T"],
                    row["sT"] or 0,
                )
                lines += [" ".join(f"{value:.17g}" for value in values)]
        (output / f"{name}.dat").write_text("\n".join(lines) + "\n")
    text = (DATA.parent / "sources/pyrope-official.eos").read_text()
    text = re.sub(
        r"(?m)^Title =.*",
        "Title =DIAGNOSTIC recovered subset; not original author fit",
        text,
    )
    text = re.sub(r"(?m)^Comment =.*\n", "", text)
    text = text.replace(
        "System =Cubic",
        "Comment =Constructed 92-channel subset; initial Hartwig errors only\nComment =Groups 1 Milani VT; 2 Hartwig VT; 3 Zou PV and KS\nComment =Source-used corrections and normalization anchors unrecovered\nSystem =Cubic",
    )
    for index, value, label in [
        (51, 1.001, "Milani VT"),
        (52, 1.0003, "Hartwig VT"),
        (53, 1, "Zou PV and KS; fixed"),
    ]:
        text = re.sub(
            r"(?m)^(Param =\s*" + str(index) + r"\s+).*",
            lambda m: m[1] + f"{value:<16.7f}(Sca {index - 50}, {label})",
            text,
        )
    (output / "initial.eos").write_text(text)
    (output / "row-map.json").write_text(json.dumps(rows, indent=2) + "\n")
    commands = [f"log {output / 'run.log'}", f"load {output / 'initial.eos'}"]
    for i, name in enumerate(("milani_vt", "hartwig_vt", "zou_pv", "zou_ks")):
        commands += [f"read {output / (name + '.dat')}"]
        if i:
            commands += ["a"]
    # Appending files automatically creates groups 1..4; consolidate Zou.
    commands += ["group", "3", "53,-92", "", "", "fit", "n"]
    commands += ["n", "y", "y", "y", "y", "y", "y", "n"]
    commands += ["y", "y", "y", "", "y", "n"]
    commands += [f"save {output / 'joint-fitted.eos'}", "y", "exit"]
    macro = "\n".join(commands) + "\n"
    (output / "commands.txt").write_text(macro)
    return macro


def run(executable, output):
    macro = prepare(output)
    process = subprocess.run(
        [str(executable)],
        input=macro,
        text=True,
        capture_output=True,
        cwd=output,
        timeout=60,
    )
    (output / "stdout.txt").write_text(process.stdout)
    (output / "stderr.txt").write_text(process.stderr)
    if process.returncode or not (output / "joint-fitted.eos").exists():
        raise ValueError(
            "Official console did not complete/save the constructed subset fit"
        )
    return compare(executable, output)


def compare(executable, output):
    """Compare saved output, retaining finite printed parameter precision."""
    stdout = (output / "stdout.txt").read_text()
    text = (output / "joint-fitted.eos").read_text()
    tokens = {
        int(i): value
        for i, value in re.findall(r"^Param =\s*(\d+)\s+(\S+)", text, re.M)
    }
    parameters = {i: float(value) for i, value in tokens.items()}
    reference = json.loads((DATA / "subset-refits.json").read_text())["cases"][
        "joint_recovered"
    ]["refit"]
    fit_names = NAMES[:4] + NAMES[4:]
    actual = dict(zip(fit_names, [parameters[i] for i in (2, 3, 11, 18, 51, 52)]))
    deltas = {key: actual[key] - reference["parameters"][key] for key in fit_names}
    fractions = {
        key: abs(deltas[key]) / reference["conditional_unscaled_widths"][key]
        for key in fit_names
    }
    rounding = dict(
        zip(
            fit_names,
            [
                0.5 * 10.0 ** Decimal(tokens[i]).as_tuple().exponent
                for i in (2, 3, 11, 18, 51, 52)
            ],
        )
    )
    chisq = float(re.findall(r"CHI\^2 \(WEIGHTED\)=\s*(\S+)", stdout)[-1])
    checked = {
        "scope": "Official console replay of a constructed recovered subset; not original Angel joint input or publication fit parity",
        "executable_source": "https://www.rossangel.com/Download/EoSFit7installer_MacOS.dmg",
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "version": re.search(r"Version:\s*(\S+)", stdout)[1],
        "n_channels": 92,
        "free_parameters": 6,
        "numerical_convergence_reported": "Least Squares converged" in stdout,
        "parameters": actual,
        "minus_independent_irls": deltas,
        "differences_over_conditional_unscaled_widths": fractions,
        "saved_parameter_half_last_digit": rounding,
        "conditional_reduced_chi_square": chisq,
        "independent_conditional_reduced_chi_square": reference[
            "conditional_reduced_chi_square"
        ],
        "comparison_tolerances": {
            "parameter_fraction_of_conditional_width": 0.01,
            "reduced_chi_square_absolute": 0.001,
        },
        "comparison_basis": "Nearby numerical fits checked at 1 percent of conditional marginal widths plus half the last printed parameter digit, consistent with console shift/esd termination, source single precision and file rounding. Scale factors are printed to only five decimals. This is not equality of algorithms/coefficients/covariance or original workflow parity.",
        "missing_error_encoding": "Some sigP/sigT are zero only to omit unavailable variance terms. Console warns; all total effective pressure variances are positive, independently checked. No physical zero error or uncertainty floor is inferred.",
    }
    checked["within_declared_comparison_tolerances"] = bool(
        all(
            abs(deltas[key])
            < 0.01 * reference["conditional_unscaled_widths"][key] + rounding[key]
            for key in fit_names
        )
        and abs(chisq - reference["conditional_reduced_chi_square"]) < 0.001
    )
    if not checked["within_declared_comparison_tolerances"]:
        raise ValueError(
            "Official/independent constructed-subset fits exceed comparison tolerances"
        )
    (output / "checks.json").write_text(json.dumps(checked, indent=2) + "\n")
    return checked


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, default=DATA / "subset-eosfit")
    args = parser.parse_args()
    run(args.executable.resolve(), args.output.resolve())
