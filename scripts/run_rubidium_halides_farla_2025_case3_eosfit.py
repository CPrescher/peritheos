#!/usr/bin/env python3
"""Direct EosFit7c Case 3 fits: fixed static anchors, free gamma0 and q.

The source includes two distinct RbBr Case 3 variants. No author fit files are
claimed and no source observations or canonical coefficients are modified.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np

from scripts.investigate_rubidium_halides_farla_2025 import AUTHOR_THETA
from scripts.reproduce_rubidium_halides_farla_2025 import DATA, PREFIX, SOURCE
from scripts.run_fes_eosfit_console import parse_console
from scripts.run_rubidium_halides_farla_2025_eosfit import (
    covariance,
    macro,
    prepare_input,
)

INPUTS = DATA / f"{PREFIX}-case3-inputs.json"
FACTS = json.loads(INPUTS.read_text(encoding="utf-8"))
VARIANTS = FACTS["variants"]
PARAMETERS = ["Gamm0", "q"]
CASES = {
    "unweighted_printed_theta": {"weighted": False},
    "weighted_printed_theta": {"weighted": True},
    "weighted_precise_theta": {"weighted": True, "precise_theta": True},
    "weighted_all_lvp_printed_theta": {"weighted": True, "include_cold_lvp": True},
    "weighted_lvp_only_printed_theta": {"weighted": True, "include_dac": False},
    "weighted_anchor_rounding_example": {
        "weighted": True,
        "rounding_anchor_example": True,
    },
}
# Illustrative fitted-DAC anchors inside their printed half-digit intervals.
# The deliberately selected RbBr Kprime=4.98 remains exactly 4.98; its example
# varies fitted V0 instead. These are not recovered author values.
ROUNDING_ANCHOR_EXAMPLES = {
    "rbcl": [36.25, 19.93, 4.989],
    "rbbr_dac": [41.74, 16.32, 5.2575],
    "rbbr_kprime_4p98": [41.549, 17.55, 4.98],
    "rbi": [51.11, 13.68, 4.9606],
}


def commands(variant, folder, options):
    source = VARIANTS[variant]
    anchors = list(source["fixed_parameters"][:3])
    if options.get("rounding_anchor_example"):
        anchors = ROUNDING_ANCHOR_EXAMPLES[variant]
    return macro(
        source["compound"],
        folder,
        **options,
        fixed_static=anchors,
        thermal_start=source["published_thermal_parameters"],
        refine_static=False,
    )


def comparison(variant, fitted):
    source = VARIANTS[variant]
    values = np.array([fitted[p]["value"] for p in PARAMETERS])
    delta = values - source["published_thermal_parameters"]
    scaled = delta / source["published_thermal_esd"]
    return {
        "published_parameters": source["published_thermal_parameters"],
        "published_esd": source["published_thermal_esd"],
        "eosfit_parameters": values.tolist(),
        "eosfit_minus_published": delta.tolist(),
        "difference_over_published_esd": scaled.tolist(),
        "maximum_absolute_difference_over_published_esd": float(max(abs(scaled))),
        "both_parameters_round_to_publication": bool(np.all(abs(delta) < 0.005)),
    }


def run(executable, output):
    output.mkdir(parents=True, exist_ok=False)
    report = {
        "audit_date": "2026-10-16",
        "source_doi": FACTS["source_doi"],
        "source_facts_sha256": hashlib.sha256(INPUTS.read_bytes()).hexdigest(),
        "executable_sha256": hashlib.sha256(executable.read_bytes()).hexdigest(),
        "version": None,
        "parameter_order": PARAMETERS,
        "fixed_parameter_order": FACTS["fixed_parameter_order"],
        "source_observation_hashes": {
            c: hashlib.sha256(
                (DATA / f"{PREFIX}-{c}-observations.csv").read_bytes()
            ).hexdigest()
            for c in SOURCE
        },
        "original_author_input_files_recovered": False,
        "author_fit_parity": False,
        "author_contact": "none",
        "qualification": FACTS["qualification"],
        "variants": {},
    }
    for variant, source in VARIANTS.items():
        compound = source["compound"]
        report["variants"][variant] = {"source": source, "cases": {}}
        for name, options in CASES.items():
            folder = output / variant / name
            folder.mkdir(parents=True)
            data, mapping = prepare_input(compound, **options)
            (folder / "input.dat").write_text(data, encoding="utf-8")
            (folder / "row-map.json").write_text(
                json.dumps(mapping, indent=2) + "\n", encoding="utf-8"
            )
            (folder / "run.mcr").write_text(
                commands(variant, folder, options), encoding="utf-8"
            )
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
                raise RuntimeError(f"Native Case 3 run failed: {folder}")
            stdout = process.stdout.decode("utf-8", errors="replace")
            case = parse_console(
                stdout,
                "rubidium",
                n_observations=len(mapping),
                parameter_names=PARAMETERS,
            )
            if not case["converged"]:
                raise RuntimeError(f"No convergence: {folder}")
            match = re.search(r"Version:\s*([\d.]+)", stdout)
            report["version"] = match.group(1) if match else "unknown"
            case["input_options"] = options
            case["observations"] = len(mapping)
            case["convergence_messages"] = stdout.count("Least Squares converged")
            case["fixed_parameters"] = list(source["fixed_parameters"])
            if options.get("rounding_anchor_example"):
                case["fixed_parameters"][:3] = ROUNDING_ANCHOR_EXAMPLES[variant]
                case["fixed_anchor_qualification"] = (
                    "Illustrative fitted-DAC anchors within published half-digit intervals, "
                    "chosen to test rounding compatibility; not recovered author values. "
                    "Deliberately selected RbBr Kprime=4.98 remains exact."
                )
            if options.get("precise_theta"):
                case["fixed_parameters"][3] = AUTHOR_THETA[compound]
            case["published_comparison"] = comparison(
                variant, case["refined_parameters_final_cycle"]
            )
            case["saved_thermal_covariance"] = covariance(
                folder / "fitted.eos", indices=(17, 18)
            ).tolist()
            case["covariance_parameter_order"] = PARAMETERS
            case["covariance_qualification"] = (
                "Native chi-squared-scaled covariance conditional on all four fixed anchors; their source uncertainties are not propagated."
            )
            case["files_sha256"] = {
                f.name: hashlib.sha256(f.read_bytes()).hexdigest()
                for f in sorted(folder.iterdir())
            }
            report["variants"][variant]["cases"][name] = case
            print(
                variant,
                name,
                case["published_comparison"]["eosfit_parameters"],
                flush=True,
            )
            (output / "manifest.json").write_text(
                json.dumps(report, indent=2) + "\n", encoding="utf-8"
            )
    report["rounding_examples_all_match_publication"] = all(
        v["cases"]["weighted_anchor_rounding_example"]["published_comparison"][
            "both_parameters_round_to_publication"
        ]
        for v in report["variants"].values()
    )
    report["runtime_reference"] = f"docs/data/{PREFIX}-eosfit7c/manifest.json"
    (output / "manifest.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    run(args.executable.resolve(), args.output.resolve())


if __name__ == "__main__":
    main()
