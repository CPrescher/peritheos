"""Replay/refit author-supplied FeS files with an external official EosFit.

The original input bytes, including repeated observations and errors, are
retained. Diagnostics do not replace the database source coefficients.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
import shutil
import subprocess
from pathlib import Path

import numpy as np

from peritheos import __version__, _rust
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from scripts.audit_morard_2026_fes import CELL_PER_MOLAR, observations
from scripts.audit_sata_2010_fes_vi import observations as cold_observations
from scripts.run_fes_eosfit_console import (
    enrich_saved_model,
    parse_console,
    saved_covariance,
)

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/fes-author-eosfit"


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_data(path):
    """Honor the actual Format header: P, sigmaP, T, sigmaT, V, sigmaV."""
    header = path.read_text(encoding="utf-8").splitlines()[0].split()
    if [word.lower() for word in header] != [
        "format",
        "1",
        "p",
        "sigp",
        "t",
        "sigt",
        "v",
        "sigv",
    ]:
        raise ValueError(f"Unexpected author input columns: {header}")
    data = np.loadtxt(path, skiprows=1)
    if data.ndim != 2 or data.shape[1] != 6 or not np.all(np.isfinite(data)):
        raise ValueError("Invalid author observations")
    return data


def summary(residual, indices):
    return {
        "rows": len(residual),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "min_gpa": float(residual.min()),
        "max_gpa": float(residual.max()),
        "max_abs_gpa": float(np.max(np.abs(residual))),
        "rows_exceeding_3_gpa": (indices[np.abs(residual) > 3] + 1).tolist(),
    }


def groups(data, residual):
    hot = data[:, 2] > 300
    masks = {
        "all": np.ones(len(data), dtype=bool),
        "thermal": hot,
        "cold": ~hot,
        "thermal_pressure_ge_40": hot & (data[:, 0] >= 40),
    }
    return {
        key: summary(residual[mask], np.flatnonzero(mask))
        for key, mask in masks.items()
    }


def initial_rows(stdout, n):
    rows = []
    for line in stdout.splitlines():
        values = line.split()
        if len(values) == 14 and values[0].isdigit() and values[1] == "1":
            rows.append({"row": int(values[0]), "residual_gpa": -float(values[9])})
    rows = rows[-n:]
    if [r["row"] for r in rows] != list(range(1, n + 1)):
        raise ValueError("Incomplete EosFit replay table")
    return rows


def run_macro(executable, folder, commands):
    folder.mkdir()
    path = folder / "run.mcr"
    path.write_text(
        "\n".join([f"log {folder / 'run.log'}", *commands]) + "\n", encoding="utf-8"
    )
    result = subprocess.run(
        [str(executable)],
        input=f"macro {path}\nexit\n",
        text=True,
        capture_output=True,
        cwd=folder,
        timeout=120,
    )
    (folder / "stdout.txt").write_text(result.stdout, encoding="utf-8")
    (folder / "stderr.txt").write_text(result.stderr, encoding="utf-8")
    if result.returncode:
        raise RuntimeError(f"EosFit failed: {folder}")
    return result.stdout


def compare_table(data):
    d = observations()
    reference = np.column_stack(
        [
            d["pressure_gpa"],
            d["sigma_pressure_gpa"],
            d["temperature_k"],
            d["sigma_temperature_k"],
            d["volume_angstrom3"] / CELL_PER_MOLAR,
            d["sigma_volume_angstrom3"] / CELL_PER_MOLAR,
        ]
    )
    hot = data[data[:, 2] > 300]
    delta = hot - reference
    # Both sources contain a repeated seventh thermal observation. Match by
    # sequence, rather than a nearest-neighbour match that collapses duplicates.
    return {
        "thermal_rows": len(hot),
        "same_thermal_pressure_temperature_order": bool(
            len(hot) == len(reference)
            and np.allclose(hot[:, :4], reference[:, :4], rtol=0, atol=5e-7)
        ),
        "max_absolute_column_differences": np.max(np.abs(delta), axis=0).tolist(),
        "volume_conversion_cell_to_cm3_mol": 0.15055,
        "max_volume_difference_from_author_rounded_conversion_cm3_mol": float(
            np.max(np.abs(hot[:, 4] - d["volume_angstrom3"] * 0.15055))
        ),
        "qualification": "Author volumes use the rounded factor 0.15055. All rows and duplicated observations are preserved.",
    }


def compare_cold_sources(data):
    literature = data[(data[:, 2] == 300) & (data[:, 1] != 0.0056)]
    matches = []
    matched = set()
    for source in cold_observations():
        if not source["included"]:
            continue
        volume = source["volume_angstrom3"] * 0.15055
        candidates = np.flatnonzero(np.abs(literature[:, 4] - volume) < 1e-8)
        if len(candidates) != 1:
            raise ValueError("Expected one match for every re-reported Sata VI row")
        i = int(candidates[0])
        matched.add(i)
        matches.append(
            {
                "sata_table_fes_row": source["table_fes_row"],
                "author_literature_row": i + 1,
                "author_pressure_gpa": float(literature[i, 0]),
                "sata_pressure_gpa": source["pressure_gpa"],
                "author_pressure_error_gpa": float(literature[i, 1]),
                "sata_pressure_error_gpa": source["error_pressure_gpa"],
                "author_volume_error_cm3_mol": float(literature[i, 5]),
                "sata_volume_error_cm3_mol": source["error_volume_angstrom3"] * 0.15055,
            }
        )
    return {
        "literature_cold_rows": len(literature),
        "sata_vi_volume_matches": matches,
        "additional_cold_rows": [
            x.tolist() for i, x in enumerate(literature) if i not in matched
        ],
        "qualification": "13 rows match the re-reported Sata VI cell volumes. Eight further cold rows are supplied in the author's combined Sata/Ohfuji input; their individual source identity is not encoded. One matched pressure is 101.2 GPa here versus 101.1 GPa in Sata Table 1. Author input uncertainties are preserved as assigned, not treated as the original table errors.",
    }


def cold_refits(executable, output, saved, name, data):
    """Fit the 300 K subset with BM3, retaining each weighting choice."""
    data = data[data[:, 2] == 300]
    input_file = output / "sources" / f"{name}-derived-cold.dat"
    np.savetxt(
        input_file,
        data,
        header="Format 1 P sigP T sigT V sigV",
        comments="",
        fmt="%.17g",
    )
    results = {}
    schemes = [("unit", ["n"] * 3), ("errors", ["y"] * 3)]
    if name == "without-own-cold":
        schemes += [("volume", ["n", "n", "y"]), ("pressure", ["y", "n", "n"])]
    for scheme, bits in schemes:
        folder = output / f"{name}-cold-{scheme}"
        commands = [
            "input",
            f"load {saved}",
            "th",
            "0",
            "x",
            f"read {input_file}",
            "fit",
            "n",
            "y",
            "y",
            "y",
            *bits,
            "",
            "y",
            "n",
            f"save {folder / 'fitted.eos'}",
            "y",
            "exit",
        ]
        stdout = run_macro(executable, folder, commands)
        case = parse_console(
            stdout,
            "combined",
            n_observations=len(data),
            parameter_names=["V0", "K0", "Kp"],
        )
        if not case["converged"]:
            raise RuntimeError(f"Refinement did not converge: {folder}")
        case["covariance_parameter_order"] = ["V0", "K0", "Kp"]
        case["saved_refined_parameter_covariance"] = saved_covariance(
            folder / "fitted.eos", "cold"
        ).tolist()
        case["cold_thermal_model"] = "No thermal term: BM3 cold fit at 300 K"
        case["weights"] = {
            "unit": "equal pressure weights",
            "errors": "supplied P/V/T uncertainties",
            "volume": "supplied volume uncertainty only",
            "pressure": "supplied pressure uncertainty only",
        }[scheme]
        results[folder.name] = case
    return results


def replay(executable, output, eos_path, input_paths):
    output.mkdir(parents=True, exist_ok=False)
    sources = output / "sources"
    sources.mkdir()
    original_paths = [eos_path, *input_paths]
    report = {
        "executable_sha256": digest(executable),
        "executable_source": "https://www.rossangel.com/Download/EoSFit7installer_MacOS.dmg",
        "peritheos_version": __version__,
        "native_extension_sha256": digest(Path(_rust.__file__)),
        "source_files": {},
        "qualification": "Author-supplied input/model bytes, replayed unchanged. Refits are diagnostic; the final author fit mask, free-parameter selection and weighting history are not encoded in these files. Reported esds and covariance are conditional EosFit results, not calibration or physical validation.",
        "cases": {},
    }
    for source in original_paths:
        target = sources / source.name
        shutil.copyfile(source, target)
        report["source_files"][source.name] = {
            "original_path": str(source),
            "sha256": digest(target),
            "bytes": target.stat().st_size,
        }
    saved = sources / eos_path.name
    text = saved.read_text(encoding="utf-8")
    parameters = {
        int(i): float(v) for i, v in re.findall(r"Param\s*=\s*(\d+)\s+(\S+)", text)
    }
    for field, expected in [("Model", 2), ("Order", 3), ("Thermal", 7)]:
        match = re.search(rf"^{field}\s*=\s*(\d+)", text, re.M)
        if match is None or int(match.group(1)) != expected:
            raise ValueError(f"Author replay requires BM3/full MGD: {field}")
    if parameters[14] != 0:
        raise ValueError("Author model is not full MGD")
    tref = float(re.search(r"^Tref\s*=\s*(\S+)", text, re.M).group(1))
    covariance_text = text.split("Variance-Covariance matrix=")[1]
    covariance = [
        [float(x) for x in re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)]
        for line in covariance_text.splitlines()
        if len(re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)) == 59
    ][:59]
    if np.shape(covariance) != (59, 59):
        raise ValueError("Author covariance matrix missing")
    report["saved_author_model"] = {
        "V0_cm3_mol": parameters[1],
        "K0_gpa": parameters[2],
        "Kprime": parameters[3],
        "theta0_k": parameters[11],
        "n": parameters[13],
        "q_compromise_flag": parameters[14],
        "gamma0": parameters[18],
        "q": parameters[19],
        "Tref_k": tref,
        "saved_date": "2024-03-12",
        "gui_version": "20210609",
        "covariance_nonzero_entries": [
            {
                "parameter_i": int(i + 1),
                "parameter_j": int(j + 1),
                "value": covariance[i][j],
            }
            for i, j in zip(*np.nonzero(covariance))
        ],
        "gamma0_saved_esd": float(np.sqrt(covariance[17][17])),
        "covariance_qualification": "Only the gamma0 variance is nonzero. This is consistent with a conditional gamma-only refinement; cold-parameter uncertainty is not supplied by this saved matrix.",
    }
    eos = MieGruneisenDebye(
        BM3(parameters[1] / 10, parameters[2], parameters[3]),
        tref,
        parameters[11],
        parameters[18],
        parameters[19],
        int(parameters[13]),
    )
    if not hasattr(eos, "_native"):
        raise RuntimeError("This comparison requires native Peritheos")
    data_sets = {}
    for original in input_paths:
        input_file = sources / original.name
        name = "without-own-cold" if "-300K" in original.stem else "all-author-rows"
        data = read_data(input_file)
        data_sets[name] = data
        folder = output / f"{name}-replay"
        commands = ["input", f"load {saved}", "x", f"read {input_file}", "list", "exit"]
        stdout = run_macro(executable, folder, commands)
        rows = initial_rows(stdout, len(data))
        residual = np.array([row["residual_gpa"] for row in rows])
        native = np.asarray(eos.pressure(data[:, 4] / 10, data[:, 2])) - data[:, 0]
        report["cases"][folder.name] = {
            "eosfit": groups(data, residual),
            "peritheos": groups(data, native),
            "max_pressure_difference_gpa": float(np.max(np.abs(native - residual))),
            "rows": rows,
            "thermal_table_comparison": compare_table(data),
            "cold_source_comparison": compare_cold_sources(data),
        }
        schemes = [
            ("gamma", "errors", ["y"] * 3),
            ("joint", "unit", ["n"] * 3),
            ("joint", "errors", ["y"] * 3),
        ]
        if name == "without-own-cold":
            schemes += [
                ("gamma", "volume", ["n", "n", "y"]),
                ("gamma", "pressure", ["y", "n", "n"]),
            ]
        for mode, scheme, bits in schemes:
            folder = output / f"{name}-{mode}-{scheme}"
            commands = ["input", f"load {saved}", "x", f"read {input_file}", "fit", "n"]
            commands += ["y" if mode == "joint" else "n"] * 3
            commands += ["n", "y", "n"]  # theta fixed, gamma free, q fixed
            commands += bits
            commands += ["", "y", "n", f"save {folder / 'fitted.eos'}", "y", "exit"]
            stdout = run_macro(executable, folder, commands)
            dataset = "combined" if mode == "joint" else "morard"
            case = parse_console(stdout, dataset, n_observations=len(data))
            if not case["converged"]:
                raise RuntimeError(f"Refinement did not converge: {folder}")
            enrich_saved_model(case, folder, dataset)
            case["weights"] = {
                "errors": "supplied P/V/T uncertainties",
                "unit": "equal pressure weights",
                "volume": "supplied volume uncertainty only",
                "pressure": "supplied pressure uncertainty only",
            }[scheme]
            residual = np.array(
                [r["residual_gpa"] for r in case["fitted_replay"]["rows"]]
            )
            case["groups"] = groups(data, residual)
            report["cases"][folder.name] = case
        report["cases"].update(cold_refits(executable, output, saved, name, data))
    a, b = data_sets["all-author-rows"], data_sets["without-own-cold"]
    report["difference_between_inputs"] = {
        "removed_from_300K_named_file": [
            x.tolist() for x in a if not any(np.array_equal(x, y) for y in b)
        ],
        "extra_in_300K_named_file": [
            x.tolist() for x in b if not any(np.array_equal(x, y) for y in a)
        ],
        "columns": [
            "P_GPa",
            "sigmaP_GPa",
            "T_K",
            "sigmaT_K",
            "V_cm3_mol",
            "sigmaV_cm3_mol",
        ],
    }
    report["artifact_sha256"] = {
        path.relative_to(output).as_posix(): digest(path)
        for path in sorted(output.rglob("*"))
        if path.is_file()
    }
    (output / "manifest.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--executable", type=Path, required=True)
    parser.add_argument("--eos", type=Path, required=True)
    parser.add_argument("--inputs", type=Path, nargs=2, required=True)
    parser.add_argument("--output", type=Path, default=OUTPUT)
    args = parser.parse_args()
    report = replay(
        args.executable.resolve(),
        args.output.resolve(),
        args.eos.resolve(),
        [p.resolve() for p in args.inputs],
    )
    for name, case in report["cases"].items():
        print(name, case.get("groups", case.get("eosfit")))
        if "refined_parameters_final_cycle" in case:
            print(case["refined_parameters_final_cycle"])


if __name__ == "__main__":
    main()
