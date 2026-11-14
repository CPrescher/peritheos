"""Inspect Ross Angel's supplied processed pyrope inputs, without replacing data.

Repository-only: python -m scripts.audit_pyrope_author_inputs --source-dir DIR.
Original attachments stay in the private archive; exported CSVs are attributed
factual transcriptions. This does not recover the preceding processing history.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from collections import Counter
from pathlib import Path

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_pyrope_angel_2022_observations import DATA
from scripts.refit_pyrope_angel_2022_expanded import expanded_rows
from scripts.refit_pyrope_angel_2022_subset import (
    arrays,
    effective_sigma,
    pressure_at_observations,
    solve,
)

PREFIX = "author-2026-11-14"
NAMES = ["K0_gpa", "K0_prime", "theta0_k", "gamma0", "scale_sinogeikin_2000_ks"]
ATTACHMENTS = [
    "pyrope_qcompMGD_March2023.eos",
    "pyrope datasets.pdf",
    "Pyrope_PKT_used.dat",
    "pyrope_25Oct2021.log",
    "Pyrope_PVT_used.dat",
]
GROUPS = [
    ("Zhang et al (1998) Helium", "zhang_helium"),
    ("Zhang et al (1998) Neon", "zhang_neon"),
    ("Zou et al (2012) J Appl", "zou_pv"),
    ("Du et al (2015)", "du_pv"),
    ("Milani et al (2015) single crystal", "milani_pv"),
    ("Chantel et al (2016): room T", "chantel_pv"),
    ("Milani et al 2015, with esd(T)", "milani_vt"),
    ("Hartwig and Galkin", "hartwig_vt"),
    ("Sinogeikin and Bass (2000)", "sinogeikin_2000_ks"),
    ("Sinogeikin and Bass (2002)", "sinogeikin_2002_ks"),
    ("Zou et al 2012b", "zou_ks"),
    ("Chantel et al 2016 (highT", "chantel_ks"),
]


def read_dat(path, kind):
    result, dataset, index = [], None, 0
    for line_number, line in enumerate(path.read_text().splitlines(), 1):
        if line.startswith("#"):
            matches = [name for prefix, name in GROUPS if line[1:].startswith(prefix)]
            if matches:
                dataset, index = matches[0], 0
        if not re.match(r"^\s*\d", line):
            continue
        if dataset is None or len(line.split()) != 6:
            raise ValueError(f"Unexpected observation at {path.name}:{line_number}")
        index += 1
        tokens = dict(zip(["T", "sT", "P", "sP", "Y", "sY"], line.split()))
        result.append(
            dict(
                source_file=path.name,
                source_line=line_number,
                dataset=dataset,
                row=index,
                kind=kind,
                **{key: float(value) for key, value in tokens.items()},
                **{key + "_token": value for key, value in tokens.items()},
            )
        )
    return result


def read_log(path):
    text = path.read_text()
    result = []
    for line in text.splitlines():
        if not re.match(r"^\s*\d+\s+[01]\s+\d+\s+(V|KS)\s", line):
            continue
        tokens = line.split()
        if len(tokens) != 17:
            raise ValueError("Unexpected author log row")
        result.append(tokens)
    if len(result) != 194:
        raise ValueError("Expected 194 original log observations")
    return text, result


def author_rows(directory):
    rows = read_dat(directory / "Pyrope_PVT_used.dat", "V") + read_dat(
        directory / "Pyrope_PKT_used.dat", "KS"
    )
    _, log = read_log(directory / "pyrope_25Oct2021.log")
    if len(rows) != len(log):
        raise ValueError("DAT/log row counts differ")
    checks = []
    for i, (row, printed) in enumerate(zip(rows, log), 1):
        if int(printed[0]) != i or printed[3] != row["kind"]:
            raise ValueError("DAT/log order differs")
        row.update(fit_row=i, use=int(printed[1]), group=int(printed[2]))
        row["scale_index"] = 4 if row["group"] == 3 else None
        # The log rounds tiny pressure errors to zero; compare printed precision,
        # retaining the full DAT tokens as the actual coordinate/error inputs.
        for key, value in zip(["P", "sP", "T", "sT", "Y", "sY"], printed[4:10]):
            half = 0.5 * 10 ** -len(value.split(".")[1])
            if abs(row[key] - float(value)) > half + 2e-6:
                raise ValueError(f"DAT/log mismatch: row {i}, {key}")
        checks.append(
            dict(
                fit_row=i,
                calculated_pressure_gpa_token=printed[10],
                calculated_value_token=printed[13],
                weight_token=printed[16],
            )
        )
    return rows, checks


def prepared(rows):
    return [
        {
            k: row[k]
            for k in [
                "dataset",
                "row",
                "kind",
                "T",
                "P",
                "Y",
                "sT",
                "sP",
                "sY",
                "scale_index",
            ]
        }
        for row in rows
        if row["use"]
    ]


def comparison(rows):
    old = expanded_rows(include_temperature=True)
    result = []
    for dataset in sorted({r["dataset"] for r in rows}):
        before = [r for r in old if r["dataset"] == dataset]
        after = [r for r in rows if r["dataset"] == dataset]
        item = dict(
            dataset=dataset, reconstructed_rows=len(before), author_rows=len(after)
        )
        # Compare row sequences only when selection/order is shared; otherwise
        # expose the full arrays rather than asserting a false row pairing.
        item["author_ranges"] = {
            key: [min(r[key] for r in after), max(r[key] for r in after)]
            for key in ["T", "sT", "P", "sP", "Y", "sY"]
        }
        item["reconstruction_ranges"] = {
            key: [min(r[key] or 0 for r in before), max(r[key] or 0 for r in before)]
            for key in ["T", "sT", "P", "sP", "Y", "sY"]
        }
        if len(before) == len(after) and dataset != "sinogeikin_2000_ks":
            item["sequence_maximum_absolute_differences"] = {
                key: max(abs((a[key] or 0) - b[key]) for a, b in zip(before, after))
                for key in ["T", "sT", "P", "sP", "Y", "sY"]
            }
            item["author_over_reconstructed_value_range"] = [
                min(b["Y"] / a["Y"] for a, b in zip(before, after)),
                max(b["Y"] / a["Y"] for a, b in zip(before, after)),
            ]
        result.append(item)
    return result


def write_csv(name, rows):
    with (DATA / name).open("w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def check_covariance(rows, fit):
    """Independent five-point derivatives and QR inverse-information check.

    Coordinate weights are held at the fitted state, matching the stated
    effective-variance least-squares convention, not a likelihood Hessian
    including parameter-dependent error-normalization terms.
    """
    c = np.array(list(fit["parameters"].values()))
    d = arrays(rows)
    sigma = effective_sigma(c, d)
    scale = np.array([170, 4.5, 770, 1.2, 1])
    covariance = np.array(fit["conditional_unscaled_covariance"])
    steps = []
    for fraction in [3e-5, 1e-4, 3e-4]:
        columns = []
        for i in range(5):
            h = scale[i] * fraction
            values = []
            for multiple in [-2, -1, 1, 2]:
                trial = c.copy()
                trial[i] += multiple * h
                values.append(pressure_at_observations(trial, d))
            columns.append(
                (values[0] - 8 * values[1] + 8 * values[2] - values[3])
                / (12 * h * sigma)
            )
        jac = np.array(columns).T * scale
        _, upper = np.linalg.qr(jac, mode="reduced")
        inverse = np.linalg.solve(upper, np.eye(5))
        checked = (inverse @ inverse.T) * np.outer(scale, scale)
        steps.append(
            dict(
                parameter_step_fraction=fraction,
                maximum_width_relative_difference=float(
                    np.max(np.abs(np.sqrt(np.diag(checked) / np.diag(covariance)) - 1))
                ),
                scaled_inverse_information_maximum_error=float(
                    np.max(
                        np.abs(
                            (jac.T @ jac) @ (checked / np.outer(scale, scale))
                            - np.eye(5)
                        )
                    )
                ),
            )
        )
    chi = fit["conditional_reduced_chi_square"]
    widths = np.sqrt(np.diag(covariance))
    scaled = widths * np.sqrt(chi)
    np.testing.assert_allclose(
        scaled, list(fit["conditional_residual_scaled_widths"].values())
    )
    residual = pressure_at_observations(c, d) - d["P"]
    np.testing.assert_allclose(np.sum((residual / sigma) ** 2) / (194 - 5), chi)
    return dict(
        method="Independent five-point pressure Jacobians at three step sizes; QR inverse information versus production three-point Jacobians and SVD",
        finite_difference_checks=steps,
        covariance_positive_definite=bool(np.linalg.eigvalsh(covariance).min() > 0),
        unscaled_widths=dict(zip(NAMES, widths.tolist())),
        residual_scaled_widths=dict(zip(NAMES, scaled.tolist())),
        variance_multiplier=chi,
        standard_error_multiplier=float(np.sqrt(chi)),
        degrees_of_freedom=189,
        correlation=(covariance / np.outer(widths, widths)).tolist(),
        qualification="Conditional local one-standard-error widths with residual rescaling; not formal total experimental uncertainty, bootstrap intervals or independent-two-sample significance tests",
    )


def run(directory):
    rows, log = author_rows(directory)
    p, covariance = source("pyrope")
    eos = (directory / ATTACHMENTS[0]).read_text()
    supplied = {
        int(k): float(v) for k, v in re.findall(r"^Param =\s*(\d+)\s+(\S+)", eos, re.M)
    }
    matrix = np.array(
        [
            [float(x) for x in re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)]
            for line in eos.split("Variance-Covariance matrix=")[1].splitlines()
            if len(re.findall(r"[-+]?\d\.\d+E[-+]\d+", line)) == 59
        ]
    )
    np.testing.assert_array_equal(matrix, covariance)
    if supplied != p:
        raise ValueError("Author EOS differs numerically from archived publication")
    data = prepared(rows)
    fits = []
    for initial in [[p[2], p[3], p[11], p[18], p[53]], [170, 4.5, 800, 1.2, 0.98]]:
        fits.append(
            solve(
                data,
                initial,
                list(range(5)),
                names=NAMES,
                scaling=np.array([170, 4.5, 770, 1.2, 1]),
                lower=np.array([50, 2, 50, 0.05, 0.9]),
                upper=np.array([350, 10, 3000, 4, 1.1]),
            )
        )
    write_csv(PREFIX + "-inputs.csv", rows)
    write_csv(PREFIX + "-log-calculations.csv", log)
    # EOS and the later final log are distinct parameter/covariance snapshots.
    text, _ = read_log(directory / "pyrope_25Oct2021.log")
    block = text.split("RESULTS AFTER FINAL REFINEMENT CYCLE")[1].split("INPUT <CR>")[0]
    reported = {}
    for label, key in zip(["K0", "Kp", "ThMGD", "Gamm0", "Sca 3"], NAMES):
        match = re.search(
            r"^\s*" + label + r"\s+1\s+(\S+)\s+(\S+)\s+(\S+)", block, re.M
        )
        reported[key] = dict(value=float(match[1]), rescaled_esd=float(match[3]))
    report = dict(
        received_on="2026-11-14",
        scope="Author-supplied processed-input audit and independent fit replay; upstream preprocessing history not fully recovered",
        attachments=[
            dict(
                file=name,
                sha256=hashlib.sha256((directory / name).read_bytes()).hexdigest(),
                bytes=(directory / name).stat().st_size,
            )
            for name in ATTACHMENTS
        ],
        source_directory=str(directory),
        n_rows=len(rows),
        counts=dict(Counter(r["dataset"] for r in rows)),
        all_rows_used=all(r["use"] == 1 for r in rows),
        group_counts=dict(Counter(str(r["group"]) for r in rows)),
        author_eos_parameters_and_full_covariance_equal_published=True,
        author_eos_byte_identical=False,
        author_eos_difference="Title includes citation; scientific parameter and full 59x59 covariance blocks are identical",
        final_log_parameters=reported,
        final_log_reduced_chi_square=2.2335,
        reconstructed_input_comparison=comparison(rows),
        independent_covariance_check=check_covariance(data, fits[0]),
        independent_fits=fits,
        original_objective="Pressure residuals, weights propagated from all supplied P/T/value coordinate errors, scales multiply calculated values; V0 fixed; K0, Kprime, theta, gamma and group3 scale free",
        original_scale_assignment="Group1: all115volumes; group2:64unscaled moduli including two Sinogeikin2000ambient points; group3:15 Sinogeikin2000pressure points",
        covariance_semantics="Author log explicitly reports residual-rescaled covariance and esds; supplied EOS is a nearby earlier snapshot, not the final log's full-precision coefficient vector",
        rights="Attachments remain private and are not distributed; factual transcriptions and comparison metadata retain author attribution, without inferring a third-party license",
        author_contacted=False,
        catalog_changed=False,
    )
    (DATA / (PREFIX + "-audit.json")).write_text(json.dumps(report, indent=2) + "\n")
    print(
        json.dumps(
            {
                "counts": report["counts"],
                "fits": [
                    {
                        k: f[k]
                        for k in [
                            "parameters",
                            "conditional_reduced_chi_square",
                            "numerically_converged",
                        ]
                    }
                    for f in fits
                ],
            },
            indent=2,
        )
    )
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    run(parser.parse_args().source_dir.resolve())
