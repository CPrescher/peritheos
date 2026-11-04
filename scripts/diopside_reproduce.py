"""Audit the published diopside equations and uncorrected source tables.

Run with ``python -m scripts.diopside_reproduce``. Independent equations and
explicit diagnostic objectives distinguish equation agreement, coefficient
agreement, and the unavailable original author regression workflows.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import re
from decimal import Decimal
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/diopside-reproduction.json"
THERMAL = np.array([109.1, 4.84, 2.32e-5, 1.88e-8, -0.0205])
STATIC = np.array([438.66, 118.0, 3.8])


def split_error(text):
    """Expand the source's parenthesized last digits with decimal arithmetic."""
    match = re.fullmatch(r"([\d.]+)(?:\((\d+)\))?", text)
    if match is None:
        raise ValueError(text)
    value, error = match.groups()
    return [
        value,
        str(Decimal(error) * Decimal(10) ** Decimal(value).as_tuple().exponent)
        if error
        else "",
    ]


def normalized_csv(stem):
    fields = [
        "lattice_a_angstrom",
        "lattice_b_angstrom",
        "lattice_c_angstrom",
        "beta_deg",
        "cell_volume_a3",
    ]
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    lines = (DATA / (stem + "-original.txt")).read_text(encoding="utf-8").splitlines()
    for index, line in enumerate(line for line in lines if not line.startswith("#")):
        tokens = line.split()
        if "zhao" in stem:
            row = [str(index + 1), tokens[0], "", tokens[1], ""]
            columns = [
                "source_row",
                "pressure_gpa",
                "pressure_error_gpa",
                "temperature_k",
                "temperature_error_k",
            ]
            values = tokens[2:]
        else:
            row = [str(index + 1), tokens[0], *split_error(tokens[1]), tokens[2]]
            columns = [
                "source_row",
                "run",
                "pressure_gpa",
                "pressure_error_gpa",
                "intensity_collected",
            ]
            values = tokens[3:]
        if len(values) != len(fields):
            raise ValueError(f"Wrong source column count on row {index + 1}")
        for name, value in zip(fields, values):
            row.extend(split_error(value))
            columns.extend(
                [
                    name,
                    name.replace("_angstrom", "_error_angstrom")
                    .replace("_deg", "_error_deg")
                    .replace("_a3", "_error_a3"),
                ]
            )
        if index == 0:
            writer.writerow(columns)
        writer.writerow(row)
    return stream.getvalue()


def bm3(volume, v0, k0, kp):
    """Independent finite-strain equation, without importing Peritheos."""
    f = ((v0 / np.asarray(volume)) ** (2.0 / 3.0) - 1.0) / 2.0
    return 3.0 * k0 * f * (1.0 + 2.0 * f) ** 2.5 * (1.0 + 1.5 * (kp - 4.0) * f)


def thermal_pressure(volume, temperature, coefficients=THERMAL):
    k0, kp, a, b, dkdt = coefficients
    temperature = np.asarray(temperature)
    delta = temperature - 300.0
    v0_t = 438.67 * np.exp(a * delta + 0.5 * b * (temperature**2 - 300.0**2))
    return bm3(volume, v0_t, k0 + dkdt * delta, kp)


def load_rows(stem):
    with (DATA / (stem + ".csv")).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def values(rows, column):
    return np.array([float(row[column]) for row in rows])


def statistics(residual):
    return {
        "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "maximum_absolute_pressure_residual_gpa": float(np.max(np.abs(residual))),
    }


def metric_differences(rows):
    return values(rows, "lattice_a_angstrom") * values(
        rows, "lattice_b_angstrom"
    ) * values(rows, "lattice_c_angstrom") * np.sin(
        np.deg2rad(values(rows, "beta_deg"))
    ) - values(rows, "cell_volume_a3")


def reproduce():
    stems = ["diopside-zhao-1998-table1", "diopside-thompson-2008-table1"]
    for stem in stems:
        if normalized_csv(stem).encode() != (DATA / (stem + ".csv")).read_bytes():
            raise ValueError(f"CSV differs from retained source notation: {stem}")
    zhao, thompson = (load_rows(stem) for stem in stems)
    v, t, p = (
        values(zhao, c) for c in ("cell_volume_a3", "temperature_k", "pressure_gpa")
    )
    residual = thermal_pressure(v, t) - p
    thermal_fits = []
    for excluded in ([], [15]):
        mask = np.array([int(r["source_row"]) not in excluded for r in zhao])
        fit = least_squares(
            lambda x: thermal_pressure(v[mask], t[mask], x) - p[mask],
            THERMAL,
            x_scale=np.abs(THERMAL),
            ftol=1e-12,
            xtol=1e-12,
            gtol=1e-12,
        )
        thermal_fits.append(
            {
                "excluded_source_rows": excluded,
                "observations": int(mask.sum()),
                "coefficients": dict(
                    zip(
                        ["K0", "K0_prime", "alpha0", "alpha1", "dK_dT"],
                        map(float, fit.x),
                    )
                ),
                "fixed_V0_a3": 438.67,
                "success": bool(fit.success),
                **statistics(fit.fun),
            }
        )
    sv, sp = (values(thompson, c) for c in ("cell_volume_a3", "pressure_gpa"))
    sigma_v = values(thompson, "cell_volume_error_a3")
    # Ambient pressure is an imposed reference in this diagnostic, rather than
    # an invented measurement with zero uncertainty. Other P errors are printed.
    sigma_p = np.array([float(r["pressure_error_gpa"] or 0.0) for r in thompson])

    def static_residual(x):
        k = -sv * (bm3(sv + 0.0001, *x) - bm3(sv - 0.0001, *x)) / 0.0002
        return (bm3(sv, *x) - sp) / np.sqrt(sigma_p**2 + (k * sigma_v / sv) ** 2)

    static_fit = least_squares(
        static_residual, STATIC, x_scale=STATIC, ftol=1e-12, xtol=1e-12, gtol=1e-12
    )
    return {
        "audit_date": "2026-10-15",
        "scope": "Independent equation validation and diagnostic refits; no executable refit records and no original-author workflow parity claim.",
        "resources": {
            path.name: hashlib.sha256(path.read_bytes()).hexdigest()
            for path in sorted(DATA.glob("diopside*"))
            if path.is_file()
        },
        "zhao_1998": {
            "observations": len(zhao),
            "duplicate_ambient_source_rows": [1, 32],
            "published_curve_all_rows": statistics(residual),
            "published_curve_without_row15_diagnostic_only": statistics(
                np.delete(residual, 14)
            ),
            "row15_printed_pressure_gpa": 6.95,
            "row15_predicted_pressure_gpa": float(thermal_pressure(431.66, 1072)),
            "source_mean_alpha_300_1300_per_k": 3.82e-5,
            "integrated_mean_alpha_300_1300_per_k": float(
                THERMAL[2] + THERMAL[3] * 800
            ),
            "metric_minus_printed_volume_a3": list(
                map(float, metric_differences(zhao))
            ),
            "diagnostic_objective": "Unweighted pressure residuals; V0 fixed to printed ambient volume, Tr=300 K. Source weights, V0 fit treatment, selection details and covariance are unreported. Row15 exclusion is a sensitivity experiment, never a source correction.",
            "diagnostic_refits": thermal_fits,
            "author_fit_reproduced": False,
            "author_fit_reproduced_scope": "Exact original numerical fit input and workflow only. Published-model and parameter parity after accounting for the documented table error are assessed separately in diopside-plot-check.json.",
            "author_workflow_parity": False,
        },
        "thompson_2008": {
            "observations": len(thompson),
            "published_curve_all_rows": statistics(bm3(sv, *STATIC) - sp),
            "metric_minus_printed_volume_a3": list(
                map(float, metric_differences(thompson))
            ),
            "diagnostic_objective": "Effective pressure variance sigma_P^2+(K*sigma_V/V)^2; ambient P treated as imposed reference. Parameter-dependent weights; source says weighted Levenberg-Marquardt but does not specify this objective.",
            "diagnostic_refit": {
                "coefficients": dict(
                    zip(["V0", "K0", "K0_prime"], map(float, static_fit.x))
                ),
                "success": bool(static_fit.success),
                **statistics(bm3(sv, *static_fit.x) - sp),
            },
            "author_workflow_parity": False,
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=REPORT)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    if args.check:
        stored = json.loads(args.output.read_text(encoding="utf-8"))
        # Optimizer last digits can vary across supported SciPy/platform builds.
        for source in ("zhao_1998", "thompson_2008"):
            assert result[source]["observations"] == stored[source]["observations"]
            for key, value in result[source]["published_curve_all_rows"].items():
                np.testing.assert_allclose(
                    value, stored[source]["published_curve_all_rows"][key], atol=1e-12
                )
        for current, previous in zip(
            result["zhao_1998"]["diagnostic_refits"],
            stored["zhao_1998"]["diagnostic_refits"],
        ):
            assert current["excluded_source_rows"] == previous["excluded_source_rows"]
            for key, value in current["coefficients"].items():
                np.testing.assert_allclose(
                    value, previous["coefficients"][key], rtol=1e-5, atol=1e-12
                )
        current = result["thompson_2008"]["diagnostic_refit"]["coefficients"]
        previous = stored["thompson_2008"]["diagnostic_refit"]["coefficients"]
        for key, value in current.items():
            np.testing.assert_allclose(value, previous[key], rtol=1e-5, atol=1e-12)
        assert result["resources"] == stored["resources"]
        print("Diopside source resources and published-curve checks pass")
    else:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
        print(args.output)


if __name__ == "__main__":
    main()
