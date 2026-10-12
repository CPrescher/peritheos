"""Check every packaged historical Au/MgO pair and diagnostic Au recalibration.

python -m scripts.audit_gold_speziale_reduction [--check]

Writes only its own report and residual CSV. Source data, published EOS records,
calibration metadata, shared ledgers and source uncertainty widths are untouched.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos import get_eos_record
from scripts.check_numerical_archive import check_csv
from scripts.reproduce_fei_2007_gold import check_saved
from scripts.reproduce_hirose_2008_gold import source_pressure as hirose2008
from scripts.reproduce_speziale_2001_mgo import pressure as mgo

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/gold-speziale-reduction-audit.json"
RESIDUALS = ROOT / "docs/data/gold-speziale-reduction-residuals.csv"
FILES = {
    "Fei2004": "gold-fei-2004-table1.csv",
    "Hirose2006": "gold-hirose-2006-table1.csv",
    "Hirose2008": "gold-hirose-2008-table1.csv",
}


def observations():
    points = []
    for source, name in FILES.items():
        with (DATA / name).open() as stream:
            rows = list(csv.DictReader(stream))
        for number, row in enumerate(rows, 1):
            if not row["mgo_volume_a3"]:
                continue
            key = (
                "mgo_speziale_pressure_gpa"
                if source == "Hirose2006"
                else "pressure_gpa"
            )
            hot = float(row["temperature_k"]) > 300
            points.append(
                {
                    "source": source,
                    "source_row": number,
                    "run": row["run"],
                    "cycle": row.get("cycle", ""),
                    "au_volume_a3": float(row["volume_a3"]),
                    "mgo_volume_a3": float(row["mgo_volume_a3"]),
                    "mgo_a_angstrom": float(row["mgo_a_angstrom"]),
                    "temperature_k": float(row["temperature_k"]),
                    "reported_pressure_gpa": float(row[key]),
                    "reported_pressure_width_gpa": float(row[key + "_uncertainty"])
                    if row[key + "_uncertainty"]
                    else None,
                    # Fei lattice values are printed to four decimal places;
                    # its CSV float representation sometimes drops trailing zeros.
                    "lattice_half_digit_angstrom": 0.00005
                    if source != "Hirose2008" or row["run"] == "3"
                    else 0.000005,
                    "pressure_half_digit_gpa": 0.005
                    if source == "Fei2004" or (source == "Hirose2008" and not hot)
                    else 0.05,
                }
            )
    seen = {}
    for point in points:
        key = tuple(
            point[k]
            for k in (
                "au_volume_a3",
                "mgo_volume_a3",
                "temperature_k",
                "reported_pressure_gpa",
            )
        )
        identity = f"{point['source']}:{point['source_row']}"
        point["duplicate_of"] = seen.get(key, "")
        seen.setdefault(key, identity)
    return points


def metrics(delta):
    return {
        "rows": len(delta),
        "rmse_gpa": float(np.sqrt(np.mean(delta**2))),
        "mae_gpa": float(np.mean(np.abs(delta))),
        "mean_residual_gpa": float(np.mean(delta)),
        "min_residual_gpa": float(np.min(delta)),
        "max_residual_gpa": float(np.max(delta)),
        "max_abs_residual_gpa": float(np.max(np.abs(delta))),
    }


def arrays(points):
    return tuple(
        np.array([p[k] for p in points])
        for k in (
            "au_volume_a3",
            "mgo_volume_a3",
            "temperature_k",
            "reported_pressure_gpa",
        )
    )


def diagnostic_fit(points, evaluator, initial, scales, names):
    v, vm, t, reported = arrays(points)
    targets = {
        "reported": reported,
        "published_speziale": mgo(vm, t),
        "local_power_candidate": mgo(vm, t, gamma_law="local_power"),
    }
    fitted = {}
    predictions = {}
    for label, target in targets.items():
        result = least_squares(
            lambda z: evaluator(v, t, z * scales) - target,
            initial / scales,
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
            max_nfev=4000,
        )
        values = result.x * scales
        predicted = evaluator(v, t, values)
        predictions[label] = predicted
        fitted[label] = {
            "parameters": dict(zip(names, values.tolist())),
            "solver_success": bool(result.success),
            "published_curve_residual": metrics(evaluator(v, t, initial) - target),
            "refit_residual": metrics(predicted - target),
            "pressure_target_shift": metrics(target - reported),
        }
    for label in ("published_speziale", "local_power_candidate"):
        fitted[label]["refitted_curve_shift_from_reported_refit"] = metrics(
            predictions[label] - predictions["reported"]
        )
    return {
        "source_rows": [f"{p['source']}:{p['source_row']}" for p in points],
        "objective": "equal pressure residual weights; source-exact weights unknown",
        "author_fit_reproduced": False,
        "fits": fitted,
    }


def reproduce():
    points = observations()
    v, vm, t, p = arrays(points)
    unique = np.array([not r["duplicate_of"] for r in points])
    groups = {
        source: np.array([r["source"] == source for r in points]) for source in FILES
    }
    for source in FILES:
        for suffix, thermal_mask in (("hot", t > 300), ("cold", t == 300)):
            mask = groups[source] & thermal_mask
            if np.any(mask):
                groups[f"{source}_{suffix}"] = mask
    groups.update(
        unique_all=unique,
        unique_hot=unique & (t > 300),
        unique_cold=unique & (t == 300),
    )
    conventions = {}
    for gamma in ("eq11", "local_power"):
        for theta in ("integrated", "eq4_local_q", "eq4_q0", "power", "fixed"):
            for tr in (300, 298, 298.15):
                label = f"{gamma}_{theta}_{tr}K"
                # The deliberately inconsistent local-q theta can become huge.
                # expm1 overflow then gives the correct limiting zero integrand.
                with np.errstate(over="ignore"):
                    predicted = mgo(vm, t, gamma_law=gamma, theta_law=theta, tr=tr)
                    a = np.array([r["mgo_a_angstrom"] for r in points])
                    da = np.array([r["lattice_half_digit_angstrom"] for r in points])
                    bounds = [
                        mgo(
                            (a + sign * da) ** 3,
                            t,
                            gamma_law=gamma,
                            theta_law=theta,
                            tr=tr,
                        )
                        for sign in (-1, 1)
                    ]
                assert np.all(np.isfinite(predicted)), label
                dp = np.array([r["pressure_half_digit_gpa"] for r in points])
                inside = (p >= np.minimum(*bounds) - dp) & (
                    p <= np.maximum(*bounds) + dp
                )
                conventions[label] = {
                    name: {
                        **metrics((predicted - p)[mask]),
                        "inside_lattice_pressure_rounding_envelope": int(
                            np.sum(inside[mask])
                        ),
                    }
                    for name, mask in groups.items()
                }
                for point, value in zip(points, predicted):
                    point[label + "_pressure_gpa"] = float(value)
    native = get_eos_record("mgo_speziale_2001_bm3_2").pressure(vm, t)
    assert np.max(np.abs(native - mgo(vm, t))) < 1e-9
    fei = [r for r in points if r["source"] == "Fei2004"]
    hot = [
        r for r in points if r["source"] == "Hirose2008" and r["temperature_k"] > 300
    ]
    cold = [
        r for r in points if r["source"] == "Hirose2008" and r["temperature_k"] == 300
    ]
    fit_specs = {
        "Hirose2008_300K": (
            cold,
            lambda v, t, z: hirose2008(v, t, fit=0, coefficients=z),
            [5.58],
            [1],
            ["K0_prime"],
        ),
        "Hirose2008_fit1": (
            hot,
            lambda v, t, z: hirose2008(v, t, fit=1, coefficients=z),
            [-0.028, 3.179e-5, 1.477e-8],
            [0.01, 1e-5, 1e-8],
            ["dK_dT", "alpha0", "alpha1"],
        ),
        "Hirose2008_fit2": (
            fei + hot,
            lambda v, t, z: hirose2008(v, t, coefficients=z),
            [1.03e-6, 3.95e-10, 1.61e-13, 3.61e-4],
            [1e-6, 1e-10, 1e-13, 1e-4],
            ["beta1", "beta2", "beta3", "a"],
        ),
    }
    refits = {
        name: diagnostic_fit(
            rows, evaluator, np.array(initial), np.array(scales), names
        )
        for name, (rows, evaluator, initial, scales, names) in fit_specs.items()
    }
    report = {
        "scope": "All packaged Fei2004/Hirose2006/Hirose2008 Au-MgO pairs; no calculated Au pressures used as targets",
        "tabulated_pairs": len(points),
        "unique_pairs": int(unique.sum()),
        "unique_hot_pairs": int(np.sum(unique & (t > 300))),
        "unique_cold_pairs": int(np.sum(unique & (t == 300))),
        "excluded": "Ten Hirose2006 Run1 rows have no MgO observations. Other Au datasets use different scales or provide volumes without historical Speziale pressure outputs; they do not identify this reducer.",
        "duplicates": [
            {k: r[k] for k in ("source", "source_row", "duplicate_of")}
            for r in points
            if r["duplicate_of"]
        ],
        "rounding_scope": "Lattice and pressure last-digit intervals only, at fixed printed T and coefficients; not measurement uncertainty or a confidence test",
        "uncertainty": "Source parenthetical widths retained without assigning sigma; no temperature errors or covariance inferred",
        "native_equation_max_difference_gpa": float(
            np.max(np.abs(native - mgo(vm, t)))
        ),
        "source_hashes": {
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in FILES.values()
        },
        "implementation_hashes": {
            name: hashlib.sha256((ROOT / "scripts" / name).read_bytes()).hexdigest()
            for name in (
                "audit_gold_speziale_reduction.py",
                "reproduce_speziale_2001_mgo.py",
                "reproduce_fei_2007_gold.py",
                "reproduce_hirose_2008_gold.py",
            )
        },
        "conventions": conventions,
        "diagnostic_au_refits": refits,
        "refit_scope": "Hirose2008 diagnostics only. Fei EOS refits are paused pending the author's reply; Fei2004 observations remain source inputs to Hirose fit2. Hirose thermal fits retain published Kprime=5.58 and ambient constraints; the RT sensitivity is separate, not propagated. Equal-pressure weights are diagnostic, with no claim of corrected absolute pressures or source-exact covariance.",
        "historical_reducer_identified": False,
        "published_coefficients_changed": False,
    }
    stream = io.StringIO()
    writer = csv.DictWriter(stream, fieldnames=list(points[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(points)
    return report, stream.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report, residuals = reproduce()
    if args.check:
        saved = json.loads(REPORT.read_text())
        check_saved(saved, report)
        # Check tiny dimensional coefficients relatively: the generic archive
        # tolerance would otherwise accept large fractional changes in beta.
        for name, refit in report["diagnostic_au_refits"].items():
            for target, fit in refit["fits"].items():
                for parameter, value in fit["parameters"].items():
                    expected = saved["diagnostic_au_refits"][name]["fits"][target][
                        "parameters"
                    ][parameter]
                    assert np.isclose(expected, value, rtol=1e-4, atol=1e-16)
        numerical = {f"{name}_pressure_gpa" for name in report["conventions"]}
        check_csv(RESIDUALS.read_text(), residuals, numerical, atol=1e-8)
    else:
        REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
        RESIDUALS.write_text(residuals)
    print(
        f"Verified {report['unique_pairs']} unique Au-MgO pairs; published EOS unchanged."
    )


if __name__ == "__main__":
    main()
