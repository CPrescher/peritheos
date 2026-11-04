"""Replay the Finger (1979) Fe2SiO4 static ambient-cell reconstruction.

Repository maintenance only. The two-point K0 diagnostic does not replace the
published modulus or recover the original fit intercept, weights or covariance.
No thermal candidate is loaded or fitted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/ferrous-ringwoodite-finger-1979-table2.csv"
RECORD_ID = "ferrous_ringwoodite_finger_1979_bm2_ambient_cell_reconstruction"
DATASET_ID = "ferrous_ringwoodite_finger_1979_table2"
REPORT = "docs/data/finger-1979-ringwoodite-static.json"
QUALIFICATION = (
    "Derived BM2 ambient-cell reconstruction at 296.15 K and 0–4 GPa: measured "
    "V0=558.7(2) A^3 is adopted and fixed; published K0=196(6) GPa and implicit "
    "K0-prime=4 are preserved. V0 error is a measurement esd, not a fit error; "
    "K0 error confidence and covariance are unspecified. The two compressed "
    "states agree within their combined coordinate error widths. The K0-only "
    "equal-pressure diagnostic is conditional on the adopted ambient volume; "
    "its error is RSS/dof-scaled and excludes V0/calibration uncertainty. This "
    "reconstructs the declared static curve, not the authors' original fit "
    "intercept, objective, weighting or pressure reduction. The separate 2011 "
    "thermal candidate remains staged and is not installed or validated here."
)


def observations():
    """Retain all three rows; convert only kbar to GPa for computation."""
    with DATA.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    p = np.array([float(r["pressure_kbar"]) / 10 for r in rows])
    v = np.array([float(r["volume_a3"]) for r in rows])
    sv = np.array([float(r["volume_esd_a3"]) for r in rows])
    return rows, p, v, sv


def bm2(volume, k0=196.0):
    eta = (558.7 / np.asarray(volume)) ** (1 / 3)
    return 1.5 * k0 * (eta**7 - eta**5)


def reproduce():
    rows, p, v, sv = observations()
    eta = (558.7 / v) ** (1 / 3)
    factors = bm2(v, 1.0)
    # Ambient measurement fixes V0; only the compressed rows constrain K0.
    f, observed = factors[1:], p[1:]
    fitted = float(f @ observed / (f @ f))
    residual = f * fitted - observed
    variance = float(residual @ residual / (len(observed) - 1))
    checks = []
    for i, row in enumerate(rows):
        bulk = 196.0 / 2 * (7 * eta[i] ** 7 - 5 * eta[i] ** 5)
        sp = float(row["pressure_esd_kbar"]) / 10 if row["pressure_esd_kbar"] else None
        # Descriptive linearized coordinate-error width, assuming no P–V
        # covariance; parameter/calibration covariance is not propagated.
        width = float(np.hypot(bulk / v[i] * sv[i], sp)) if sp is not None else None
        delta = float(bm2(v[i]) - p[i])
        checks.append(
            {
                "source_row": i + 1,
                "pressure_gpa": float(p[i]),
                "volume_a3": float(v[i]),
                "pressure_error_gpa": sp,
                "volume_esd_a3": float(sv[i]),
                "calculated_pressure_gpa": float(bm2(v[i])),
                "pressure_difference_gpa": delta,
                "combined_coordinate_error_gpa": width,
                "within_combined_coordinate_error": abs(delta) <= width
                if width is not None
                else None,
                "role": "adopted ambient anchor" if i == 0 else "compressed state",
            }
        )
    return {
        "record_identifier": RECORD_ID,
        "source_sha256": hashlib.sha256(DATA.read_bytes()).hexdigest(),
        "temperature_ref_k": 296.15,
        "validity_pressure_gpa": [0.0, 4.0],
        "source_checks": checks,
        "diagnostic_K0_gpa": fitted,
        "diagnostic_K0_standard_error_gpa": float(np.sqrt(variance / (f @ f))),
        "source_compressed_pressure_rms_gpa": float(
            np.sqrt(np.mean((bm2(v[1:]) - p[1:]) ** 2))
        ),
        "diagnostic_compressed_pressure_rms_gpa": float(np.sqrt(np.mean(residual**2))),
        "author_fit_reproduced": False,
        "qualification": QUALIFICATION,
    }


def ledger_outcome(record):
    if record["identifier"] != RECORD_ID:
        raise ValueError("Wrong record for the Finger (1979) static audit")
    result = reproduce()
    fitted = result["diagnostic_K0_gpa"]
    return {
        "status": "source_reconstruction",
        "reason": QUALIFICATION,
        "qualification": QUALIFICATION,
        "reconstruction_kind": "derived_ambient_cell_static_curve",
        "dataset_identifiers": [DATASET_ID],
        "observations": 2,
        "selection": "Two compressed Table 2 Fe2SiO4 rows; the separate measured ambient row fixes V0.",
        "fit_kind": "conditional_K0_diagnostic",
        "objective": "unweighted pressure residuals at measured volumes",
        "weighting": "equal pressure weights; source coordinate errors retained separately",
        "free_parameters": ["K0"],
        "fixed_parameters": ["V0", "K0_prime"],
        "fixed_values": {"V0": 558.7, "K0_prime": 4.0},
        "degrees_of_freedom": 1,
        "solver_success": True,
        "rmse_gpa": result["diagnostic_compressed_pressure_rms_gpa"],
        "published_rmse_gpa": result["source_compressed_pressure_rms_gpa"],
        "author_fit_reproduced": False,
        "parameters": [
            {
                "parameter": "K0",
                "published": 196.0,
                "published_error": 6.0,
                "refit": fitted,
                "refit_error": result["diagnostic_K0_standard_error_gpa"],
                "difference": fitted - 196.0,
                "relative_difference": abs(fitted - 196.0) / 196.0,
                "within_reported_error": abs(fitted - 196.0) <= 6.0,
                "within_combined_2sigma": None,
                "similar": True,
            }
        ],
        "reproduction_report": REPORT,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument(
        "--check", action="store_true", help="check the stored static audit"
    )
    mode.add_argument(
        "--write", action="store_true", help="write the static audit only"
    )
    args = parser.parse_args()
    result = reproduce()
    path = ROOT / REPORT
    if args.check:
        if json.loads(path.read_text()) != result:
            raise SystemExit("Finger (1979) static audit is stale")
        print("Finger (1979) static source checks match the stored audit")
    else:
        payload = (
            json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n"
        )
        if args.write:
            path.write_text(payload, encoding="utf-8")
        else:
            print(payload, end="")


if __name__ == "__main__":
    main()
