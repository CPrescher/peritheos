"""Repository-only diagnostic: split BM3 fits with the archived starting cell.

Run with ``python -m scripts.check_prescher_2012_starting_volume``.
Preserves existing recovery reports and catalog coefficients. The stored cell
error has no recovered confidence convention and is used only for sensitivity.
Use ``--v0 155.08 --output <separate-report.json>`` for an additional fixed-V0
trial; its uncertainty is not propagated.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import scipy

from scripts.recover_prescher_2012_fe3c import ROOT, TABLE, fit, polynomial_pressure_fit
from scripts.recover_prescher_2012_fe3c_splits import TARGETS, compare_fit

CELL = ROOT / "peritheos/data/research/cementite/prescher-2012-starting-cell.json"
OWN = ROOT / "peritheos/data/research/cementite/prescher-2012-own-shelx-cells.json"
PREVIOUS = ROOT / "docs/data/prescher-2012-fe3c-recovery.json"
REPORT = ROOT / "docs/data/prescher-2012-starting-volume-split-check.json"


def linear_fixed_volume_check(p, v, v0, sigma):
    """Independent weighted linear solution in K0 and K0*(K0_prime-4)."""
    eta = (v0 / v) ** (1 / 3)
    a = 1.5 * (eta**7 - eta**5)
    design = np.column_stack([a, 0.75 * a * (eta**2 - 1)]) / sigma
    beta = np.linalg.lstsq(design, p / sigma, rcond=None)[0]
    k, q = beta
    transform = np.array([[1, 0], [-q / k**2, 1 / k]])
    covariance = transform @ np.linalg.inv(design.T @ design) @ transform.T
    return np.array([k, 4 + q / k]), covariance


def replay(requested_fixed_v0=None):
    cell = json.loads(CELL.read_text())
    previous = json.loads(PREVIOUS.read_text())
    v0 = cell["cell"]["volume_a3"]
    width = cell["stored_cell_errors"]["volume_a3"]
    old_v0 = previous["source_drawn_curve_reconstruction"]["parameters"]["V0"]
    table = np.genfromtxt(TABLE, delimiter=",", names=True)
    working = np.array([r["volume_a3"] for r in json.loads(OWN.read_text())["rows"]])
    p, v = table["pressure_gpa"], table["volume_a3"]
    definitions = [
        (
            "low_previous_protocol",
            "low",
            p < 24,
            v,
            0.2,
            "Final Table S1 rows 1-10; includes 22.6 GPa, as in previous conditional replay",
        ),
        (
            "high_previous_protocol",
            "high",
            p >= 22,
            working,
            0.1,
            "Archived own/SHELX p10-p19 volumes paired by run index with final Table S1 pressures; working p16 differs from final table",
        ),
        (
            "low_literal_below_22",
            "low",
            p < 22,
            v,
            0.2,
            "Final Table S1 rows 1-9, strictly below 22 GPa",
        ),
        (
            "high_final_table",
            "high",
            p >= 22,
            v,
            0.1,
            "Final Table S1 rows 10-19, P>=22 GPa",
        ),
    ]
    cases = []
    for name, side, mask, volumes, sigma, description in definitions:
        ps, vs = p[mask], volumes[mask]
        bases = {
            "previous_drawn_curve_fixed": old_v0,
            "archived_starting_cell_fixed": v0,
        }
        if requested_fixed_v0 is not None:
            bases["requested_ambient_volume_fixed"] = requested_fixed_v0
        if name.endswith("previous_protocol"):
            bases |= {
                "starting_cell_minus_stored_width_fixed": v0 - width,
                "starting_cell_plus_stored_width_fixed": v0 + width,
                "free_V0_initialized_at_starting_cell": None,
            }
        for basis, fixed in bases.items():
            result = fit(
                ps,
                vs,
                "pressure_widths",
                fixed,
                np.full(len(ps), sigma),
                initial=[v0, 160, 6],
            )
            if fixed is not None:
                check, covariance = linear_fixed_volume_check(ps, vs, fixed, sigma)
                values = np.array([result["parameters"][k] for k in ("K0", "K0_prime")])
                np.testing.assert_allclose(values, check, rtol=1e-6, atol=1e-6)
                np.testing.assert_allclose(
                    result["unscaled_inverse_normal_matrix"],
                    covariance,
                    rtol=2e-5,
                    atol=1e-7,
                )
                verification = "Independent weighted linear BM3 solution and transformed covariance agree"
            else:
                check = polynomial_pressure_fit(ps, vs)
                np.testing.assert_allclose(
                    list(result["parameters"].values()), check, rtol=2e-5, atol=2e-5
                )
                verification = "Independent three-term pressure polynomial gives the same free-V0 optimum"
            cases.append(
                {
                    "selection": name,
                    "side": side,
                    "selection_description": description,
                    "table_row_indices_1based": (np.flatnonzero(mask) + 1).tolist(),
                    "pressure_gpa": ps.tolist(),
                    "volume_a3": vs.tolist(),
                    "V0_basis": basis,
                    "initial_V0_a3": v0,
                    "assumed_absolute_pressure_sigma_gpa": sigma,
                    "verification": verification,
                    **result,
                    "comparisons": compare_fit(result, TARGETS[side], True),
                }
            )
    return {
        "format": "peritheos.prescher-2012-starting-volume-split-check",
        "date": "2026-10-30",
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "purpose": "Requested diagnostic of the recovered starting-cell volume; no catalog or Studio changes",
        "input_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in (TABLE, OWN, CELL, PREVIOUS)
        },
        "archived_starting_volume_a3": v0,
        "stored_volume_error_a3_unknown_convention": width,
        "previous_adopted_volume_a3": old_v0,
        "requested_fixed_volume_a3": requested_fixed_v0,
        "targets_K0_K0prime_errorK0_errorK0prime": TARGETS,
        "qualifications": [
            "Only V0 is changed in each fixed-V0 comparison; observations, row selection and pressure weights are preserved.",
            "Absolute pressure sigmas 0.2 GPa (low) and 0.1 GPa (high) are the previous conditional assumptions, not recovered author settings. Volume errors are omitted.",
            "Absolute conditional one-sigma errors and residual-scaled errors are both retained. Source error confidence is unknown.",
            "The archived V0 error is not propagated as an established Gaussian standard deviation. Plus/minus that stored width is a sensitivity check only; fixed-V0 errors condition on exact V0.",
            "Using the recovered volume solely as an optimizer starting value is a separate free-V0 check, with three fitted EOS parameters.",
            "Pressure/temperature and the original EOS role of the starting-cell entry remain unconfirmed; no ambient-pressure observation is invented or added.",
            "Printed agreement uses the same half-up intervals as the earlier split recovery, separately for central coefficients and errors.",
        ],
        "cases": cases,
        "author_workflow_reproduced": False,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--v0", type=float, help="Additional fixed reference volume in A^3"
    )
    parser.add_argument("--output", type=Path, help="Separate output report")
    args = parser.parse_args()
    if args.v0 is not None and (
        not np.isfinite(args.v0) or args.v0 <= 0 or args.output is None
    ):
        parser.error("--v0 requires a positive finite value and a separate --output")
    output = args.output or REPORT
    if args.v0 is not None and output.resolve() == REPORT.resolve():
        parser.error("Use a separate --output to preserve the previous diagnostic")
    report = replay(args.v0)
    output.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(output)
    for case in report["cases"]:
        c = case["comparisons"]["input_widths_as_absolute_1sigma"]
        print(
            case["selection"],
            case["V0_basis"],
            c["K0_K0prime_errorK0_errorK0prime"],
            c["matches_printed_digits"],
        )
