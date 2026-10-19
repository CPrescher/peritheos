"""Replay recoverable Sakai calibration branches without fitting calculated S9 grids.

These conditional subset fits omit unrecovered original selections/weights and
some historical constraints. Sequential fits carry fitted upstream coefficients;
published-upstream controls instead hold the Table 2 calibrants fixed. Neither
mode replaces catalog coefficients or establishes original-workflow parity.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos.eos.rt import RydbergStacey
from scripts.reproduce_sakai_2025_gold import pressure, read_rows
from scripts.reproduce_sakai_2025_nine_materials import PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
ARCHIVE = DATA / "sakai-2025-au-workbook-cells.json"
PAIRS = DATA / "sakai-2025-network-paired-volumes.csv"
OUTPUT = ROOT / "docs/data/sakai-2025-network-reconstruction.json"
# Table, target, calibrant, target-volume column, calibrant-volume column.
# S6 has Fe a/c between Fe volume and W volume: W is J, never D.
BRANCHES = (
    (1, "rhenium", "copper", "D", "B"),
    (3, "tungsten", "copper", "D", "B"),
    (2, "platinum", "rhenium", "B", "D"),
    (4, "gold", "tungsten", "B", "D"),
    (5, "molybdenum", "tungsten", "B", "D"),
    (6, "iron", "tungsten", "B", "J"),
    (7, "mgo", "platinum", "D", "B"),
)


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def workbook_pairs():
    archive = json.loads(ARCHIVE.read_text(encoding="utf-8"))
    rows = []
    for table, target, calibrant, tc, cc in sorted(BRANCHES):
        name = next(n for n in archive["sheets"] if n.startswith(f"Table S{table},"))
        cells = archive["sheets"][name]["cells"]
        run = str(table)
        for i in sorted({int(k.lstrip("ABCDEFGHIJKLMNOPQRSTUVWXYZ")) for k in cells}):
            if i < 5:
                continue
            marker = cells.get(f"A{i}", {}).get("value")
            if marker and "#" in str(marker):
                run = str(marker)
            tv, cv = [cells.get(f"{col}{i}", {}) for col in (tc, cc)]
            if not isinstance(tv.get("value"), (int, float)):
                continue
            if not isinstance(cv.get("value"), (int, float)):
                raise ValueError(f"Missing paired volume: {name}:{i}")
            row = {
                "source_table": f"S{table}",
                "source_sheet": name,
                "source_row": i,
                "run": run,
                "target": target,
                "calibrant": calibrant,
            }
            for role, col, cell in (("target", tc, tv), ("calibrant", cc, cv)):
                error = cells.get(f"{chr(ord(col) + 1)}{i}", {})
                row.update(
                    {
                        f"{role}_volume_a3": cell["value"],
                        f"{role}_volume_source_cell": f"{col}{i}",
                        f"{role}_volume_source_token": cell["source_token"],
                        f"{role}_error_a3": error.get("value"),
                        f"{role}_error_source_token": error.get("source_token"),
                    }
                )
            rows.append(row)
    return rows


def historical_nacl_pairs():
    """Only existing measured pairs; omit ambiguous Sata peak-combination inputs."""
    rows = []
    for i, r in enumerate(read_rows("nacl-b2-sakai-2011-table1.csv"), 1):
        rows.append(
            {
                "source": "Sakai2011 Table1",
                "source_row": i,
                "target_volume_a3": float(r["volume_a3"]),
                "calibrant": "platinum",
                "calibrant_volume_a3": float(r["pt_volume_a3"]),
            }
        )
    groups = defaultdict(dict)
    for r in read_rows("dorfman-2012-tables-s1-s6-cocompression.csv"):
        groups[(r["run"], r["observation_index"])][r["material"]] = r
    for key, group in groups.items():
        if "NaCl B2" not in group:
            continue
        # One equality per measured pair; Ne has no Sakai Table 2 scale.
        for symbol, material in (("Pt", "platinum"), ("MgO", "mgo"), ("Au", "gold")):
            if symbol in group:
                rows.append(
                    {
                        "source": "Dorfman2012 " + key[0],
                        "source_row": key[1],
                        "target_volume_a3": float(
                            group["NaCl B2"]["volume_a3_conventional_cell"]
                        ),
                        "calibrant": material,
                        "calibrant_volume_a3": float(
                            group[symbol]["volume_a3_conventional_cell"]
                        ),
                    }
                )
    return rows


def fit_stage(target, rows, calibrants, *, anchor=False):
    """Equal pressure residuals, conditional on fixed upstream scales."""
    v = np.array([r["target_volume_a3"] for r in rows])
    p = np.array(
        [
            float(r["reconstructed_pressure_gpa"])
            if "reconstructed_pressure_gpa" in r
            else float(pressure(r["calibrant_volume_a3"], calibrants[r["calibrant"]]))
            for r in rows
        ]
    )
    published = PARAMETERS[target]
    free = target in ("iron", "nacl_b2")
    if anchor:
        v = np.append(v, 27.25)
        p = np.append(p, 32.3)

    def parameters(x):
        return (
            (float(x[0]), float(x[1]), float(x[2]), published[3])
            if free
            else (published[0], published[1], float(x[0]), published[3])
        )

    x0 = published[:3] if free else [published[2]]
    bounds = (
        ([v.max() * 1.000001, 1e-8, 0], [published[0] * 3, 2000, 20])
        if free
        else ([0], [20])
    )
    fit = least_squares(
        lambda x: pressure(v, parameters(x)) - p,
        x0,
        bounds=bounds,
        x_scale="jac",
        max_nfev=10000,
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    fitted = parameters(fit.x)
    source = pressure(v, published) - p
    refit = pressure(v, fitted) - p
    dof = len(v) - len(fit.x)
    rank = int(np.linalg.matrix_rank(fit.jac))
    if rank != len(fit.x) or np.any(fit.active_mask):
        raise ValueError(f"{target}: singular or bound-limited conditional fit")
    covariance = np.linalg.pinv(fit.jac.T @ fit.jac) * (refit @ refit) / dof
    native_error = float(
        np.max(np.abs(RydbergStacey(*fitted).pressure(v) - pressure(v, fitted)))
    )
    if not fit.success or native_error > 1e-8:
        raise ValueError(f"{target}: failed fit or native cross-check")
    result = {
        "target": target,
        "measured_pairs": len(rows),
        "observations": len(v),
        "anchor": {
            "pressure_gpa": 32.3,
            "volume_a3": 27.25,
            "role": "Matsui2012 reference point from Sakai Table S8; one equal-residual constraint, not a new paired measurement",
        }
        if anchor
        else None,
        "published_parameters": list(published),
        "fitted_parameters": list(fitted),
        "free_parameters": ["V0", "K0", "K0_prime"] if free else ["K0_prime"],
        "conditional_parameter_errors": np.sqrt(np.diag(covariance)).tolist(),
        "conditional_covariance": covariance.tolist(),
        "degrees_of_freedom": dof,
        "covariance_convention": "Residual-scaled, equal-pressure conditional covariance; omits fixed-input, shared-calibrant and measurement-error covariance. Source err. fields do not define these errors.",
        "published_rms_gpa": float(np.sqrt(np.mean(source**2))),
        "refit_rms_gpa": float(np.sqrt(np.mean(refit**2))),
        "published_max_abs_gpa": float(np.max(np.abs(source))),
        "refit_max_abs_gpa": float(np.max(np.abs(refit))),
        "optimizer_success": bool(fit.success),
        "jacobian_rank": rank,
        "active_bounds": fit.active_mask.tolist(),
        "native_max_difference_gpa": native_error,
        "author_workflow_reproduced": False,
        "source_uncertainty_reproduced": False,
        "reductions": [
            {
                **r,
                "calculated_pressure_gpa": float(pp),
                "published_residual_gpa": float(ss),
                "refit_residual_gpa": float(rr),
            }
            for r, pp, ss, rr in zip(rows, p, source, refit)
        ],
    }
    return result


def reproduce():
    from scripts.sakai_2025_historical_inputs import FILES, historical_rows

    pairs = workbook_pairs()
    modes = {}
    for mode in ("published_upstream", "sequential_measured_subsets"):
        scales = dict(PARAMETERS)
        stages = []
        for table, target, calibrant, _, _ in BRANCHES:
            rows = [
                {
                    "source": r["source_sheet"],
                    "source_row": r["source_row"],
                    "target_volume_a3": r["target_volume_a3"],
                    "calibrant": calibrant,
                    "calibrant_volume_a3": r["calibrant_volume_a3"],
                }
                for r in pairs
                if r["source_table"] == f"S{table}"
            ]
            stage = fit_stage(target, rows, scales)
            stage["upstream_parameters"] = {calibrant: list(scales[calibrant])}
            stages.append(stage)
            if mode == "sequential_measured_subsets":
                scales[target] = tuple(stage["fitted_parameters"])
        nacl_rows = historical_nacl_pairs()
        nacl = fit_stage("nacl_b2", nacl_rows, scales, anchor=True)
        nacl["upstream_parameters"] = {
            m: list(scales[m]) for m in {r["calibrant"] for r in nacl_rows}
        }
        stages.append(nacl)
        modes[mode] = stages
    # The source's W -> Re auxiliary constraint changes the ordering once
    # earlier measured inputs are included. NaCl must precede the Fe branch.
    scales = dict(PARAMETERS)
    stages = []
    order = (3, 1, 2, 7, 4, 5, 0, 6)
    for table in order:
        if table == 0:
            target = "nacl_b2"
            rows = historical_nacl_pairs()
        else:
            _, target, calibrant, _, _ = next(b for b in BRANCHES if b[0] == table)
            rows = [
                {
                    "source": r["source_sheet"],
                    "source_row": r["source_row"],
                    "target_volume_a3": r["target_volume_a3"],
                    "calibrant": calibrant,
                    "calibrant_volume_a3": r["calibrant_volume_a3"],
                }
                for r in pairs
                if r["source_table"] == f"S{table}"
            ]
            rows += historical_rows(target)
        stage = fit_stage(target, rows, scales, anchor=table == 0)
        stage["upstream_parameters"] = {
            m: list(scales[m])
            for m in sorted({r["calibrant"] for r in rows})
            if m in scales
        }
        stages.append(stage)
        scales[target] = tuple(stage["fitted_parameters"])
    modes["sequential_with_recovered_history"] = stages
    sources = [
        ARCHIVE,
        DATA / "nacl-b2-sakai-2011-table1.csv",
        DATA / "dorfman-2012-tables-s1-s6-cocompression.csv",
    ]
    sources += [DATA / f for f in FILES if DATA / f not in sources]
    sources += [ROOT / "docs/data/sakai-2025-gold-reconstruction.json"]
    return {
        "audit_date": "2026-10-19",
        "residual_sign": "Pmodel(Vtarget) - Pcalibrant(Vcalibrant)",
        "temperature_k": 300,
        "workbook_measured_pairs": len(pairs),
        "source_hashes": {str(p.relative_to(ROOT)): digest(p) for p in sources},
        "scope": "Eight conditional subset calibration stages; published Cu anchor held fixed. No S9 output is fitted. Historical inputs outside the named pairs are omitted, so coefficient mismatch here is not a failure of the full published fit.",
        "missing_original_workflow": [
            "Original weighting, covariance and unrounded regression state",
            "The first two modes omit earlier-study constraints. The third includes only the explicitly recovered historical rows; original row selection remains unverified.",
            "Original Sakai2018 Pt-Re observations, unrecovered Fe historical sources, and Sata2002 NaCl peak combination remain outside the historical control.",
            "NaCl Sata2002 peak-combination prescription; Matsui reference point enters as one unweighted constraint, not a hard equality",
        ],
        "modes": modes,
        "full_network_parity_established": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    if args.check:
        from scripts.validate_primary_eos_refits import check_ledger_values

        saved = json.loads(OUTPUT.read_text(encoding="utf-8"))
        if saved["source_hashes"] != result["source_hashes"]:
            raise ValueError("Source hashes differ")
        for mode, stages in result["modes"].items():
            for original, current in zip(saved["modes"][mode], stages):
                if original["published_parameters"] != current["published_parameters"]:
                    raise ValueError("Published source parameters differ")
        check_ledger_values(saved, result, calculated=True, absolute_tolerance=1e-6)
    else:
        with PAIRS.open("w", newline="", encoding="utf-8") as stream:
            writer = csv.DictWriter(
                stream, fieldnames=list(workbook_pairs()[0]), lineterminator="\n"
            )
            writer.writeheader()
            writer.writerows(workbook_pairs())
        OUTPUT.write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    for mode, stages in result["modes"].items():
        print(mode)
        for s in stages:
            print(
                s["target"],
                s["observations"],
                s["fitted_parameters"],
                s["refit_rms_gpa"],
            )


if __name__ == "__main__":
    main()
