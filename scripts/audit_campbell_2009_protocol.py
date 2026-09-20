#!/usr/bin/env python3
"""Separate pressure-scale and regression-objective sensitivities for Fe/FeO.

No input is selected on the basis of proximity to the published coefficients.
Dewaele recalibration is a hypothesis, not Campbell's recovered procedure.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

if __package__:
    from .reconstruct_seagle_2008_fcc import CELL_TO_MOLAR, fcc_pressure
    from .reconstruct_seagle_2008_pressures import hcp_pressure
    from .reproduce_campbell_2009_buffers import (
        SOURCE,
        current_data,
        pressure,
        read_rows,
        seagle_nominal_data,
        volume,
    )
    from .reproduce_iron_source_papers import pressure as iron_pressure
else:
    from reconstruct_seagle_2008_fcc import CELL_TO_MOLAR, fcc_pressure
    from reconstruct_seagle_2008_pressures import hcp_pressure
    from reproduce_campbell_2009_buffers import (
        SOURCE,
        current_data,
        pressure,
        read_rows,
        seagle_nominal_data,
        volume,
    )
    from reproduce_iron_source_papers import pressure as iron_pressure

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/campbell-2009-protocol-audit.json"
PARAMETERS = ["K0_gpa", "gamma0", "q"]


def calibrant_pressure(v, t, scale):
    """Keep original Seagle and hypothetical later iron scales distinct."""
    if scale == "seagle_2006":
        return hcp_pressure(v, t)
    if scale != "dewaele_2006_hypothesis":
        raise ValueError("Unknown pressure scale")
    document = json.loads((ROOT / "peritheos/data/materials/iron.eosmat").read_text())
    record = next(
        r
        for r in document["eos_records"]
        if r["identifier"] == "iron_dewaele_2006_vinet_thermal"
    )
    return iron_pressure(record, v, t)


def dataset(
    material, scale="current_only", include_nominal=False, reconstruct_fcc=False
):
    """Return observations, row IDs, and covariance(P,T) from shared temperature.

    Current-study and nominal rows have unknown P-T covariance, assumed zero.
    Reconstructed P shares the measured T, giving cov(P,T)=dPcal/dT*sigmaT².
    No unreported calibrant-volume uncertainty is filled in.
    """
    data = current_data(material)
    current_rows = read_rows("feo-campbell-2009-table-s2-pvt.csv")
    ids = [
        f"campbell_s2:{r['source_row']}"
        for r in current_rows
        if material == "feo" or r["iron_phase"] == "fcc"
    ]
    covariance_pt = np.zeros(data.shape[1])
    if scale == "current_only":
        return data, ids, covariance_pt
    added, _ = seagle_nominal_data(material)
    rows = read_rows("feo-seagle-2008-supplement-volume-temperature.csv")
    if material == "fe_fcc":
        if scale != "nominal_only" and not reconstruct_fcc:
            raise ValueError("Original individual fcc calibration remains unavailable")
        rows = [r for r in rows if r["iron_phase"] == "fcc"]
        if not reconstruct_fcc:
            return (
                np.concatenate([data, added], axis=1),
                ids + [f"seagle:{r['source_row']}" for r in rows],
                np.zeros(data.shape[1] + len(rows)),
            )
    keep = []
    covariances = []
    for i, row in enumerate(rows):
        if row["iron_phase"] == "hcp" or (
            reconstruct_fcc and row["iron_phase"] == "fcc"
        ):

            def marker_pressure(v, t):
                return (
                    fcc_pressure(v, t)
                    if row["iron_phase"] == "fcc"
                    else calibrant_pressure(v, t, scale)
                )

            v, t = float(row["iron_unit_cell_volume_a3"]), added[1, i]
            added[2, i] = marker_pressure(v, t)
            dpdt = (marker_pressure(v, t + 0.01) - marker_pressure(v, t - 0.01)) / 0.02
            dpdv = (marker_pressure(v + 1e-4, t) - marker_pressure(v - 1e-4, t)) / 2e-4
            sv = row["iron_unit_cell_volume_uncertainty_a3"]
            added[5, i] = (
                np.hypot(dpdv * float(sv), dpdt * added[4, i]) if sv else np.nan
            )
            keep.append(i)
            covariances.append(float(dpdt * added[4, i] ** 2))
        elif include_nominal:
            keep.append(i)
            covariances.append(0.0)
    return (
        np.concatenate([data, added[:, keep]], axis=1),
        ids + [f"seagle:{rows[i]['source_row']}" for i in keep],
        np.concatenate([covariance_pt, covariances]),
    )


def fitted_volumes(p, t, parameters):
    """Invert the actual EOS at each observation, not a linearized residual."""
    return np.array([volume(pi, ti, parameters) for pi, ti in zip(p, t)])


def effective_pressure_variance(material, data, covariance_pt, covariance_pv=None):
    """Residual variance with shared measured variables, at published parameters."""
    v, t, p, sv, st, sp = data
    source = np.array(SOURCE[material])
    dpdv = (pressure(v + 1e-4, t, source) - pressure(v - 1e-4, t, source)) / 2e-4
    dpdt = (pressure(v, t + 0.01, source) - pressure(v, t - 0.01, source)) / 0.02
    variance_p = sp**2 + (dpdv * sv) ** 2 + (dpdt * st) ** 2 - 2 * dpdt * covariance_pt
    if covariance_pv is not None:
        variance_p -= 2 * dpdv * covariance_pv
    return variance_p, dpdv


def fit_protocol(
    material, data, covariance_pt, objective, weighting, covariance_pv=None
):
    """Fit identical free parameters and bounds for every stated protocol."""
    v, t, p, sv, st, sp = data
    source = np.array(SOURCE[material])
    variance_p, dpdv = effective_pressure_variance(
        material, data, covariance_pt, covariance_pv
    )
    if weighting == "unweighted":
        weights = np.ones(len(p))
    elif weighting == "pressure_uncertainty" and objective == "pressure":
        weights = sp
    elif weighting == "effective_variance":
        weights = np.sqrt(variance_p)
        if objective == "volume":
            weights /= abs(dpdv)
    else:
        raise ValueError("Unsupported objective/weight combination")
    if not np.all(np.isfinite(weights) & (weights > 0)):
        raise ValueError("Weighted fits require finite positive input errors")

    def candidate(x):
        parameters = source.copy()
        parameters[[1, 4, 5]] = x
        return parameters

    def residuals(x):
        parameters = candidate(x)
        residual = (
            pressure(v, t, parameters) - p
            if objective == "pressure"
            else fitted_volumes(p, t, parameters) - v
        )
        return residual / weights

    solution = least_squares(
        residuals,
        source[[1, 4, 5]],
        bounds=([20, 0.01, 0.01], [400, 10, 10]),
        x_scale="jac",
        ftol=1e-10,
        xtol=1e-10,
        gtol=1e-10,
    )
    parameters = candidate(solution.x)
    return {
        "observations": len(p),
        "objective": objective,
        "objective_units": "GPa" if objective == "pressure" else "cm3/mol",
        "weighting": weighting,
        "weights_fixed_at_published_parameters": weighting == "effective_variance",
        "parameters": dict(zip(PARAMETERS, solution.x.tolist())),
        "active_parameter_bounds": [
            PARAMETERS[i] for i, active in enumerate(solution.active_mask) if active
        ],
        "sum_squared_scaled_residuals": float(np.sum(solution.fun**2)),
        "rmse_pressure_gpa": float(
            np.sqrt(np.mean((pressure(v, t, parameters) - p) ** 2))
        ),
        "rmse_volume_cm3_mol": float(
            np.sqrt(np.mean((fitted_volumes(p, t, parameters) - v) ** 2))
        ),
        "published_rmse_pressure_gpa": float(
            np.sqrt(np.mean((pressure(v, t, source) - p) ** 2))
        ),
        "published_rmse_volume_cm3_mol": float(
            np.sqrt(np.mean((fitted_volumes(p, t, source) - v) ** 2))
        ),
        "solver_success": bool(solution.success),
        "solver_message": solution.message,
    }


def reproduce():
    groups = {}
    for material, scale, nominal, full_fcc in [
        ("fe_fcc", "current_only", False, False),
        ("fe_fcc", "nominal_only", True, False),
        ("feo", "current_only", False, False),
        ("feo", "seagle_2006", False, False),
        ("feo", "seagle_2006", True, False),
        ("feo", "dewaele_2006_hypothesis", False, False),
        ("feo", "dewaele_2006_hypothesis", True, False),
        ("fe_fcc", "funamori_boehler", False, True),
        ("feo", "seagle_2006", False, True),
        ("feo", "seagle_2006", True, True),
        ("feo", "dewaele_2006_hypothesis", False, True),
        ("feo", "dewaele_2006_hypothesis", True, True),
    ]:
        data, ids, cov = dataset(material, scale, nominal, full_fcc)
        cov_pv = np.zeros(data.shape[1])
        if material == "fe_fcc" and full_fcc:
            # Fe is both the pressure marker and the fitted phase here. Include
            # covariance(P,V) as well as covariance(P,T) for its Seagle rows.
            v, t, _, sv, _, _ = data[:, 15:]
            cell = v / CELL_TO_MOLAR
            dpdv = (
                (fcc_pressure(cell + 1e-4, t) - fcc_pressure(cell - 1e-4, t))
                / 2e-4
                / CELL_TO_MOLAR
            )
            cov_pv[15:] = dpdv * sv**2
        mask = np.all(np.isfinite(data), axis=0) & np.all(data[3:] > 0, axis=0)
        runs = {}
        for selected, subset, subset_cov, subset_cov_pv in [
            ("all_rows", data, cov, cov_pv),
            ("complete_errors", data[:, mask], cov[mask], cov_pv[mask]),
        ]:
            for objective in ["pressure", "volume"]:
                runs[f"{selected}_{objective}_unweighted"] = fit_protocol(
                    material, subset, subset_cov, objective, "unweighted", subset_cov_pv
                )
                if selected == "complete_errors":
                    runs[f"{selected}_{objective}_effective_variance"] = fit_protocol(
                        material,
                        subset,
                        subset_cov,
                        objective,
                        "effective_variance",
                        subset_cov_pv,
                    )
            if selected == "complete_errors":
                runs[f"{selected}_pressure_uncertainty"] = fit_protocol(
                    material, subset, subset_cov, "pressure", "pressure_uncertainty"
                )
        name = f"{material}_{scale}" + ("_plus_nominal" if nominal else "")
        if full_fcc:
            name += "_with_reconstructed_fcc"
        groups[name] = {
            "material": material,
            "seagle_pressure_scale": scale,
            "uses_nominal_seagle_pressures": nominal,
            "uses_reconstructed_fcc_pressures": full_fcc,
            "shared_marker_volume_covariance_included": material == "fe_fcc"
            and full_fcc,
            "row_ids": ids,
            "weighted_row_ids": [key for key, keep in zip(ids, mask) if keep],
            "missing_error_row_ids": [key for key, keep in zip(ids, mask) if not keep],
            "fits": runs,
        }
    original, ids, _ = dataset("feo", "seagle_2006")
    alternative, _, _ = dataset("feo", "dewaele_2006_hypothesis")
    shift = alternative[2, 25:] - original[2, 25:]
    return {
        "format": "peritheos.campbell-2009-protocol-audit",
        "format_version": 1,
        "qualification": "Conditional sensitivity matrix, not recovery of Campbell's unpublished pressure reduction, row selection or regression procedure. Dewaele is an explicitly hypothetical replacement of the earlier hcp scale. Published coefficient values are unchanged; the Fe/FeO source records are deferred from execution.",
        "uncertainty_assumptions": "Measurement-only reconstructed pressure errors. Effective-variance fits include shared temperature covariance for iron marker-derived pressures, and shared pressure-volume covariance when fcc Fe is both marker and fitted phase. Covariance is unknown and assumed zero for current-study and nominal rows. Other cross-variable covariance and EOS systematics are unavailable. Complete-error fits drop missing errors and have matched unweighted controls. Pressure-error-only weighting ignores covariance by definition. Source confidence conventions are not inferred.",
        "source_parameters": {
            m: dict(zip(PARAMETERS, np.array(SOURCE[m])[[1, 4, 5]].tolist()))
            for m in ["fe_fcc", "feo"]
        },
        "source_error_intervals": {
            "fe_fcc": dict(zip(PARAMETERS, [3, 0.04, 0.6])),
            "feo": dict(zip(PARAMETERS, [1.3, 0.04, 0.3])),
        },
        "dewaele_minus_seagle_hcp_pressure": {
            "observations": len(shift),
            "mean_gpa": float(np.mean(shift)),
            "range_gpa": [float(min(shift)), float(max(shift))],
            "rows": [
                {"source_row": key, "difference_gpa": float(dp)}
                for key, dp in zip(ids[25:], shift)
            ],
        },
        "groups": groups,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != encoded:
            raise SystemExit("Campbell protocol audit is stale")
    else:
        OUTPUT.write_text(encoded)
    print(f"Checked {len(result['groups'])} data/pressure-scale selections")


if __name__ == "__main__":
    main()
