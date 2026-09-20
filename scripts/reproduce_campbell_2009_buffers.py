#!/usr/bin/env python3
"""Independent Campbell (2009) BM3-Debye reproduction and conditional refits.

The numerical reference uses cm3/mol, a vectorized Gauss-Legendre Debye
integral, and scipy inversion; it does not call Peritheos to compute expected
pressures, fit coefficients, or the published buffer-table checkpoints.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
from scipy.constants import Avogadro, R
from scipy.optimize import brentq, least_squares
from scipy.special import roots_legendre

if __package__:
    from . import reconstruct_seagle_2008_pressures as seagle_pressure
else:
    import reconstruct_seagle_2008_pressures as seagle_pressure

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/campbell-2009-reproduction.json"
# Table 1: V0 (cm3/mol), K0 (GPa), K0', theta0 (K), gamma0, q, atom count.
SOURCE = {
    "fe_fcc": (7.076, 133.0, 5.0, 470.0, 1.95, 1.6, 1.0),
    "feo": (12.256, 146.9, 4.0, 380.0, 1.42, 1.3, 2.0),
    "nickel": (6.587, 179.0, 4.3, 415.0, 2.50, 1.0, 1.0),
    "nickel_oxide_b1": (10.973, 190.0, 5.4, 480.0, 1.80, 1.0, 2.0),
}
RECORDS = {key: f"{key}_campbell_2009_bm3_mgd" for key in SOURCE}
S2 = "feo_campbell_2009_table_s2_pvt"
S3 = "ni_nio_campbell_2009_table_s3_pvt"
SEAGLE = "feo_seagle_2008_supplement_volume_temperature"
SELECTION = (
    "All Table S2 rows for FeO; only its 15 explicitly fcc rows for Fe. "
    "For Ni/NiO, retain all 92 heated/MAP Table S3 rows (T>295 K); "
    "exclude its nine 293-295 K quenched DAC rows from the main diagnostic. "
    "This exclusion follows the Section 2.2 stress warning but is an audit "
    "choice: the authors do not publish a regression row mask. All 101 rows "
    "are retained in the dataset and tested in a sensitivity refit."
)
QUALIFICATION = (
    "Conditional reproduction, not an exact recovery of the unpublished "
    "regression protocol. Campbell specifies least squares but no residual "
    "variable, weights, covariance, or final row mask. The main diagnostic "
    "uses equal pressure-residual weights and simultaneous fitting of the "
    "three Table 1 free parameters. Fe/FeO main fits use current-study Table "
    "S2 only: Seagle contributes 81 recoverable V-T pairs but no individual "
    "pressures or NaCl observations. Separate sensitivities reconstruct 65 "
    "hcp-calibrated pressures using benchmark-checked Seagle (2006), compare "
    "its conflicting reference volumes, and optionally retain 16 unresolved "
    "rows as nominal-isobar estimates. A separate protocol audit adds 14 "
    "conditional fcc pressures from Funamori/Boehler/Basinski, leaving only "
    "two missing-iron-volume rows without individual pressures. The exact "
    "historical fcc reduction and combined regression remain unresolved. "
    "These are derived pressures, not original row-level pressure data. "
    "The published Fe/FeO records are deferred from the executable catalog "
    "because the combined-data mismatch remains unresolved. "
    "No refitted coefficients replace the source "
    "records. See literature-reproductions/campbell-2009-buffers.md."
)
NODES, WEIGHTS = roots_legendre(48)
NODES = (NODES + 1.0) / 2.0
WEIGHTS = WEIGHTS / 2.0


def pressure(volume, temperature, parameters):
    """Equation (5), in GPa, for molar volume in cm3/mol and temperature in K."""
    v0, k0, kp, theta0, gamma0, q, n = parameters
    volume = np.asarray(volume, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    ratio = volume / v0
    gamma = gamma0 * ratio**q
    theta = theta0 * ratio ** (-gamma)

    def energy(t):
        x = theta / t
        z = x[..., None] * NODES
        debye = 3 * np.sum(WEIGHTS * NODES**2 * z / np.expm1(z), axis=-1)
        return 3 * n * R * t * debye

    cold = (
        1.5
        * k0
        * (ratio ** (-7 / 3) - ratio ** (-5 / 3))
        * (1 + 0.75 * (kp - 4) * (ratio ** (-2 / 3) - 1))
    )
    return cold + gamma * (energy(temperature) - energy(295.0)) / volume / 1000


def volume(pressure_gpa, temperature, parameters):
    """Independently invert the monotone branch surrounding the source data."""
    return brentq(
        lambda v: float(pressure(v, temperature, parameters)) - pressure_gpa,
        parameters[0] * 0.45,
        parameters[0] * 1.5,
        xtol=1e-12,
    )


def read_rows(name):
    with (DATA / name).open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def current_data(material, include_quenched=False):
    """Return unchanged source observations and their three reported errors."""
    iron = material in {"fe_fcc", "feo"}
    rows = read_rows(
        "feo-campbell-2009-table-s2-pvt.csv"
        if iron
        else "ni-nio-campbell-2009-table-s3-pvt.csv"
    )
    if material == "fe_fcc":
        rows = [r for r in rows if r["iron_phase"] == "fcc"]
    if not include_quenched:
        rows = [r for r in rows if float(r["sample_temperature_k"]) > 295]
    prefix = {
        "fe_fcc": "iron",
        "feo": "feo",
        "nickel": "nickel",
        "nickel_oxide_b1": "nio",
    }[material]
    names = [
        f"{prefix}_molar_volume_cm3_mol",
        "sample_temperature_k",
        "reported_pressure_gpa",
        f"{prefix}_molar_volume_uncertainty_cm3_mol",
        "sample_temperature_uncertainty_k",
        "reported_pressure_uncertainty_gpa",
    ]
    return np.array([[float(row[name]) for name in names] for row in rows]).T


def seagle_nominal_data(material):
    """Explicit sensitivity-only reconstruction from Seagle (2008) Table 1.

    Match each paired V-T row to the nearest Table 1 FeO thermal-expansion
    line. This identifies the seven separated volume clusters without using
    any Campbell fitted coefficients to manufacture calibration pressures.
    """
    # P, sigmaP, V(1620 K), alpha from printed Table 1.
    series = np.array(
        [
            [20, 4, 76.23, 2.1e-5],
            [27, 4, 73.31, 1.6e-5],
            [50, 4, 67.79, 1.2e-5],
            [55, 4, 67.24, 0.9e-5],
            [72, 4, 63.09, 1.0e-5],
            [89, 5, 60.48, 0.50e-5],
            [93, 5, 59.76, 0.54e-5],
        ]
    )
    rows = read_rows("feo-seagle-2008-supplement-volume-temperature.csv")
    values = []
    assignments = []
    for row in rows:
        t = float(row["temperature_k"])
        vo = float(row["feo_unit_cell_volume_a3"])
        predicted = series[:, 2] * np.exp(series[:, 3] * (t - 1620))
        index = int(np.argmin(abs(predicted - vo)))
        p, sp = series[index, :2]
        assignments.append(
            {"source_row": int(row["source_row"]), "nominal_pressure_gpa": float(p)}
        )
        if material == "fe_fcc" and row["iron_phase"] != "fcc":
            continue
        prefix = "iron" if material == "fe_fcc" else "feo"
        v = row[f"{prefix}_unit_cell_volume_a3"]
        sv = row[f"{prefix}_unit_cell_volume_uncertainty_a3"]
        # Missing Fe volume errors remain absent: only unweighted combined fits.
        if not v:
            continue
        scale = Avogadro * 1e-24 / 4
        values.append(
            [float(v) * scale, t, p, float(sv) * scale if sv else np.nan, 150.0, sp]
        )
    return np.array(values).T, assignments


def seagle_reconstructed_data(reference="density_text", include_nominal=False):
    """FeO observations using independently calibrated hcp Fe when available."""
    data, assignments = seagle_nominal_data("feo")
    rows = seagle_pressure.reconstructed_rows(assignments)
    keep = []
    column = (
        "reconstructed_pressure_gpa"
        if reference == "density_text"
        else "table_2_reference_pressure_gpa"
    )
    if reference not in {"density_text", "table_2"}:
        raise ValueError("Unknown Seagle reference-volume convention")
    for i, row in enumerate(rows):
        if row[column] != "":
            data[2, i] = row[column]
            uncertainty = row["propagated_measurement_uncertainty_gpa"]
            # No calibrated uncertainty/covariance is available for weighted fits.
            data[5, i] = float(uncertainty) if uncertainty != "" else np.nan
            keep.append(i)
        elif include_nominal:
            keep.append(i)
    return data[:, keep]


def fit(
    material,
    weighting="unweighted",
    include_quenched=False,
    seagle=False,
    seagle_reference=None,
    include_nominal=False,
):
    """Simultaneously fit only Table 1's free coefficients; keep adopted inputs."""
    data = current_data(material, include_quenched)
    if seagle_reference is not None:
        if material != "feo" or seagle or weighting != "unweighted":
            raise ValueError(
                "Reconstructed Seagle inputs require an unweighted FeO fit"
            )
        data = np.concatenate(
            [data, seagle_reconstructed_data(seagle_reference, include_nominal)], axis=1
        )
    if seagle:
        additional, _ = seagle_nominal_data(material)
        data = np.concatenate([data, additional], axis=1)
    v, t, p, sv, st, sp = data
    source = np.array(SOURCE[material])
    iron = material in {"fe_fcc", "feo"}
    indices = [1, 4, 5] if iron else [1, 2, 4]
    names = (
        ["rt_eos.K0", "gamma0", "q"]
        if iron
        else ["rt_eos.K0", "rt_eos.K0_prime", "gamma0"]
    )
    scale = np.ones(len(p))
    if weighting == "pressure_uncertainty":
        scale = sp
    elif weighting == "effective_variance":
        dpdv = (pressure(v + 1e-4, t, source) - pressure(v - 1e-4, t, source)) / 2e-4
        dpdt = (pressure(v, t + 0.01, source) - pressure(v, t - 0.01, source)) / 0.02
        scale = np.sqrt(sp**2 + (dpdv * sv) ** 2 + (dpdt * st) ** 2)

    def residuals(x):
        candidate = source.copy()
        candidate[indices] = x
        return (pressure(v, t, candidate) - p) / scale

    result = least_squares(
        residuals,
        source[indices],
        bounds=([20, 0.01, 0.01], [400, 10, 10]),
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
    )
    dof = len(p) - len(indices)
    variance = float(np.sum(result.fun**2) / dof)
    covariance = np.linalg.pinv(result.jac.T @ result.jac) * variance
    errors = np.sqrt(np.diag(covariance))
    params = dict(zip(names, result.x.tolist()))
    out = {
        "observations": len(p),
        "free_parameters": names,
        "parameters": params,
        "standard_errors": dict(zip(names, errors.tolist())),
        "covariance": covariance.tolist(),
        "covariance_convention": "Residual-scaled local inverse-Jacobian covariance; audit-only, not source errors.",
        "weighting": weighting,
        "objective": "pressure_residuals",
        "weights_fixed_at_published_parameters": weighting == "effective_variance",
        "include_quenched": include_quenched,
        "seagle_nominal_isobars": seagle,
        "observed_pressure_range_gpa": [float(p.min()), float(p.max())],
        "observed_temperature_range_k": [float(t.min()), float(t.max())],
        "published_rmse_gpa": float(
            np.sqrt(np.mean((pressure(v, t, source) - p) ** 2))
        ),
        "rmse_gpa": float(np.sqrt(np.mean((result.fun * scale) ** 2))),
        "chi_square": float(np.sum(result.fun**2)),
        "degrees_of_freedom": dof,
        "reduced_chi_square": variance,
        "solver_success": bool(result.success),
        "solver_message": result.message,
    }
    if seagle_reference is not None:
        out["active_parameter_bounds"] = [
            names[i] for i, active in enumerate(result.active_mask) if active
        ]
        out["covariance_convention"] += (
            " Boundary-active fits do not support symmetric parameter intervals."
        )
        out["seagle_pressure_reconstruction"] = {
            "reference_volume": seagle_reference,
            "hcp_individual_rows": 65,
            "nominal_only_rows": 16 if include_nominal else 0,
            "pressure_uncertainties_used_as_weights": False,
        }
    return out


def fit_result(material):
    """Adapter for the repository's common coefficient-comparison gate."""
    out = fit(material)
    result = SimpleNamespace(
        free_parameters=out["free_parameters"],
        parameters=out["parameters"],
        standard_errors=out["standard_errors"],
    )
    return out, result


def buffer_checks():
    """Check independent official tabulated outputs, never fit them as observations."""
    rows = read_rows("campbell-2009-tables-s4-s5-buffer-results.csv")
    checks = []
    for table, metal, oxide, p, t in [
        ("Table S4", "fe_fcc", "feo", 10, 1500),
        ("Table S4", "fe_fcc", "feo", 10, 2400),
        ("Table S5", "nickel", "nickel_oxide_b1", 10, 1500),
        ("Table S5", "nickel", "nickel_oxide_b1", 60, 1500),
        ("Table S5", "nickel", "nickel_oxide_b1", 60, 2400),
    ]:
        values = {
            float(r["pressure_gpa"]): float(r["log10_fugacity"])
            for r in rows
            if r["source_table"] == table and float(r["temperature_k"]) == t
        }
        # Five-point derivative, h=1 GPa; compare with three-point to quantify h error.
        derivative = (
            values[p - 2] - 8 * values[p - 1] + 8 * values[p + 1] - values[p + 2]
        ) / 12
        delta_source = derivative * R * t / (0.8686 * 1000)
        three = (values[p + 1] - values[p - 1]) / 2 * R * t / (0.8686 * 1000)
        delta_model = volume(p, t, SOURCE[oxide]) - volume(p, t, SOURCE[metal])
        checks.append(
            {
                "source_table": table,
                "pressure_gpa": p,
                "temperature_k": t,
                "source_delta_volume_cm3_mol": delta_source,
                "model_delta_volume_cm3_mol": delta_model,
                "difference_cm3_mol": delta_model - delta_source,
                "three_vs_five_point_difference_cm3_mol": three - delta_source,
                "tolerance_cm3_mol": 0.02,
            }
        )
    return checks


def reproduce():
    fits = {}
    for material, record in RECORDS.items():
        main = fit(material)
        main["selection"] = SELECTION
        main["qualification"] = QUALIFICATION
        main["dataset_identifiers"] = [S2 if material in {"fe_fcc", "feo"} else S3]
        main["sensitivity"] = {
            w: fit(material, weighting=w)
            for w in ["pressure_uncertainty", "effective_variance"]
        }
        if material in {"fe_fcc", "feo"}:
            main["sensitivity"]["with_seagle_nominal_isobars"] = fit(
                material, seagle=True
            )
            if material == "feo":
                for reference in ["density_text", "table_2"]:
                    for nominal in [False, True]:
                        key = f"seagle_hcp_{reference}" + (
                            "_plus_nominal" if nominal else ""
                        )
                        main["sensitivity"][key] = fit(
                            material,
                            seagle_reference=reference,
                            include_nominal=nominal,
                        )
        else:
            main["sensitivity"]["all_101_including_quenched"] = fit(
                material, include_quenched=True
            )
        fits[record] = main
    return {
        "format": "peritheos.campbell-2009-reproduction",
        "format_version": 1,
        "record_refits": fits,
        "buffer_checkpoints": buffer_checks(),
        "seagle_nominal_isobar_assignments": seagle_nominal_data("feo")[1],
        "seagle_pressure_calibration_checks": seagle_pressure.calibration_checks(),
        "seagle_pressure_series_checks": seagle_pressure.series_checks(
            seagle_pressure.reconstructed_rows(seagle_nominal_data("feo")[1])
        ),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
    rows = seagle_pressure.reconstructed_rows(
        result["seagle_nominal_isobar_assignments"]
    )
    pressure_csv = seagle_pressure.encode_rows(rows)
    if args.check:
        if not OUTPUT.exists() or OUTPUT.read_text() != encoded:
            raise SystemExit("Campbell reproduction is stale")
        if (
            not seagle_pressure.OUTPUT.exists()
            or seagle_pressure.OUTPUT.read_text() != pressure_csv
        ):
            raise SystemExit("Seagle reconstructed pressures are stale")
    else:
        OUTPUT.write_text(encoded)
        seagle_pressure.OUTPUT.write_text(pressure_csv)
    for record, row in result["record_refits"].items():
        print(record, row["observations"], row["parameters"])


if __name__ == "__main__":
    main()
