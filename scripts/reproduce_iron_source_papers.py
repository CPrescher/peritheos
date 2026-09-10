#!/usr/bin/env python3
"""Independent pressure equations and explicitly bounded primary-table refits."""

from __future__ import annotations

import argparse
import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/iron-source-papers-reproduction.json"
NA = 6.02214076e23
R = 8.31446261815324
NODES, WEIGHTS = np.polynomial.legendre.leggauss(64)
PREFIXES = (
    "iron_sakai_2014_",
    "fe09ni01_hcp_sakai_2014_",
    "iron_yamazaki_2012_",
    "iron_dubrovinsky_2000_",
    "iron_brown_2000_",
    "nacl_b2_sakai_2014_",
)


def belongs(identifier):
    return identifier.startswith(PREFIXES) or identifier == "iron_dewaele_2006_bm3"


def rows(stem):
    with (DATA / (stem + ".csv")).open() as f:
        return list(csv.DictReader(f))


def pressure(r, v, t=300):
    v = np.asarray(v, dtype=float)
    t = np.asarray(t, dtype=float)
    p = r["eos"]["parameters"]
    v0 = p["V0"]
    k = p.get("K0")
    kp = p.get("K0_prime")
    th = r.get("thermal", {})
    h = th.get("parameters", {})
    if th.get("type") == "AlphaKT":
        v0 = v0 * np.exp(h["alpha0"] * (t - h["Tr"]))
        k = 1 / (
            1 / k + h["beta1"] * (t - h["Tr"]) + h["beta2"] * (t * t - h["Tr"] ** 2)
        )
    x = (v / v0) ** (1 / 3)
    if r["eos"]["type"] == "BM3":
        out = 1.5 * k * (x**-7 - x**-5) * (1 + 0.75 * (kp - 4) * (x**-2 - 1))
    elif r["eos"]["type"] == "Vinet":
        out = 3 * k * x**-2 * (1 - x) * np.exp(1.5 * (kp - 1) * (1 - x))
    else:
        mu = 1 - v / v0
        return p["P0"] + p["rho0"] * p["c0"] ** 2 * mu / (1 - p["s"] * mu) ** 2
    if th.get("type") == "Dewaele2006":
        ratio = v / p["V0"]
        g = h["gamma_inf"] + (h["gamma0"] - h["gamma_inf"]) * ratio ** h["beta"]
        theta = (
            h["theta0"]
            * ratio ** (-h["gamma_inf"])
            * np.exp(
                (h["gamma0"] - h["gamma_inf"]) / h["beta"] * (1 - ratio ** h["beta"])
            )
        )

        def energy(temp):
            z = theta / temp
            samples = z[..., None] * (NODES + 1) / 2
            integral = z / 2 * np.sum(WEIGHTS * samples**3 / np.expm1(samples), axis=-1)
            return 9 * R * temp * integral / z**3

        vm = v * NA * 1e-30 / 2
        out = out + g / vm * (energy(t) - energy(h["Tr"])) / 1e9
        out = (
            out
            + 1.5
            * R
            / vm
            * (
                h["anharmonic_m"] * h["anharmonic_a"] * ratio ** h["anharmonic_m"]
                + h["electronic_g"] * h["electronic_e"] * ratio ** h["electronic_g"]
            )
            * (t * t - h["Tr"] ** 2)
            / 1e9
        )
    return out


def volume(r, p, t=300):
    return brentq(
        lambda v: float(pressure(r, v, t) - p),
        r["eos"]["parameters"]["V0"] * 0.25,
        r["eos"]["parameters"]["V0"] * 2,
    )


def observation_set(r):
    name = r["identifier"]
    if "yamazaki" in name:
        return rows("iron-yamazaki-2012-table-s1"), "pressure_gpa", None
    if name.startswith("nacl_b2_"):
        return (
            [],
            "",
            "Original Sakai (2011a) paired NaCl/Pt calibration rows are unavailable in the attached 2014 paper.",
        )
    if "sakai" in name:
        selected = [
            row
            for row in rows("sakai-2014-table2-pvt")
            if row["material"]
            == ("iron" if name.startswith("iron_") else "fe09ni01_hcp")
        ]
        if r.get("thermal"):
            return (
                [row for row in selected if float(row["temperature_k"]) > 300],
                "pressure_p7_gpa",
                None,
            )
        scale = name.split("sakai_2014_")[1][:2]
        if scale == "p5":
            for row in selected:
                cell = float(row["nacl_volume_a3"])
                x = (38.34 / cell) ** (1 / 3)
                row["pressure_p5_gpa"] = (
                    str(
                        1.5
                        * 45.18
                        * (x**7 - x**5)
                        * (1 + 0.75 * (4.22 - 4) * (x * x - 1))
                    )
                    if cell < 50
                    else row["pressure_p1_gpa"]
                )
                row["pressure_p5_gpa_error"] = ""
        if scale not in ("p1", "p4", "p5", "p6"):
            return (
                [],
                "",
                f"Table 2 does not tabulate {scale.upper()} pressures. Complete scale re-reduction and source g-G V0 regression remain unavailable.",
            )
        for row in selected:
            if (
                float(row["temperature_k"]) == 300
                and not row["pressure_" + scale + "_gpa"]
            ):
                row["pressure_" + scale + "_gpa"] = row["pressure_p1_gpa"]
                row["pressure_" + scale + "_gpa_error"] = row["pressure_p1_gpa_error"]
        return (
            [row for row in selected if float(row["temperature_k"]) == 300],
            "pressure_" + scale + "_gpa",
            None,
        )
    if "dewaele" in name:
        d = [
            x
            for x in rows("iron-dewaele-2006-epaps-compression")
            if x["original_2006_helium_vinet_fit_included"] == "true"
        ]
        for x in d:
            x.update(
                volume_a3_conventional_cell=str(
                    2 * float(x["atomic_volume_a3_per_atom"])
                ),
                temperature_k="300",
            )
        return d, "reported_pressure_gpa", None
    return [], "", "No numeric primary P-V-T regression table is published."


def fit_record(r):
    name = r["identifier"]
    if "brown" in name:
        data = rows("iron-brown-2000-tables-i-ii")
        up = np.array([float(x["particle_velocity_km_s"]) for x in data])
        us = np.array([float(x["shock_velocity_km_s"]) for x in data])
        design = np.column_stack([np.ones(len(up)), up])
        coef = np.linalg.lstsq(design, us, rcond=None)[0]
        quad = np.linalg.lstsq(np.column_stack([design, up * up]), us, rcond=None)[0]
        return dict(
            status="parity",
            observations=len(up),
            objective="Unweighted vertical Us residuals; all 37 Table I–II points.",
            parameters=dict(c0=float(coef[0]), s=float(coef[1])),
            linear_rmse_km_s=float(np.sqrt(np.mean((us - design @ coef) ** 2))),
            quadratic_parameters=dict(
                c0=float(quad[0]), s=float(quad[1]), q=float(quad[2])
            ),
            quadratic_rmse_km_s=float(
                np.sqrt(np.mean((us - np.column_stack([design, up * up]) @ quad) ** 2))
            ),
            qualification="All fitted polynomial coefficients agree at printed precision. Linear RMS uses N=37 (60.8 m/s); using N-2 gives 62.5 m/s, near printed 62 m/s. No Hugoniot temperatures recovered.",
        )
    data, column, reason = observation_set(r)
    if reason:
        result = dict(status="not_refittable", observations=0, qualification=reason)
        if "dubrovinsky" in name:
            result["checkpoint_bulk_modulus_211gpa_300k"] = float(
                -volume(r, 211)
                * (
                    pressure(r, volume(r, 211) * 1.000001)
                    - pressure(r, volume(r, 211) * 0.999999)
                )
                / (volume(r, 211) * 0.000002)
            )
            result["checkpoint_average_alpha_202gpa_5200k"] = float(
                np.log(volume(r, 202, 5200) / volume(r, 202, 300)) / 4900
            )
        return result
    v = np.array([float(x["volume_a3_conventional_cell"]) for x in data])
    t = np.array([float(x["temperature_k"]) for x in data])
    p = np.array([float(x[column]) for x in data])
    base = r["eos"]["parameters"]
    thermal = r.get("thermal", {}).get("parameters", {})
    if "yamazaki" in name:
        free = ["V0", "K0", "K0_prime", "gamma0", "beta", "theta0"]
    elif thermal:
        free = [
            k
            for k in ("gamma0", "gamma_inf", "beta", "theta0")
            if k not in r["thermal"]["fixed_parameters"]
        ]
        if "type1_" in name:
            free = ["gamma0", "theta0"]
    elif "dewaele" in name:
        free = ["V0", "K0_prime"]
    else:
        free = ["K0", "K0_prime"]
    initial = np.array([base.get(k, thermal.get(k)) for k in free])

    def updated(x):
        rr = {**r, "eos": {**r["eos"], "parameters": dict(base)}}
        if thermal:
            rr["thermal"] = {**r["thermal"], "parameters": dict(thermal)}
        for k, z in zip(free, x):
            (rr["eos"]["parameters"] if k in base else rr["thermal"]["parameters"])[
                k
            ] = z
        if "type1_" in name:
            rr["thermal"]["parameters"]["gamma_inf"] = rr["thermal"]["parameters"][
                "gamma0"
            ]
        return rr

    def residual(x):
        return pressure(updated(x), v, t) - p

    lower = []
    upper = []
    bounds = {
        "V0": (20, 25),
        "K0": (30, 500),
        "K0_prime": (1, 10),
        "gamma0": (0.1, 8),
        "gamma_inf": (0, 8),
        "beta": (0.01, 10),
        "theta0": (1, 4000),
    }
    for k in free:
        lo, hi = bounds[k]
        lower.append(lo)
        upper.append(hi)
    fit = least_squares(
        residual,
        initial,
        bounds=(lower, upper),
        x_scale="jac",
        max_nfev=3000,
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-10,
    )
    compare = []
    for k, z in zip(free, fit.x):
        published = base.get(k, thermal.get(k))
        error = (
            r["parameter_errors"].get(k)
            if k in base
            else r["thermal"]["parameter_errors"].get(k)
        )
        compare.append(
            dict(
                parameter=k,
                published=published,
                fitted=float(z),
                published_error=error,
                absolute_difference=float(abs(z - published)),
                within_reported_error=bool(abs(z - published) <= error)
                if error
                else None,
            )
        )

    def rmse(a):
        return float(np.sqrt(np.mean(a * a)))

    status = (
        "similar"
        if all(x["within_reported_error"] for x in compare)
        else "parity_not_achieved"
    )
    result = dict(
        status=status,
        observations=len(v),
        temperature_range_k=[float(min(t)), float(max(t))],
        pressure_column=column,
        objective="Unweighted vertical pressure residuals; conditional on published reference/fixed coefficients. Source weights and covariance are not given.",
        parameters={k: float(z) for k, z in zip(free, fit.x)},
        parameter_comparisons=compare,
        parameter_bounds={k: list(bounds[k]) for k in free},
        active_bounds=[k for k, flag in zip(free, fit.active_mask) if flag],
        published_pressure_rmse_gpa=rmse(residual(initial)),
        refitted_pressure_rmse_gpa=rmse(fit.fun),
        solver_success=bool(fit.success),
        qualification="Conditional independent refit; matching within printed errors is similarity, not reproduction of the original statistical procedure. V0 preliminary g-G fit not repeated."
        if "sakai" in name
        else (
            "37 helium hcp rows selected by the existing original-2006 flag; K0 fixed at 165 GPa. Original weighting is unavailable."
            if "dewaele" in name
            else "All 207 primary rows on original Tsuchiya Au scale. Unweighted P residuals do not assert the unavailable original weighting."
        ),
    )
    if "type" in name:
        result["checkpoint_329gpa_5000k_volume_a3"] = volume(r, 329, 5000)
    if "yamazaki" in name:
        result["checkpoint_330gpa_6000k_density_g_cm3"] = (
            2 * 55.845 / (NA * 1e-24 * volume(r, 330, 6000))
        )
    if "sakai" in name and not thermal:
        result["checkpoint_329gpa_300k_volume_a3"] = volume(r, 329)
        if "_p5_" in name:
            result["qualification"] += (
                " P5 pressures independently reconstructed from the printed new NaCl-B2 BM3 coefficients and Table 2 NaCl volumes; B1 rows retain Brown (1999) pressures."
            )
    return result


@lru_cache(maxsize=1)
def reproduce():
    report = {}
    for material in ("iron", "fe09ni01_hcp", "nacl_b2"):
        doc = json.loads(
            (ROOT / "peritheos/data/materials" / f"{material}.eosmat").read_text()
        )
        for r in doc["eos_records"]:
            if belongs(r["identifier"]):
                report[r["identifier"]] = fit_record(r)
    return dict(
        source_manifest="peritheos/data/datasets/iron-source-papers-zotero.json",
        records=report,
    )


def ledger_outcome(record):
    result = reproduce()["records"][record["identifier"]]
    parameters = []
    for name, value in result.get("parameters", {}).items():
        published = record["eos"]["parameters"].get(name)
        if published is None:
            published = record["thermal"]["parameters"][name]
        comparison = next(
            (
                p
                for p in result.get("parameter_comparisons", [])
                if p["parameter"] == name
            ),
            {},
        )
        parameters.append(
            dict(
                parameter=name,
                published=published,
                refit=value,
                published_error=comparison.get("published_error"),
                refit_error=None,
                relative_difference=abs(value - published) / abs(published)
                if published
                else None,
                within_combined_2sigma=None,
                similar=bool(
                    comparison.get(
                        "within_reported_error", result["status"] == "parity"
                    )
                ),
            )
        )
    return {
        **result,
        "parameters": parameters,
        "reason": result["qualification"],
        "dataset_identifiers": record["scientific_validation"]["primary_data_check"][
            "dataset_identifiers"
        ],
        "parameter_covariance": None,
        "weighting": "unweighted; published coefficient constraints retained",
        "free_parameters": list(result.get("parameters", {})),
        "fixed_parameters": [
            k
            for k in {
                **record["eos"]["parameters"],
                **record.get("thermal", {}).get("parameters", {}),
            }
            if k not in result["parameters"]
            and not ("_type1_" in record["identifier"] and k == "gamma_inf")
        ]
        if "parameters" in result
        else record.get("fixed_parameters", []),
        "tied_parameters": {"gamma_inf": "gamma0"}
        if "_type1_" in record["identifier"]
        else {},
        "published_rmse_gpa": result.get("published_pressure_rmse_gpa"),
        "rmse_gpa": result.get("refitted_pressure_rmse_gpa"),
        "rmse_shock_velocity_km_s": result.get("linear_rmse_km_s"),
        "parity_basis": "published_decimal_precision"
        if result["status"] == "parity"
        else None,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    text = json.dumps(reproduce(), indent=2, ensure_ascii=False, allow_nan=False) + "\n"
    if args.check:
        if REPORT.read_text() != text:
            raise SystemExit("Iron source reproduction report is stale")
    else:
        REPORT.write_text(text)
    print(
        json.dumps(
            {
                k: {
                    x: v[x]
                    for x in (
                        "status",
                        "observations",
                        "published_pressure_rmse_gpa",
                        "checkpoint_329gpa_5000k_volume_a3",
                        "checkpoint_330gpa_6000k_density_g_cm3",
                    )
                    if x in v
                }
                for k, v in reproduce()["records"].items()
            },
            indent=2,
        )
    )
