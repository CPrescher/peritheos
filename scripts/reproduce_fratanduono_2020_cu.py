"""Attempt the Cu isentrope-to-298 K reduction using cited primary inputs.

Run ``python -m scripts.reproduce_fratanduono_2020_cu`` to regenerate the
qualified audit. No output replaces the published catalog parameters.
Optional ``--extract-supplement PDF`` requires pdfplumber and transcribes the
source table and vector curves before reconstruction.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, quad
from scipy.interpolate import PchipInterpolator
from scipy.optimize import least_squares

from peritheos.eos.rt import Vinet3

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data"
PREFIX = "fratanduono-2020-cu"
RHO0 = 8.939
TR = 298.0
THETA0 = 343.5
GAMMA0 = 2.0
GAMMA_INF = 1.41
GAMMA_EXPONENT = 13.6  # Kraus eta, unrelated to Vinet3 eta.
SWITCH_DENSITY = RHO0 / 0.64
R_SPECIFIC = 8.31446261815324 / 0.063546  # J/(kg K)
PUBLISHED_S = np.array([138.9, 6.05, 2.53, 1.34])
PUBLISHED_T = np.array([133.6, 6.29, 2.06, 1.65])
PUBLISHED_T_ERRORS = np.array([0.8, 0.8, 0.4, 0.6])


def _write_csv(path, header, rows):
    with path.open("w", newline="", encoding="utf-8") as stream:
        writer = csv.writer(stream, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)


def extract_supplement(pdf: Path) -> None:
    """Extract the exact audited PDF; fail rather than select a different curve."""
    import re

    import pdfplumber

    digest = hashlib.sha256(pdf.read_bytes()).hexdigest()
    record = json.loads((ROOT / "peritheos/data/materials/copper.eosmat").read_text())
    source = next(
        r
        for r in record["eos_records"]
        if r["identifier"] == "copper_fratanduono_2020_vinet3_298k"
    )
    expected = source["scientific_validation"]["primary_source_check"][
        "source_artifacts"
    ][1]["sha256"]
    if digest != expected:
        raise ValueError("Supplement PDF hash differs from the visually audited source")
    with pdfplumber.open(pdf) as document:
        text = (
            document.pages[6]
            .extract_text()
            .replace("TABLES", "TABLE S")
            .split("TABLE S1.")[1]
            .split("TABLE S2.")[0]
        )
        matches = re.findall(r"([\d.]+)\(±([\d.]+)\)\s+([\d.]+)", text)
        rows = sorted((float(p), float(e), float(rho)) for p, e, rho in matches)
        if len(rows) != 100:
            raise ValueError(f"Expected 100 Table S1 rows, found {len(rows)}")
        _write_csv(
            DATA / f"{PREFIX}-table-s1-isentrope.csv",
            ["pressure_gpa", "pressure_uncertainty_gpa", "density_g_cm3"],
            rows,
        )
        page = document.pages[4]
        # Tick marks, not text-label centers. Affine calibration uses all major ticks.
        x_ticks = sorted(
            line["x0"]
            for line in page.lines
            if 200 < line["x0"] < 475
            and abs(line["top"] - 326.146373) < 0.001
            and line["bottom"] - line["top"] > 6
        )
        y_ticks = sorted(
            (
                line["top"]
                for line in page.lines
                if abs(line["x0"] - 186.969206) < 0.001 and 80 < line["top"] < 325
            ),
            reverse=True,
        )
        if len(x_ticks) != 5 or len(y_ticks) != 6:
            raise ValueError("Unexpected Figure S2 axis ticks")
        ax, bx = np.polyfit(x_ticks, [10, 15, 20, 25, 30], 1)
        ay, by = np.polyfit(y_ticks, [0, 500, 1000, 1500, 2000, 2500], 1)
        for name, color, count in [
            ("hugoniot", (0.4, 0.4, 0.4), 903),
            ("blue", (0.0, 0.0, 1.0), 1000),
        ]:
            curves = [
                c
                for c in page.curves
                if c["stroking_color"] == color
                and len(c["pts"]) == count
                and (name == "hugoniot" or c["dash"][0] == [])
            ]
            if len(curves) != 1:
                raise ValueError(f"Ambiguous {name} vector curve")
            pts = np.asarray(curves[0]["pts"])
            _write_csv(
                DATA / f"{PREFIX}-figure-s2-{name}.csv",
                ["density_g_cm3", "pressure_gpa"],
                zip(pts[:, 0] * ax + bx, pts[:, 1] * ay + by),
            )
        metadata = {
            "source_url": "https://journals.aps.org/prl/supplemental/10.1103/PhysRevLett.124.015701",
            "source_sha256": digest,
            "table": "S1, page S7; all 100 rows sorted by pressure, no rescaling",
            "figure": "S2, page S5; gray dashed Our Hugoniot Fit and solid blue curve",
            "x_ticks_pdf_points": x_ticks,
            "y_ticks_pdf_points": y_ticks,
            "density_affine": [ax, bx],
            "pressure_affine": [ay, by],
            "max_tick_residual_density_g_cm3": float(
                np.max(np.abs(np.asarray(x_ticks) * ax + bx - np.arange(10, 31, 5)))
            ),
            "max_tick_residual_pressure_gpa": float(
                np.max(np.abs(np.asarray(y_ticks) * ay + by - np.arange(0, 2501, 500)))
            ),
            "qualification": "Vector digitization retains graphical quantization. These are plotted curves, not original fit observations; neither curve is relabeled as measured 298 K data.",
        }
        (DATA / f"{PREFIX}-extraction.json").write_text(
            json.dumps(metadata, indent=2) + "\n"
        )


def kraus_gamma(density):
    return (
        GAMMA_INF
        + (GAMMA0 - GAMMA_INF) * (RHO0 / np.asarray(density)) ** GAMMA_EXPONENT
    )


def kraus_temperature_factor(density):
    """Analytic integral of gamma d(log rho), also theta/theta0."""
    x = RHO0 / np.asarray(density)
    return x ** (-GAMMA_INF) * np.exp(
        (GAMMA0 - GAMMA_INF) / GAMMA_EXPONENT * (1 - x**GAMMA_EXPONENT)
    )


def debye_energy(temperature, theta):
    """Thermal energy J/kg; zero-point energy cancels at equal density."""
    z = theta / temperature
    integral = quad(lambda t: t**3 / np.expm1(t), 0, z, epsabs=1e-12)[0]
    return 9 * R_SPECIFIC * temperature * integral / z**3


def thermal_correction(density, gamma, factor):
    energy = np.array(
        [
            debye_energy(TR * f, THETA0 * f) - debye_energy(TR, THETA0 * f)
            for f in np.atleast_1d(factor)
        ]
    )
    return np.asarray(gamma) * np.asarray(density) * 1000 * energy / 1e9


def reconstruct(grid_points=8193):
    """Return source rows and qualified candidate reductions, without extrapolation."""
    table = np.loadtxt(
        DATA / f"{PREFIX}-table-s1-isentrope.csv", delimiter=",", skiprows=1
    )
    hug = np.loadtxt(
        DATA / f"{PREFIX}-figure-s2-hugoniot.csv", delimiter=",", skiprows=1
    )
    # Preserve printed ambient rho=8.938 in the source file, but exclude that
    # rounded row from fitting; thermodynamic integration uses prescribed rho0.
    compressed = table[table[:, 0] > 0]
    pressure, _, density = compressed.T
    isentrope = PchipInterpolator(
        np.r_[RHO0, density], np.r_[0.0, pressure], extrapolate=False
    )
    hugoniot = PchipInterpolator(hug[:, 0], hug[:, 1], extrapolate=False)
    selected = density <= hug[-1, 0]
    rho = density[selected]
    grid = np.linspace(RHO0, hug[-1, 0], grid_points)
    # P [GPa] / rho [g/cm3] is 1e6 J/kg; common factors cancel in gamma.
    es = PchipInterpolator(
        grid, cumulative_trapezoid(isentrope(grid) / grid**2, grid, initial=0)
    )
    high_grid = np.linspace(SWITCH_DENSITY, hug[-1, 0], grid_points)
    ph = hugoniot(high_grid)
    eh = 0.5 * ph * (1 / RHO0 - 1 / high_grid)
    gh = (ph - isentrope(high_grid)) / (high_grid * (eh - es(high_grid)))
    high_factor = kraus_temperature_factor(SWITCH_DENSITY) * np.exp(
        cumulative_trapezoid(gh / high_grid, high_grid, initial=0)
    )
    low = rho <= SWITCH_DENSITY
    gamma = kraus_gamma(rho)
    factor = kraus_temperature_factor(rho)
    gamma[~low] = np.interp(rho[~low], high_grid, gh)
    factor[~low] = np.interp(rho[~low], high_grid, high_factor)
    correction = thermal_correction(rho, gamma, factor)
    return table, {
        "density": rho,
        "isentrope": pressure[selected],
        "gamma": gamma,
        "factor": factor,
        "correction": correction,
        "isotherm": pressure[selected] - correction,
        "low": low,
        "hugoniot_max_density": float(hug[-1, 0]),
    }


def fit_diagnostic(density, pressure, initial, sigma=None):
    """Diagnostic fit; no inferred source weighting or coefficient covariance."""
    x = RHO0 / density
    scale = 1.0 if sigma is None else sigma
    fit = least_squares(
        lambda a: (Vinet3(1.0, *a).pressure(x) - pressure) / scale,
        initial,
        bounds=([1.0, -np.inf, -np.inf, -np.inf], np.inf),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
        max_nfev=5000,
    )
    return {
        "parameters": dict(zip(["K0", "eta", "beta", "psi"], map(float, fit.x))),
        "pressure_rmse_gpa": float(np.sqrt(np.mean((fit.fun * scale) ** 2))),
        "objective": "unweighted pressure residuals"
        if sigma is None
        else "Table S1 isentrope pressure-error weighted residuals; source fitting protocol unconfirmed",
        "solver_success": bool(fit.success),
        "parameter_covariance": None,
    }


def reproduce():
    table, r = reconstruct()
    low = r["low"]
    ps = Vinet3(1.0, *PUBLISHED_S)
    pt = Vinet3(1.0, *PUBLISHED_T)
    x = RHO0 / r["density"]
    raw = table[table[:, 0] > 0]
    blue = np.loadtxt(DATA / f"{PREFIX}-figure-s2-blue.csv", delimiter=",", skiprows=1)
    bpressure = np.interp(raw[:, 2], blue[:, 0], blue[:, 1])
    observed_correction = raw[:, 0] - pt.pressure(RHO0 / raw[:, 2])
    result = {
        "status": "partial_thermal_reconstruction_without_parameter_parity",
        "record_identifier": "copper_fratanduono_2020_vinet3_298k",
        "sources": {
            "fratanduono": "https://doi.org/10.1103/PhysRevLett.124.015701",
            "kraus": "https://doi.org/10.1103/PhysRevB.93.134105",
            "kraus_pdf_sha256": "6602e5b65ada728b75dd542e90abb6c4d92dded978ed2c946496325a64a83a3b",
            "kraus_locations": [
                "Section III.B.1 Eqs. (11)-(12)",
                "Section III.B.3 theta0=343.5 K",
                "Section III.B.4 Eq. (16) and Table I",
            ],
        },
        "inputs": {
            "rho0_g_cm3": RHO0,
            "temperature_k": TR,
            "theta0_k": THETA0,
            "gamma0": GAMMA0,
            "gamma_infinity": GAMMA_INF,
            "gamma_density_exponent": GAMMA_EXPONENT,
            "molar_mass_g_mol": 63.546,
        },
        "selection": {
            "source_table_rows": len(table),
            "compressed_source_rows": len(raw),
            "kraus_branch_rows": int(low.sum()),
            "candidate_hugoniot_branch_rows": int((~low).sum()),
            "unreconstructed_high_density_rows": int(len(raw) - len(low)),
            "max_hugoniot_density_g_cm3": r["hugoniot_max_density"],
            "reconstructed_source_isentrope_range_gpa": [
                float(r["isentrope"][0]),
                float(r["isentrope"][-1]),
            ],
        },
        "qualifications": [
            "A refit of the published 298 K EOS is not possible with the current information. Fratanduono explicitly describes stress-to-isentrope corrections, but does not explicitly specify the Debye isentrope-to-298 K calculation used here. This is an assumed candidate reduction, not recovered original 298 K data; pressure agreement establishes neither coefficient parity nor a superior EOS.",
            "Debye correction follows cited Kraus Eq. (16). Carrying theta0 and the Debye caloric model into the 2020 reduction is an explicitly documented inheritance assumption, not a separately printed 2020 parameter.",
            "Fratanduono explicitly uses the Kraus Altshuler gamma law only for rho0/rho >= 0.64.",
            "Above that threshold, a candidate gamma is inferred from the vector-digitized Figure S2 Our Hugoniot Fit and the tabulated isentrope using Kraus Eq. (11). Whether the plotted Hugoniot exactly matches the corrected internal input is unconfirmed.",
            "No Hugoniot or gamma extrapolation is used. The figure ends at about 22.55 g/cm3, leaving the remaining isentrope rows unreconstructed.",
            "Ambient Table S1 density 8.938 is preserved, but excluded from fitting. Integration is anchored at the Table I reference 8.939 with zero pressure.",
            "Figure S2 blue-curve legend says isotherm, but its caption and numerical agreement with Table S1 identify isentrope; it is not independent 298 K fit input.",
            "Source Table S1 errors describe the isentrope and are not assigned to thermally reconstructed pressures. Diagnostic 298 K fits are unweighted; no source weights or covariance are inferred.",
            "Neither parameter parity nor the original full 298 K fit is reproduced. Published catalog coefficients remain unchanged.",
        ],
        "figure_blue_check": {
            "rmse_vs_table_s1_gpa": float(
                np.sqrt(np.mean((bpressure - raw[:, 0]) ** 2))
            ),
            "negative_implied_thermal_correction_rows": int(
                np.sum(observed_correction < 0)
            ),
            "minimum_table_minus_published_isotherm_gpa": float(
                observed_correction.min()
            ),
            "interpretation": "Separate rounded fitted curves cannot be subtracted to recover the original thermal correction; some table states fall below the published isotherm fit.",
        },
        "fits": {
            "source_isentrope_unweighted": fit_diagnostic(
                raw[:, 2], raw[:, 0], PUBLISHED_S
            ),
            "source_isentrope_reported_error_weighted": fit_diagnostic(
                raw[:, 2], raw[:, 0], PUBLISHED_S, raw[:, 1]
            ),
            "kraus_branch_298k_diagnostic": fit_diagnostic(
                r["density"][low], r["isotherm"][low], PUBLISHED_T
            ),
            "candidate_combined_298k_diagnostic": fit_diagnostic(
                r["density"], r["isotherm"], PUBLISHED_T
            ),
        },
        "curve_checks": {
            "published_isentrope_rmse_vs_table_gpa": float(
                np.sqrt(np.mean((ps.pressure(RHO0 / raw[:, 2]) - raw[:, 0]) ** 2))
            ),
            "kraus_branch_298k_rmse_vs_published_gpa": float(
                np.sqrt(np.mean((r["isotherm"][low] - pt.pressure(x[low])) ** 2))
            ),
            "candidate_combined_298k_rmse_vs_published_gpa": float(
                np.sqrt(np.mean((r["isotherm"] - pt.pressure(x)) ** 2))
            ),
        },
    }
    for name in ["kraus_branch_298k_diagnostic", "candidate_combined_298k_diagnostic"]:
        values = np.array(list(result["fits"][name]["parameters"].values()))
        result["fits"][name]["within_printed_parameter_errors"] = dict(
            zip(
                ["K0", "eta", "beta", "psi"],
                map(bool, abs(values - PUBLISHED_T) <= PUBLISHED_T_ERRORS),
            )
        )
    rows = list(
        zip(
            r["density"],
            r["isentrope"],
            r["gamma"],
            TR * r["factor"],
            THETA0 * r["factor"],
            r["correction"],
            r["isotherm"],
            np.where(low, "kraus_inherited", "candidate_digitized_hugoniot"),
        )
    )
    return result, rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--extract-supplement", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.extract_supplement:
        extract_supplement(args.extract_supplement)
    result, rows = reproduce()
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    target = DATA / f"{PREFIX}-thermal-reduction.json"
    if args.check:
        if target.read_text() != text:
            raise SystemExit("Cu reduction audit is stale")
    else:
        target.write_text(text)
        _write_csv(
            DATA / f"{PREFIX}-reconstructed-298k.csv",
            [
                "density_g_cm3",
                "source_isentrope_pressure_gpa",
                "gamma",
                "isentrope_temperature_k",
                "debye_temperature_k",
                "thermal_pressure_correction_gpa",
                "candidate_298k_pressure_gpa",
                "qualification",
            ],
            rows,
        )
    print(json.dumps(result["selection"]))


if __name__ == "__main__":
    main()
