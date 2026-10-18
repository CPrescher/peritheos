#!/usr/bin/env python3
"""Audit the staged Fe2SiO4 ringwoodite candidate without changing the catalog.

SciPy quadrature and explicit SI conversions provide an implementation
independent of Peritheos. Conditional fits use recovered subsets, never claim
recovery of the original combined-data fit, and do not overwrite coefficients.
"""

from __future__ import annotations

import argparse
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares
from scipy.special import bernoulli

from peritheos import Material, load_eosmat
from scripts.ferrous_ringwoodite_recover import DATA, PREFIX, read_csv

ROOT = Path(__file__).resolve().parents[1]
STAGE = ROOT / "output/material-candidates/ferrous-ringwoodite"
CARD = STAGE / "ferrous_ringwoodite.eosmat"
V0_MOLAR = 42.03
Z = 8
N = 7
AVOGADRO = 6.02214076e23
R = 8.31446261815324
CELL_TO_CM3_MOL = AVOGADRO * 1e-24 / Z
V0_CELL = V0_MOLAR / CELL_TO_CM3_MOL
THERMAL_ID = "ferrous_ringwoodite_armentrout_kavner_2011_bm2_mgd_candidate"
STATIC_ID = "ferrous_ringwoodite_finger_1979_bm2_ambient_cell_reconstruction"


def bm2(volume_cm3_mol: float, k0: float = 202.0, *, v0: float = V0_MOLAR) -> float:
    eta = (v0 / volume_cm3_mol) ** (1.0 / 3.0)
    return 1.5 * k0 * (eta**7 - eta**5)


def theta(volume_cm3_mol: float, gamma0: float, q: float, law: str) -> float:
    x = volume_cm3_mol / V0_MOLAR
    gamma = gamma0 * x**q
    if law == "variable_exponent":
        return 685.0 * x ** (-gamma)
    return 685.0 * np.exp(gamma0 / q * (1 - x**q))


def debye_energy(temperature: float, theta_k: float) -> float:
    """J/mol of Fe2SiO4 formula units; zero-point energy cancels at fixed V."""
    y = theta_k / temperature
    integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11)[0]
    return 9 * N * R * temperature * integral / y**3


def debye_energy_series(temperature: float, theta_k: float, terms: int) -> float:
    """Author dissertation Eq 4, truncated after B_(2*terms), J/mol.

    Valid for abs(theta/T) < 2*pi. The zero-point term is omitted, as in
    quadrature; it cancels in the same-volume reference-temperature difference.
    The original author truncation order is unknown.
    """
    y = theta_k / temperature
    if abs(y) >= 2 * np.pi:
        raise ValueError("Bernoulli Debye expansion requires abs(theta/T) < 2*pi")
    numbers = bernoulli(2 * terms)
    d = -3 * y / 8 + 3 * sum(
        numbers[2 * k] * y ** (2 * k) / ((2 * k + 3) * factorial(2 * k))
        for k in range(terms + 1)
    )
    return float(3 * N * R * temperature * d)


def thermal_equation_checks() -> dict:
    """Check thermodynamic identities and the author's independent Debye series."""
    identities = []
    for x in (0.8, 0.9, 1.0, 1.05):
        v = V0_MOLAR * x
        h = 1e-5
        expected = 1.08 * x**2
        implied = {}
        for law in ("integrated_gruneisen", "variable_exponent"):
            implied[law] = float(
                -(
                    np.log(theta(v * np.exp(h), 1.08, 2, law))
                    - np.log(theta(v * np.exp(-h), 1.08, 2, law))
                )
                / (2 * h)
            )
        identities.append({"volume_ratio": x, "gamma": expected, **implied})
    groups = observation_groups()
    series_checks = []
    for terms in (2, 4, 6, 10, 20):
        differences = []
        for group in groups.values():
            for v, t in zip(group["volume"], group["temperature"]):
                th = theta(v, 1.08, 2, "integrated_gruneisen")
                delta = (
                    debye_energy_series(t, th, terms)
                    - debye_energy_series(300, th, terms)
                    - debye_energy(t, th)
                    + debye_energy(300, th)
                )
                differences.append(abs(1.08 * (v / V0_MOLAR) ** 2 / v * delta / 1000))
        series_checks.append(
            {
                "highest_Bernoulli_index": 2 * terms,
                "maximum_pressure_difference_gpa": float(max(differences)),
            }
        )
    by_measurement = {
        r["measurement"]: r
        for r in read_csv(DATA / f"{PREFIX}armentrout-2011-table-s1.csv")
    }
    residual_checks = []
    for point in read_csv(DATA / f"{PREFIX}armentrout-2012-figure9-checkpoints.csv"):
        row = by_measurement[point["measurement"]]
        calculated = float(row["pressure_gpa"]) - pressure(
            float(row["derived_sample_molar_volume_cm3_mol"]),
            float(row["assumed_mean_sample_temperature_k"]),
        )
        digitized = float(point["digitized_residual_gpa"])
        tolerance = float(point["residual_tolerance_gpa"])
        residual_checks.append(
            {
                "measurement": point["measurement"],
                "measured_minus_model_gpa": calculated,
                "digitized_source_residual_gpa": digitized,
                "difference_gpa": calculated - digitized,
                "tolerance_gpa": tolerance,
                "passes": bool(abs(calculated - digitized) <= tolerance),
            }
        )
    return {
        "primary_formalism": "Dorogokupets-Dewaele 2007 Eqs 14-17: gamma_infinity=0, beta=q gives the integrated theta law; author dissertation Os Eqs 1-5 explicitly uses this limit and a 300 K reference.",
        "thermodynamic_identity": "gamma = -d ln(theta) / d ln(V)",
        "identity_checks": identities,
        "variable_exponent_implied_gamma": "gamma0*x^q*(1+q*ln(x)); inconsistent with gamma0*x^q away from x=1 when q is nonzero",
        "printed_integral_issue": "DD2007 Eq5 and dissertation prose on p10 print the inverse Debye prefactor. The normalized Debye function requires 3*(T/theta)^3. The dissertation Bernoulli Eq4 independently has the correct high-T limit and agrees with quadrature. No author code or truncation order recovered.",
        "series_checks": series_checks,
        "figure9_residual_checks": residual_checks,
        "figure9_scope": "Two isolated digitized published residual symbols, not recovered observation rows and never used in fits. Same systematic hot-data mismatch is present in the author's own figure; agreement does not establish the full author workflow.",
    }


def pressure(
    volume_cm3_mol: float,
    temperature: float,
    *,
    gamma0: float = 1.08,
    q: float = 2.0,
    law: str = "integrated_gruneisen",
    k0: float = 202.0,
) -> float:
    gamma = gamma0 * (volume_cm3_mol / V0_MOLAR) ** q
    th = theta(volume_cm3_mol, gamma0, q, law)
    # J / cm^3 = 10^6 Pa = 10^-3 GPa.
    thermal = (
        gamma
        / volume_cm3_mol
        * (debye_energy(temperature, th) - debye_energy(300.0, th))
        / 1000
    )
    return float(bm2(volume_cm3_mol, k0) + thermal)


def volume_ratio(pressure_gpa: float, temperature_k: float, **configuration) -> float:
    return float(
        brentq(
            lambda v: pressure(v, temperature_k, **configuration) - pressure_gpa,
            V0_MOLAR * 0.8,
            V0_MOLAR * 1.2,
            xtol=1e-11,
        )
        / V0_MOLAR
    )


def observation_groups() -> dict[str, dict]:
    new = read_csv(DATA / f"{PREFIX}armentrout-2011-table-s1.csv")
    accepted = [r for r in new if r["included_by_table_s1"] == "1"]
    liu = read_csv(DATA / f"{PREFIX}liu-2008-table1.csv")
    return {
        "armentrout_2011_accepted_14": {
            "volume": [
                float(r["derived_sample_molar_volume_cm3_mol"]) for r in accepted
            ],
            "temperature": [
                float(r["assumed_mean_sample_temperature_k"]) for r in accepted
            ],
            "pressure": [float(r["pressure_gpa"]) for r in accepted],
            "pressure_uncertainty": [
                float(r["pressure_reported_uncertainty_gpa"]) for r in accepted
            ],
        },
        "liu_2008_all_29": {
            # M=203.77 g/mol is a declared conversion choice, not a fitted quantity.
            # Three-place printed densities limit recovery of unrounded XRD volumes.
            "volume": [203.77 / float(r["density_g_cm3"]) for r in liu],
            "temperature": [float(r["temperature_k"]) for r in liu],
            "pressure": [float(r["pressure_gpa"]) for r in liu],
        },
    }


def residuals(group: dict, gamma0: float = 1.08, **configuration) -> np.ndarray:
    return np.array(
        [
            pressure(v, t, gamma0=gamma0, **configuration) - p
            for v, t, p in zip(group["volume"], group["temperature"], group["pressure"])
        ]
    )


def fit_gamma(group: dict, *, weighted: bool = False) -> dict:
    scale = np.array(group["pressure_uncertainty"]) if weighted else 1.0
    fit = least_squares(
        lambda x: residuals(group, float(x[0])) / scale,
        [1.08],
        bounds=(0.1, 3.0),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    source_residual = residuals(group)
    fitted_residual = residuals(group, float(fit.x[0]))
    return {
        "rows": len(group["pressure"]),
        "gamma0": float(fit.x[0]),
        "difference_from_published_gamma0": float(fit.x[0] - 1.08),
        "within_published_0_06_error_width": bool(abs(fit.x[0] - 1.08) <= 0.06),
        "published_pressure_rmse_gpa": float(np.sqrt(np.mean(source_residual**2))),
        "fitted_pressure_rmse_gpa": float(np.sqrt(np.mean(fitted_residual**2))),
        "maximum_published_absolute_pressure_residual_gpa": float(
            np.max(abs(source_residual))
        ),
        "objective": "sum((Pmodel-Psource)/reported_pressure_error)^2"
        if weighted
        else "sum((Pmodel-Psource)^2)",
        "fixed": {
            "V0_cm3_mol": 42.03,
            "K0_gpa": 202,
            "K0_prime": 4,
            "theta0_k": 685,
            "q": 2,
            "n": 7,
            "Tr_k": 300,
        },
        "scope": "Conditional recovered-subset refit; not the original combined fit. Errors/covariance are not reconstructed.",
    }


def reproduce() -> dict:
    groups = observation_groups()
    checks = read_csv(DATA / f"{PREFIX}armentrout-2011-figure1-checkpoints.csv")
    curve_checks = []
    for row in checks:
        p, t = float(row["pressure_gpa"]), float(row["temperature_k"])
        observed = float(row["digitized_volume_ratio"])
        predicted = volume_ratio(p, t)
        curve_checks.append(
            {
                "pressure_gpa": p,
                "temperature_k": t,
                "digitized_volume_ratio": observed,
                "calculated_volume_ratio": predicted,
                "difference": predicted - observed,
                "tolerance": float(row["volume_ratio_tolerance"]),
                "passes": bool(
                    abs(predicted - observed) <= float(row["volume_ratio_tolerance"])
                ),
                "variable_exponent_volume_ratio": volume_ratio(
                    p, t, law="variable_exponent"
                ),
            }
        )

    combined = {
        key: groups["armentrout_2011_accepted_14"][key] + groups["liu_2008_all_29"][key]
        for key in ("volume", "temperature", "pressure")
    }
    conditional = {name: fit_gamma(group) for name, group in groups.items()}
    conditional["armentrout_14_pressure_error_weighted"] = fit_gamma(
        groups["armentrout_2011_accepted_14"], weighted=True
    )
    conditional["armentrout_14_plus_liu_29"] = fit_gamma(combined)

    # Diagnostic addition of ALL seven positive-temperature Suzuki calculated
    # rows >=100 C, including extrapolations. Rebase its 20 C expansion to 300 K
    # by a declared interpolation of the 20/100 C table entries. This is not an
    # assertion about which Suzuki values the 2011 authors used.
    suzuki = read_csv(DATA / f"{PREFIX}suzuki-1979-table2-calculated.csv")
    y300 = np.interp(26.85, [20.0, 100.0], [0.0, 0.179]) / 100
    with_suzuki = {key: list(values) for key, values in combined.items()}
    for row in suzuki:
        if float(row["temperature_c"]) < 100:
            continue
        with_suzuki["volume"].append(
            V0_MOLAR * (1 + float(row["volume_expansion_percent"]) / 100) / (1 + y300)
        )
        with_suzuki["temperature"].append(float(row["temperature_c"]) + 273.15)
        with_suzuki["pressure"].append(0.0)
    conditional["armentrout_14_liu_29_suzuki_7_including_extrapolations"] = fit_gamma(
        with_suzuki
    )

    all_new = read_csv(DATA / f"{PREFIX}armentrout-2011-table-s1.csv")
    figure_temperatures = read_csv(
        DATA / f"{PREFIX}armentrout-2011-figure-s1-temperatures.csv"
    )
    by_measurement = {r["measurement"]: r for r in all_new}
    temperature_checks = []
    for row in figure_temperatures:
        observation = by_measurement[row["measurement"]]
        mean_temperature = float(observation["assumed_mean_sample_temperature_k"])
        source_temperature = int(row["published_figure_s1_temperature_k"])
        temperature_checks.append(
            {
                "measurement": row["measurement"],
                "table_s1_arithmetic_mean_k": mean_temperature,
                "figure_s1_temperature_k": source_temperature,
                "rounded_mean_matches": bool(
                    int(np.floor(mean_temperature + 0.5)) == source_temperature
                ),
            }
        )
    raw_density = np.array([float(r["source_density_column"]) for r in all_new])
    physical_volume = np.array(
        [float(r["derived_sample_molar_volume_cm3_mol"]) for r in all_new]
    )
    # Numerical clue only; do not reinterpret a mislabeled source column silently.
    inverse_product = raw_density * physical_volume
    checkpoints = []
    # Explicit opt-in because exact source equation mapping remains qualified.
    material = Material.from_eosmat(load_eosmat(CARD), require_primary_validation=False)
    record = material.get_eos_record(THERMAL_ID)
    for p, t in [(0, 300), (5, 1000), (10, 1500), (24, 2000), (30, 2000)]:
        ratio = volume_ratio(p, t)
        cell = ratio * V0_CELL
        checkpoints.append(
            {
                "pressure_gpa": p,
                "temperature_k": t,
                "volume_a3": cell,
                "native_volume_a3": float(record.volume(p, t)),
                "independent_pressure_gpa": pressure(cell * CELL_TO_CM3_MOL, t),
                "native_pressure_gpa": float(record.pressure(cell, t)),
            }
        )

    static_record = material.get_eos_record(STATIC_ID)
    static_rows = read_csv(DATA / f"{PREFIX}finger-1979-table2.csv")
    static_checks = []
    linear_factors = []
    observed_pressures = []
    for row in static_rows:
        observed_p = float(row["pressure_kbar"]) / 10
        observed_v = float(row["volume_a3"])
        calculated_p = bm2(observed_v, 196.0, v0=558.7)
        calculated_v = brentq(
            lambda v: bm2(v, 196.0, v0=558.7) - observed_p,
            500.0,
            560.0,
            xtol=1e-11,
        )
        if observed_p > 0.1:
            linear_factors.append(calculated_p / 196.0)
            observed_pressures.append(observed_p)
        eta = (558.7 / observed_v) ** (1 / 3)
        bulk = 196.0 / 2 * (7 * eta**7 - 5 * eta**5)
        volume_pressure_error = bulk / observed_v * float(row["volume_esd_a3"])
        pressure_error = float(row["pressure_esd_kbar"] or 0) / 10
        # Only measured P and V uncertainties, not fitted K0/V0 errors or
        # invented independence/covariance of all fit parameters.
        combined_coordinate_error = np.hypot(volume_pressure_error, pressure_error)
        static_checks.append(
            {
                "observed_pressure_gpa": observed_p,
                "observed_volume_a3": observed_v,
                "calculated_pressure_gpa": calculated_p,
                "calculated_volume_a3": calculated_v,
                "native_pressure_gpa": float(
                    static_record.pressure(observed_v, 296.15)
                ),
                "native_volume_a3": float(static_record.volume(observed_p, 296.15)),
                "pressure_difference_gpa": calculated_p - observed_p,
                "combined_coordinate_error_gpa": float(combined_coordinate_error),
                "within_combined_coordinate_error": bool(
                    abs(calculated_p - observed_p) <= combined_coordinate_error
                ),
            }
        )
    factors = np.array(linear_factors)
    conditional_k0 = float(factors @ np.array(observed_pressures) / (factors @ factors))

    return {
        "decision": "staged source-validated static ambient-cell reconstruction; separate 2011 thermal mapping remains qualified, without complete author fit",
        "static_1979_reconstruction": {
            "record_kind": "derived",
            "temperature_ref_k": 296.15,
            "published_K0_gpa": 196.0,
            "published_K0_error_gpa": 6.0,
            "adopted_measured_V0_a3": 558.7,
            "adopted_V0_measurement_esd_a3": 0.2,
            "source_K0_prime_fixed": 4.0,
            "conditional_unweighted_compressed_two_row_K0_gpa": conditional_k0,
            "author_fit_reproduced": False,
            "source_checks": static_checks,
        },
        "source_equation_assessment": {
            "published_coefficients_checked": True,
            "explicit_ringwoodite_equations_in_2011_article_or_supplement": False,
            "supporting_author_formalism": "Recovered and visually checked: Armentrout 2012 dissertation, preceding Os chapter Eqs 1-5. The cited Dorogokupets-Dewaele 2007 Eqs 14-17 independently imply the integrated theta law in the constant-q limit. Applying the Os implementation to ringwoodite remains an inference; exact author code is unavailable.",
            "figure1_all_checkpoints_pass": all(r["passes"] for r in curve_checks),
            "theta_law_uniquely_identified_by_figure": False,
            "exact_author_workflow_reproduced": False,
        },
        "thermal_equation_followup": thermal_equation_checks(),
        "units": {
            "Z": Z,
            "atoms_per_formula": N,
            "V0_cm3_mol": V0_MOLAR,
            "V0_cell_a3": V0_CELL,
            "cell_a3_to_cm3_mol": CELL_TO_CM3_MOL,
            "native_volume_is_one_tenth_cm3_mol": True,
        },
        "primary_data": {
            "new_rows": 18,
            "new_included_rows": 14,
            "liu_rows": 29,
            "suzuki_calculated_rows": 10,
            "finger_structure_compression_rows": 3,
            "excluded_measurements": [
                r["measurement"] for r in all_new if r["included_by_table_s1"] == "0"
            ],
            "accepted_new_pressure_range_gpa": [
                min(groups["armentrout_2011_accepted_14"]["pressure"]),
                max(groups["armentrout_2011_accepted_14"]["pressure"]),
            ],
            "accepted_new_assumed_mean_temperature_range_k": [
                min(groups["armentrout_2011_accepted_14"]["temperature"]),
                max(groups["armentrout_2011_accepted_14"]["temperature"]),
            ],
        },
        "inconsistent_printed_density_column": {
            "maximum_difference_from_molar_volume_cm3_mol": float(
                np.max(abs(raw_density - physical_volume))
            ),
            "product_with_lattice_derived_molar_volume_mean": float(
                np.mean(inverse_product)
            ),
            "product_standard_deviation": float(np.std(inverse_product)),
            "interpretation": "Neither molar volume nor Fe2SiO4 mass density; retained verbatim, never used for fitting.",
        },
        "source_static_tradeoff": {
            "K0_at_Kprime_4_gpa": 214.3 - 3.2 * 4,
            "stored_preferred_K0_gpa": 202,
            "difference_gpa": 202 - (214.3 - 3.2 * 4),
            "note": "Rounded 201.5 is compatible with preferred 202(4); do not replace the preferred value.",
        },
        "figure1_checkpoints": curve_checks,
        "figure_s1_temperature_checkpoints": temperature_checks,
        "conditional_gamma_refits": conditional,
        "independent_native_checkpoints": checkpoints,
        "unresolved": [
            "Full Sato 1977, Wilburn-Bassett 1976, Nestola 2010 and Plymate-Stout 1994 observation tables were not recovered.",
            "The 2011 article includes Suzuki 1979 in its thermal-data list; the dissertation ringwoodite chapter omits it. Suzuki Fe data are modeled/extrapolated from Mao 1969, not new Fe measurements.",
            "The 2011 ringwoodite article/supplements do not print the exact cold equation, gamma/theta equations, weighting, full input normalization, row mask, covariance or numerical Debye truncation. Figure S1 DOES independently confirm rounded two-surface mean temperatures for all nine cycle2 measurements.",
            "Brown1999 primary TablesI/II are recovered and all18 paired NaCl readings replayed in ferrous_ringwoodite_brown_replay.py. The two published representations differ by about0.12GPa locally; stated midpoint-temperature replay reduces but does not remove hot residuals. Exact original temperature/error reduction remains unresolved. No other NaCl EOS is substituted.",
            "2011 does not establish sample H content. Wang et al 2025 discuss possible hydration as an explanation for the low thermal expansion; this is a hypothesis, not a demonstrated correction.",
        ],
    }


def plot_residuals(path: Path) -> None:
    """Plot the fixed published coefficients; the displayed points are not a fit."""
    import matplotlib.pyplot as plt

    group = observation_groups()["armentrout_2011_accepted_14"]
    source = read_csv(DATA / f"{PREFIX}armentrout-2012-figure9-checkpoints.csv")
    temperatures = dict(
        zip(
            [
                r["measurement"]
                for r in read_csv(DATA / f"{PREFIX}armentrout-2011-table-s1.csv")
                if r["included_by_table_s1"] == "1"
            ],
            group["temperature"],
        )
    )
    with plt.rc_context(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    ):
        fig, ax = plt.subplots(figsize=(8.2, 4.7), layout="constrained")
        ax.axhline(0, color="#687684", linewidth=1, linestyle="--")
        ax.scatter(
            group["temperature"],
            -residuals(group),
            color="#174d78",
            s=42,
            label="Independent replay: 14 accepted Table S1 rows",
            zorder=3,
        )
        ax.errorbar(
            [temperatures[r["measurement"]] for r in source],
            [float(r["digitized_residual_gpa"]) for r in source],
            yerr=[float(r["residual_tolerance_gpa"]) for r in source],
            fmt="s",
            markerfacecolor="none",
            markeredgecolor="#c25d1c",
            ecolor="#c25d1c",
            markersize=9,
            capsize=3,
            label="Digitized dissertation Fig. 9; graphical tolerance",
            zorder=4,
        )
        ax.set(
            xlabel="Mean sample temperature (K)",
            ylabel="Measured pressure minus model pressure (GPa)",
            title="Fe₂SiO₄ ringwoodite: the published hot-data residual trend",
        )
        ax.set_ylim(-2.9, 0.85)
        ax.grid(axis="y", color="#dbe2e8", linewidth=0.6)
        ax.legend(loc="lower left", frameon=False, fontsize=9)
        fig.text(
            0.13,
            -0.045,
            "Published coefficients held fixed; integrated constant-q MGD. "
            "No refit. Source residual symbols are never regression inputs.",
            fontsize=8,
            color="#46535e",
        )
        fig.savefig(path, dpi=180, bbox_inches="tight")
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--residual-plot", type=Path)
    args = parser.parse_args()
    serialized = json.dumps(reproduce(), indent=2, allow_nan=False) + "\n"
    if args.output:
        args.output.write_text(serialized, encoding="utf-8")
    else:
        print(serialized, end="")
    if args.residual_plot:
        plot_residuals(args.residual_plot)


if __name__ == "__main__":
    main()
