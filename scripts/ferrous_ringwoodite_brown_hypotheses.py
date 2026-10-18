#!/usr/bin/env python3
"""Test temperature, atom-count and energy-reference hypotheses without corrections.

All objectives are explicitly unweighted pressure differences in GPa. Fits to
Brown's calculated grid characterize numerical signatures, not experiments or
physical validation. Neither original tables nor EOS coefficients are changed.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from math import factorial
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares
from scipy.special import bernoulli

from scripts.ferrous_ringwoodite_brown_replay import (
    MOLAR_MASS,
    TABLE2,
    Brown1999,
    specific_volume,
)
from scripts.ferrous_ringwoodite_recover import DATA, PREFIX, read_csv
from scripts.ferrous_ringwoodite_reproduce import STAGE, V0_MOLAR, R, bm2, pressure


def zero_point_pressure(model: Brown1999, volume: float) -> float:
    """Full Debye zero-point pressure for n=2, GPa; diagnostic only."""
    _, gamma, theta = model.parameters(volume)
    return gamma * (9 / 8 * 2 * R * theta) / (volume * MOLAR_MASS) / 1000


def series_energy(temperature: float, theta: float, *, linear: float = -3 / 8) -> float:
    """Independent Bernoulli series, n=2; vary only the linear coefficient."""
    y = theta / temperature
    if abs(y) >= 2 * np.pi:
        raise ValueError("Debye Bernoulli series requires abs(theta/T)<2*pi")
    numbers = bernoulli(16)
    d3 = (
        1
        + linear * y
        + 3
        * sum(
            numbers[2 * k] * y ** (2 * k) / ((2 * k + 3) * factorial(2 * k))
            for k in range(1, 9)
        )
    )
    return float(6 * R * temperature * d3)


def statistics(values) -> dict:
    a = np.asarray(values, dtype=float)
    return {
        "count": len(a),
        "mean_gpa": float(a.mean()),
        "rms_gpa": float(np.sqrt(np.mean(a**2))),
        "maximum_abs_gpa": float(max(abs(a))),
    }


def brown_table_diagnostics(model: Brown1999) -> dict:
    cells = []
    centered_e, centered_p = [], []
    series_standard_error, series_altered_error = [], []
    for row in read_csv(TABLE2):
        v = float(row["specific_volume_cm3_g"])
        p0, gamma, theta = model.parameters(v)
        volume_cells = []
        for t in range(300, 1201, 100):
            if not row[f"pressure_{t}k_gpa"]:
                continue
            target = float(row[f"pressure_{t}k_gpa"])
            thermal = model.thermal(v, t)
            zp = zero_point_pressure(model, v)
            cells.append([v, t, target - p0, thermal, zp])
            volume_cells.append([thermal, target])
            factor = gamma / (v * MOLAR_MASS) / 1000
            series_standard_error.append(factor * series_energy(t, theta) - thermal)
            series_altered_error.append(
                factor * series_energy(t, theta, linear=-1 / 4) - thermal - zp / 3
            )
        if len(volume_cells) >= 2:
            x, y = np.asarray(volume_cells).T
            centered_e.extend(x - x.mean())
            centered_p.extend(y - y.mean())
    a = np.asarray(cells)
    target, thermal, zp = a[:, 2], a[:, 3], a[:, 4]
    e, p = np.asarray(centered_e), np.asarray(centered_p)
    n_from_slopes = 2 * float(e @ p / (e @ e))
    n_only = 2 * float(thermal @ target / (thermal @ thermal))
    f_at_n2 = float(zp @ (target - thermal) / (zp @ zp))
    constant_pressure = float(np.mean(target - thermal))
    energy_factor = np.array(
        [model.parameters(v)[1] / (v * MOLAR_MASS) / 1000 for v in a[:, 0]]
    )
    constant_energy = float(
        energy_factor @ (target - thermal) / (energy_factor @ energy_factor)
    )
    design = np.column_stack([thermal, zp])
    scale, f = np.linalg.lstsq(design, target, rcond=None)[0]
    train = a[:, 1] <= 700
    cv = np.linalg.lstsq(design[train], target[train], rcond=None)[0]
    required_n = []
    v = 0.3709
    row = next(r for r in read_csv(TABLE2) if float(r["specific_volume_cm3_g"]) == v)
    for t in (300, 600, 1200):
        required_n.append(
            {
                "volume_cm3_g": v,
                "temperature_k": t,
                "required_n_if_no_offset": 2
                * (float(row[f"pressure_{t}k_gpa"]) - model.parameters(v)[0])
                / model.thermal(v, t),
            }
        )
    return {
        "scope": "All571 populated calculated pressure cells; equal weight per cell. Pressure differences are prediction minus TableII, GPa; no original experimental refit.",
        "wrong_n_only": {
            str(n): statistics(n / 2 * thermal - target) for n in (1, 2, 4, 7, 8)
        },
        "n_from_within_volume_temperature_slopes": n_from_slopes,
        "n_only_fitted_with_fixed_table1_p0": {
            "n": n_only,
            **statistics(n_only / 2 * thermal - target),
        },
        "required_n_by_temperature_at_one_volume": required_n,
        "fixed_one_third_zero_point": statistics(thermal + zp / 3 - target),
        "fixed_full_zero_point": statistics(thermal + zp - target),
        "fitted_constant_pressure_offset": {
            "offset_gpa": constant_pressure,
            **statistics(thermal + constant_pressure - target),
        },
        "fitted_constant_molar_energy_offset": {
            "offset_j_mol": constant_energy,
            **statistics(thermal + constant_energy * energy_factor - target),
        },
        "fitted_zero_point_fraction_at_n2": {
            "fraction": f_at_n2,
            **statistics(thermal + f_at_n2 * zp - target),
        },
        "joint_n_and_offset_diagnostic": {
            "n": float(2 * scale),
            "fraction_of_n2_zero_point": float(f),
            **statistics(scale * thermal + f * zp - target),
        },
        "temperature_holdout": {
            "training_temperature_max_k": 700,
            "test_temperature_min_k": 800,
            "fitted_n": float(2 * cv[0]),
            "fitted_fraction_of_n2_zero_point": float(cv[1]),
            "test": statistics(
                cv[0] * thermal[~train] + cv[1] * zp[~train] - target[~train]
            ),
        },
        "independent_series_checks": {
            "normalized_series_minus_quadrature_max_gpa": float(
                max(abs(np.asarray(series_standard_error)))
            ),
            "altered_linear_term_minus_added_zpe_third_max_gpa": float(
                max(abs(np.asarray(series_altered_error)))
            ),
        },
        "interpretation": "TableII temperature slopes support n=2. An added E_ZP/3, E_ZP=9*n*R*theta/8, reproduces the offset closely and is algebraically equivalent to changing D3's linear coefficient from-3/8 to-1/4. This is a numerical signature, not recovered Brown code or a justified correction. Eq1/TableI define a zero-K baseline, so adding zero-point pressure needs a consistent redefinition of that baseline.",
    }


def temperature_diagnostics(model: Brown1999) -> dict:
    original = read_csv(DATA / f"{PREFIX}armentrout-2011-table-s1.csv")
    results = []
    for row in original:
        name = row["measurement"]
        t = float(row["assumed_mean_sample_temperature_k"])
        tm = (t + 300) / 2
        heated = "heat" in name and "preheat" not in name
        half = t / 2 if heated else 300
        v = specific_volume(float(row["nacl_b1_a_angstrom"]))
        p0 = model.parameters(v)[0]

        def th(z):
            return model.thermal(v, z)

        p_source = float(row["pressure_gpa"])
        predicted = {
            "table1_n2_midpoint": p0 + th(tm),
            "table1_n2_half_kelvin_heating": p0 + th(half),
            "table1_n1_full_sample": p0 + th(t) / 2,
            "table1_n1_midpoint": p0 + th(tm) / 2,
            "table1_n4_midpoint": p0 + th(tm) * 2,
            "table1_n7_midpoint": p0 + th(tm) * 3.5,
            "table1_n8_midpoint": p0 + th(tm) * 4,
            "table1_n2_midpoint_minus_E300_global": p0 + th(tm) - th(300),
            "table1_n2_midpoint_minus_E300_heating_only": p0
            + th(tm)
            - (th(300) if heated else 0),
            "table1_n2_half_celsius_heating": p0
            + th(t / 2 + 273.15 / 2 if heated else 300),
            "table1_n2_half_kelvin_plus_ZPE_third": p0
            + th(half)
            + zero_point_pressure(model, v) / 3,
            "table2_n2_midpoint": model.table2(v, tm, extend=True),
            "table2_n2_half_kelvin_heating": model.table2(v, half, extend=True),
        }
        for hypothesis, p in predicted.items():
            results.append(
                {
                    "measurement": name,
                    "included_by_table_s1": int(row["included_by_table_s1"]),
                    "hypothesis": hypothesis,
                    "source_pressure_gpa": p_source,
                    "predicted_pressure_gpa": p,
                    "predicted_minus_source_gpa": p - p_source,
                }
            )
    summaries = {}
    for hypothesis in predicted:
        rows = [
            r
            for r in results
            if r["hypothesis"] == hypothesis and r["included_by_table_s1"]
        ]
        summaries[hypothesis] = statistics(
            [r["predicted_minus_source_gpa"] for r in rows]
        )
    accepted = [r for r in original if r["included_by_table_s1"] == "1"]
    volumes = [specific_volume(float(r["nacl_b1_a_angstrom"])) for r in accepted]
    t = np.array([float(r["assumed_mean_sample_temperature_k"]) for r in accepted])
    p_source = np.array([float(r["pressure_gpa"]) for r in accepted])
    p0 = np.array([model.parameters(v)[0] for v in volumes])

    def energy(temperatures):
        return np.array([model.thermal(v, t0) for v, t0 in zip(volumes, temperatures)])

    functions = {
        "n_at_stated_midpoint": (
            lambda x: p0 + x[0] / 2 * energy((t + 300) / 2),
            [2.0],
        ),
        "temperature_fraction_without_intercept": (
            lambda x: p0 + energy(x[0] * t),
            [0.5],
        ),
        "offset_with_fixed_half_slope_k": (lambda x: p0 + energy(t / 2 + x[0]), [0.0]),
        "temperature_fraction_and_intercept_k": (
            lambda x: p0 + energy(x[0] * t + x[1]),
            [0.5, 0.0],
        ),
    }
    fitted = {}
    for name, (function, initial) in functions.items():
        fit = least_squares(lambda x: function(x) - p_source, initial)
        predictions = function(fit.x)
        heldout = []
        for i in range(len(t)):
            mask = np.arange(len(t)) != i
            loo = least_squares(lambda x: (function(x) - p_source)[mask], initial)
            heldout.append(function(loo.x)[i] - p_source[i])
        fitted[name] = {
            "parameters": fit.x.tolist(),
            "training": statistics(predictions - p_source),
            "leave_one_row_out": statistics(heldout),
        }
    return {
        "scope": "All18 rows retained; all objectives/summaries use14 accepted rows only, unweighted GPa residuals, no fitted EOS coefficients. No temperature uncertainty/covariance or significance inferred.",
        "hypothesis_definitions": "Names with half_kelvin use Tsample/2 for heating only; cold preheat/quench states stay300K. n1_full_sample uses one atom per formula with the formula molar volume, intentionally a mismatched basis. minus_E300 subtracts a300K thermal energy from an already zero-K P0 and is deliberately an error hypothesis. TableII reconstructions keep its offset; TableI cases keep its printed baseline. Celsius halving means(Tsample-273.15)/2+273.15. ZPE/3 cases are diagnostics, not corrections.",
        "accepted_statistics": summaries,
        "conditional_unweighted_fits": fitted,
        "rows": results,
        "interpretation": "Half-Kelvin n2 is favored among tested fixed-temperature/n alternatives, but the extra300K-reference subtraction is comparably close. The pressure table does not uniquely identify the historical workflow. A fitted noninteger n is a diagnostic scaling parameter, not a revised atom count. Cold excluded rows disfavor applying extra300K subtraction globally, but cannot rule out a separate heating routine. No source code recovered.",
    }


def ringwoodite_n_curve_check() -> dict:
    """Hold published parameters fixed and check n against calculated Figure1 curves."""
    thermal_points = [
        row
        for row in read_csv(DATA / f"{PREFIX}armentrout-2011-figure1-checkpoints.csv")
        if float(row["temperature_k"]) > 300
    ]
    result = {}
    for n in (1, 2, 4, 7):
        points = []
        for row in thermal_points:
            t, p = float(row["temperature_k"]), float(row["pressure_gpa"])
            root = (
                brentq(
                    lambda v: bm2(v) + n / 7 * (pressure(v, t) - bm2(v)) - p,
                    V0_MOLAR * 0.8,
                    V0_MOLAR * 1.2,
                )
                / V0_MOLAR
            )
            points.append(root - float(row["digitized_volume_ratio"]))
        result[str(n)] = {
            "thermal_checkpoints": len(points),
            "within_volume_ratio_tolerance_0_002": sum(abs(d) <= 0.002 for d in points),
            "maximum_absolute_volume_ratio_difference": max(abs(d) for d in points),
        }
    return {
        "scope": "Nine digitized calculated thermal curves, not experimental observations. Published gamma0=1.08,q=2,theta0=685K and static parameters held fixed. Conditional check under integrated theta reconstruction, not recovered author code or an independent physical measurement of n.",
        "by_n": result,
    }


def investigate() -> dict:
    model = Brown1999()
    return {
        "audit_date": "2026-10-16",
        "brown_source_doi": "10.1063/1.371596",
        "normalization": "Per mole of NaCl: n=2 and V_m=v*58.44cm3/mol. Per mole of conventional B1 cells: n=8 and volume4*V_m gives exactly the same pressure. Wrong n means an inconsistent atom/volume basis, not merely use of8. Fe2SiO4 similarly has n=7 per formula or56 per conventional cell with8times the formula molar volume.",
        "brown_tables": brown_table_diagnostics(model),
        "table_s1": temperature_diagnostics(model),
        "ringwoodite_n_curve_check": ringwoodite_n_curve_check(),
        "source_search": "Brown's public institutional Publications/Papers lists, four tool/example ZIPs plusExamples.zip, and public jmichaelb repository trees were checked. No1999 NaCl Debye calculator found; author2020 supplement only corroborates literalTableI P0/V/K/Kprime. No author contacted.",
        "historical_temperature_convention_confirmed": False,
        "brown_offset_cause_confirmed": False,
        "source_pressures_or_eos_coefficients_changed": False,
    }


def plot(result: dict, path: Path) -> None:
    import matplotlib.pyplot as plt

    model = Brown1999()
    row = next(r for r in read_csv(TABLE2) if r["specific_volume_cm3_g"] == "0.3709")
    ts = np.arange(300, 1201, 100)
    observed = np.array([float(row[f"pressure_{t}k_gpa"]) for t in ts])
    with plt.rc_context(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    ):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.8), layout="constrained")
        for n in (1, 2, 4):
            ps = np.array(
                [
                    model.parameters(0.3709)[0] + n / 2 * model.thermal(0.3709, t)
                    for t in ts
                ]
            )
            axes[0].plot(ts, ps - observed, marker="o", label=f"n={n}")
        ps = np.array(
            [
                model.table1(0.3709, t) + zero_point_pressure(model, 0.3709) / 3
                for t in ts
            ]
        )
        axes[0].plot(
            ts,
            ps - observed,
            marker="s",
            label="n=2, add ZPE/3 diagnostic",
            c="#247a67",
        )
        axes[0].axhline(0, c="grey", ls="--", lw=1)
        axes[0].set(
            xlabel="Temperature (K)",
            ylabel="Calculated - Table II pressure (GPa)",
            title="Brown at v=0.3709 cm3/g\nAtom count versus energy offset",
        )
        axes[0].legend(fontsize=8)
        keys = [
            "table1_n2_midpoint",
            "table1_n2_half_kelvin_heating",
            "table1_n1_full_sample",
            "table1_n2_midpoint_minus_E300_global",
            "table2_n2_half_kelvin_heating",
        ]
        labels = [
            "n=2, stated midpoint",
            "n=2, half sample Kelvin",
            "n=1, full sample Kelvin",
            "n=2, midpoint, subtract E(300)",
            "Table II, half sample Kelvin",
        ]
        rms = [result["table_s1"]["accepted_statistics"][k]["rms_gpa"] for k in keys]
        axes[1].barh(
            labels, rms, color=["#8a6f5e", "#247a67", "#637fa2", "#a56821", "#5f9b8f"]
        )
        axes[1].invert_yaxis()
        axes[1].set(
            xlabel="Unweighted pressure replay RMS (GPa)",
            title="Table S1 pressure replay\n14 accepted rows, fixed hypotheses",
        )
        fig.savefig(path, dpi=180)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=STAGE / f"{PREFIX}brown-hypotheses.json"
    )
    parser.add_argument(
        "--rows", type=Path, default=STAGE / f"{PREFIX}brown-hypotheses.csv"
    )
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    result = investigate()
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    stream = io.StringIO()
    rows = result["table_s1"]["rows"]
    writer = csv.DictWriter(stream, fieldnames=list(rows[0]), lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)
    args.rows.write_text(stream.getvalue())
    if args.plot:
        plot(result, args.plot)
    print(
        json.dumps(
            {
                "brown_n_from_slopes": result["brown_tables"][
                    "n_from_within_volume_temperature_slopes"
                ],
                "brown_zpe_third": result["brown_tables"]["fixed_one_third_zero_point"],
                "table_s1": result["table_s1"]["accepted_statistics"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
