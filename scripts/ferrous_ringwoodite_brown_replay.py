#!/usr/bin/env python3
"""Replay the Fe-ringwoodite NaCl readings using Brown's published tables.

Both published representations are retained because Table I + the normalized
Debye energy does not exactly reproduce Table II. No source pressure or EOS
coefficient is overwritten. This is a calibration diagnostic, not a refit.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import CubicSpline
from scipy.optimize import brentq

from scripts.ferrous_ringwoodite_recover import DATA, PREFIX, read_csv
from scripts.ferrous_ringwoodite_reproduce import AVOGADRO, STAGE, R, pressure

TABLE1 = DATA / f"{PREFIX}brown-1999-table1.csv"
TABLE2 = DATA / f"{PREFIX}brown-1999-table2.csv"
MOLAR_MASS = 58.44  # g/mol of NaCl; explicit conversion choice, not a fit.
TEMPERATURES = np.arange(300, 1201, 100)


def specific_volume(a_angstrom: float, molar_mass: float = MOLAR_MASS) -> float:
    """B1 conventional cell: four NaCl formula units; result cm^3/g."""
    return a_angstrom**3 * 1e-24 * AVOGADRO / (4 * molar_mass)


def vibrational_energy(temperature: float, theta: float) -> float:
    """Thermal Debye energy, J/mol NaCl (two atoms); excludes zero point."""
    if temperature <= 0:
        raise ValueError("Temperature must be positive")
    y = theta / temperature
    integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11)[0]
    return float(18 * R * temperature * integral / y**3)


class Brown1999:
    def __init__(self):
        rows = read_csv(TABLE1)
        self.volume = np.array([float(r["specific_volume_cm3_g"]) for r in rows])
        # MATLAB SPLINE's default not-a-knot boundary convention.
        self.splines = {
            key: CubicSpline(self.volume, [float(r[key]) for r in rows])
            for key in ("static_pressure_gpa", "gamma", "debye_temperature_k")
        }
        # Expanded-volume rows have leading blank cells: retain them in the
        # source CSV, but use only the complete rectangular part here.
        complete = [r for r in read_csv(TABLE2) if r["pressure_300k_gpa"]]
        self.table2_volume = np.array(
            [float(r["specific_volume_cm3_g"]) for r in complete]
        )
        self.table2_spline = CubicSpline(
            self.table2_volume,
            [[float(r[f"pressure_{t}k_gpa"]) for t in TEMPERATURES] for r in complete],
            axis=0,
        )

    def parameters(self, volume: float) -> tuple[float, float, float]:
        if not self.volume[0] <= volume <= self.volume[-1]:
            raise ValueError("Specific volume outside Brown Table I")
        return tuple(float(f(volume)) for f in self.splines.values())

    def thermal(self, volume: float, temperature: float) -> float:
        _, gamma, theta = self.parameters(volume)
        # J/cm^3 -> 10^-3 GPa; volume*M is cm^3/mol.
        return (
            gamma
            * vibrational_energy(temperature, theta)
            / (volume * MOLAR_MASS)
            / 1000
        )

    def table1(self, volume: float, temperature: float) -> float:
        return self.parameters(volume)[0] + self.thermal(volume, temperature)

    def table2(
        self, volume: float, temperature: float, *, extend: bool = False
    ) -> float:
        if not self.table2_volume[0] <= volume <= self.table2_volume[-1]:
            raise ValueError("Specific volume outside complete Brown Table II grid")
        boundary = float(np.clip(temperature, 300, 1200))
        if temperature != boundary and not extend:
            raise ValueError(
                "Temperature outside Brown Table II; explicit extension required"
            )
        interpolated = float(
            CubicSpline(TEMPERATURES, self.table2_spline(volume))(boundary)
        )
        # Extension is continuous at the nearest table boundary. It uses the
        # published gamma/theta, not an unconstrained cubic extrapolation in T.
        return (
            interpolated
            + self.thermal(volume, temperature)
            - self.thermal(volume, boundary)
        )


def source_pdf_check(source_dir: Path) -> dict:
    """Compare every Table I cell and every populated Table II cell to PDF text."""
    result = subprocess.run(
        ["pdftotext", "-layout", str(source_dir / "brown-1999.pdf"), "-"],
        capture_output=True,
        check=True,
        text=True,
    )
    text = result.stdout.replace("⫺", "-")
    first = text.split("TABLE I. Properties")[1].split("model 关Eqs.")[0]
    pattern = (
        r"\b(0\.\d{4})\s+(-?\d+\.\d{2})\s+(\d+\.\d{2})"
        r"\s+(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d{3})\b"
    )
    table1 = re.findall(pattern, first)
    second = text.split("TABLE II. Equation")[1].split("\f")[0]
    table2 = []
    for line in second.splitlines():
        fields = line.split()
        if (
            len(fields) >= 3
            and re.fullmatch(r"0\.\d{4}", fields[0])
            and re.fullmatch(r"-?0\.\d{4}", fields[1])
        ):
            table2.append(fields[:2] + [""] * (12 - len(fields)) + fields[2:])
    first_match = [list(r) for r in table1] == [
        list(r.values()) for r in read_csv(TABLE1)
    ]
    second_match = table2 == [list(r.values()) for r in read_csv(TABLE2)]
    if len(table1) != 21 or len(table2) != 63 or not first_match or not second_match:
        raise ValueError("Brown primary table transcription mismatch")
    return {
        "table1_rows": len(table1),
        "table1_numeric_cells": len(table1) * 6,
        "table2_rows": len(table2),
        "table2_populated_pressure_cells": sum(
            bool(v) for row in table2 for v in row[2:]
        ),
        "all_source_cells_and_blanks_match": first_match and second_match,
    }


def replay() -> dict:
    model = Brown1999()
    checks = []
    thermal_differences = []
    for r in read_csv(TABLE2):
        v = float(r["specific_volume_cm3_g"])
        for t in TEMPERATURES:
            entry = r[f"pressure_{t}k_gpa"]
            if not entry:
                continue
            difference = model.table1(v, float(t)) - float(entry)
            checks.append(
                {
                    "specific_volume_cm3_g": v,
                    "temperature_k": int(t),
                    "table1_minus_table2_gpa": difference,
                }
            )
            if r["pressure_300k_gpa"]:
                delta = model.thermal(v, float(t)) - model.thermal(v, 300)
                thermal_differences.append(
                    delta - (float(entry) - float(r["pressure_300k_gpa"]))
                )
    rows = []
    for r in read_csv(DATA / f"{PREFIX}armentrout-2011-table-s1.csv"):
        a = float(r["nacl_b1_a_angstrom"])
        sd_a = float(r["nacl_b1_a_sd_angstrom"])
        v = specific_volume(a)
        t = float(r["assumed_mean_sample_temperature_k"])
        tm = float(r["assumed_nacl_temperature_k"])
        source = float(r["pressure_gpa"])
        ring = pressure(float(r["derived_sample_molar_volume_cm3_mol"]), t)
        p1 = model.table1(v, tm)
        p2 = model.table2(v, tm, extend=True)
        h = 1e-5
        dp_da = (
            model.table2(specific_volume(a + h), tm, extend=True)
            - model.table2(specific_volume(a - h), tm, extend=True)
        ) / (2 * h)
        dp_dt = (
            model.table2(v, tm + 0.01, extend=True)
            - model.table2(v, tm - 0.01, extend=True)
        ) / 0.02
        required_source_t1 = brentq(lambda z: model.table1(v, z) - source, 1, 4000)
        required_source_t2 = brentq(
            lambda z: model.table2(v, z, extend=True) - source, 1, 4000
        )
        required_ring_t = (
            brentq(lambda z: model.table2(v, z, extend=True) - ring, 1, 4000)
            if model.table2(v, 1, extend=True)
            <= ring
            <= model.table2(v, 4000, extend=True)
            else None
        )
        half_t = (
            t / 2
            if "heat" in r["measurement"] and "preheat" not in r["measurement"]
            else 300
        )
        rows.append(
            {
                "measurement": r["measurement"],
                "included_by_table_s1": int(r["included_by_table_s1"]),
                "exclusion_footnote": r["exclusion_footnote"],
                "nacl_specific_volume_cm3_g": v,
                "sample_temperature_k": t,
                "midpoint_marker_temperature_k": tm,
                "marker_above_published_table2_temperature_grid": int(tm > 1200),
                "source_pressure_gpa": source,
                "source_pressure_uncertainty_gpa": float(
                    r["pressure_reported_uncertainty_gpa"]
                ),
                "brown_table1_pressure_gpa": p1,
                "brown_table2_pressure_gpa": p2,
                "table1_minus_source_gpa": p1 - source,
                "table2_minus_source_gpa": p2 - source,
                "fixed_ringwoodite_model_pressure_gpa": ring,
                "source_minus_ringwoodite_gpa": source - ring,
                "table1_minus_ringwoodite_gpa": p1 - ring,
                "table2_minus_ringwoodite_gpa": p2 - ring,
                "table2_lattice_only_pressure_sd_gpa": abs(dp_da) * sd_a,
                "table2_dp_dt_gpa_k": dp_dt,
                "marker_temperature_matching_source_table1_k": required_source_t1,
                "marker_temperature_matching_source_table2_k": required_source_t2,
                "marker_temperature_matching_ringwoodite_table2_k": required_ring_t,
                "marker_fraction_matching_ringwoodite_table2": (required_ring_t - 300)
                / (t - 300)
                if t != 300 and required_ring_t is not None
                else None,
                "half_sample_temperature_diagnostic_k": half_t,
                "half_sample_temperature_table1_minus_source_gpa": model.table1(
                    v, half_t
                )
                - source,
                "half_sample_temperature_table2_minus_source_gpa": model.table2(
                    v, half_t, extend=True
                )
                - source,
                "table2_whole_sample_temperature_minus_ringwoodite_gpa": model.table2(
                    v, t, extend=True
                )
                - ring,
            }
        )
    accepted = [r for r in rows if r["included_by_table_s1"]]
    fields = [
        "table1_minus_source_gpa",
        "table2_minus_source_gpa",
        "source_minus_ringwoodite_gpa",
        "table1_minus_ringwoodite_gpa",
        "table2_minus_ringwoodite_gpa",
        "half_sample_temperature_table1_minus_source_gpa",
        "half_sample_temperature_table2_minus_source_gpa",
    ]
    statistics = {}
    for field in fields:
        values = np.array([r[field] for r in accepted])
        statistics[field] = {
            "mean": float(np.mean(values)),
            "rms": float(np.sqrt(np.mean(values**2))),
            "minimum": float(min(values)),
            "maximum": float(max(values)),
        }
    return {
        "audit_date": "2026-10-16",
        "source_doi": "10.1063/1.371596",
        "scope": "All18 TableS1 readings retained,14 accepted; original pressures, mask and fixed2011 EOS coefficients unchanged. No coefficient fitting.",
        "methods": {
            "conversion": "v=a^3*1e-24*N_A/(4*M), cm3/g; B1 Z=4, M=58.44g/mol, N_A=6.02214076e23/mol; no assumed ambient lattice parameter.",
            "table1": "Brown Eqs1-4, TableI p5804, MATLAB SPLINE prescription: not-a-knot cubic splines of P0,gamma,theta; P=P0+gamma*E_D/(v*M)/1000, E_D=18RT*(T/theta)^3*integral_0^(theta/T)[z^3/(exp(z)-1)]dz, J/molNaCl. No zero-point term added to the stated zero-K pressure.",
            "table2": "Published TableII p5805 directly interpolated by not-a-knot cubic splines in specific volume and temperature. Above1200K explicitly extend from1200K with TableI gamma/theta Debye energy increment. Below300K inverse-temperature diagnostics extend analogously from300K. No volume extrapolation.",
            "temperature": "Tsample=(Tu+Td)/2; Tmarker=(Tsample+300)/2 from2011 paragraph4. Cycle2 sample means confirmed byFigureS1; first-cycle averaging remains inferred.",
            "uncertainty": "Only lattice SD propagated by local dP/da. No rowwise temperature uncertainty, calibration covariance or joint sample/marker covariance recovered. Source pressure parentheses have unspecified confidence. These are descriptive central-value comparisons, not significance tests.",
            "half_temperature_diagnostic": "For heated rows only, Tmarker=Tsample/2,150K below the stated midpoint. Numerical hypothesis only; neither adopted nor claimed as author workflow. Quenched/preheat states kept at300K.",
            "inverse_temperatures": "Solve Brown(VNaCl,Tmarker)=Psource or fixedPmodel over1-4000K; null means no root in that interval. Fractions are null for300K sample states. These are model-conditional diagnostics, not measured gradients. Whole-sample-T scenario is also a diagnostic, not a proposed correction.",
        },
        "constants": {
            "nacl_molar_mass_g_mol": MOLAR_MASS,
            "atoms_per_formula": 2,
            "formula_units_per_cell": 4,
            "gas_constant_j_mol_k": R,
        },
        "primary_table_consistency": {
            "table1_rows": 21,
            "table2_rows": 63,
            "table2_pressure_cells_checked": len(checks),
            "maximum_abs_table1_minus_table2_gpa": max(
                abs(c["table1_minus_table2_gpa"]) for c in checks
            ),
            "maximum_abs_300k_anchored_thermal_increment_difference_gpa": float(
                max(abs(np.array(thermal_differences)))
            ),
            "assessment": "TableI plus normalized Debye energy fails to reproduce TableII at its printed precision: a mostly temperature-independent pressure offset remains. Thermal increments agree within0.016GPa over complete rows. Cause unresolved; preserve both, do not repair either published table or claim exact author calculator recovery.",
            "checks": checks,
        },
        "accepted_14_statistics_gpa": statistics,
        "table2_temperature_extension_measurements": [
            r["measurement"]
            for r in rows
            if r["marker_above_published_table2_temperature_grid"]
        ],
        "calibration_explains_hot_residuals": False,
        "author_reduction_reproduced": False,
        "rows": rows,
    }


def plot(result: dict, path: Path) -> None:
    import matplotlib.pyplot as plt

    rows = [r for r in result["rows"] if r["included_by_table_s1"]]
    t = [r["sample_temperature_k"] for r in rows]
    with plt.rc_context(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    ):
        fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), layout="constrained")
        for ax in axes:
            ax.axhline(0, c="#697680", ls="--", lw=1)
            ax.set_xlabel("Mean sample temperature (K)")
        axes[0].scatter(
            t,
            [r["table1_minus_source_gpa"] for r in rows],
            label="Table I + Debye",
            c="#9a4d1b",
            marker="s",
        )
        axes[0].scatter(
            t,
            [r["table2_minus_source_gpa"] for r in rows],
            label="Table II + flagged extension",
            c="#20786b",
        )
        axes[0].set_ylabel("Brown replay - published pressure (GPa)")
        axes[0].set_title("Calibration replay at the stated marker midpoint")
        axes[1].scatter(
            t,
            [r["source_minus_ringwoodite_gpa"] for r in rows],
            label="Published Table S1 pressures",
            c="#284b7b",
        )
        axes[1].scatter(
            t,
            [r["table1_minus_ringwoodite_gpa"] for r in rows],
            label="Brown Table I pressures",
            c="#9a4d1b",
            marker="s",
        )
        axes[1].scatter(
            t,
            [r["table2_minus_ringwoodite_gpa"] for r in rows],
            label="Brown Table II pressures",
            c="#20786b",
        )
        axes[1].set_ylabel("Pressure - fixed ringwoodite EOS (GPa)")
        axes[1].set_title("The high-temperature trend remains")
        for ax in axes:
            ax.legend(fontsize=8, loc="lower left")
        fig.suptitle(
            "14 accepted rows; published 2011 EOS coefficients unchanged", fontsize=11
        )
        fig.savefig(path, dpi=180)
        plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output", type=Path, default=STAGE / f"{PREFIX}brown-replay.json"
    )
    parser.add_argument(
        "--rows", type=Path, default=STAGE / f"{PREFIX}brown-replay.csv"
    )
    parser.add_argument("--plot", type=Path)
    parser.add_argument("--source-dir", type=Path)
    args = parser.parse_args()
    result = replay()
    if args.source_dir:
        result["source_pdf_transcription_check"] = source_pdf_check(args.source_dir)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n")
    stream = io.StringIO()
    writer = csv.DictWriter(
        stream, fieldnames=list(result["rows"][0]), lineterminator="\n"
    )
    writer.writeheader()
    writer.writerows(result["rows"])
    args.rows.write_text(stream.getvalue())
    if args.plot:
        plot(result, args.plot)
    print(json.dumps(result["accepted_14_statistics_gpa"], indent=2))


if __name__ == "__main__":
    main()
