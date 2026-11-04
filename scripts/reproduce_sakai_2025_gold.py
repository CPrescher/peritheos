"""Conditional Sakai (2025) Au numerical parity, separate from Table S9.

No author weights/covariance are inferred. Native and independent evaluators
are checked against each other. Fits below use explicitly assumed equal
pressure or volume residuals, and never update production EOS coefficients.
"""

from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from peritheos.eos.rt import RydbergStacey
from scripts.reproduce_sakai_2025_nine_materials import PARAMETERS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"


def read_rows(filename):
    with (DATA / filename).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def pressure(volume, parameters):
    """Article Equation 1, evaluated independently of the production model."""
    v0, k0, kp, kinf = parameters
    x = np.cbrt(np.asarray(volume, dtype=float) / v0)
    return (
        3
        * k0
        * x ** (-3 * kinf)
        * (1 - x)
        * np.exp((1.5 * kp - 3 * kinf + 0.5) * (1 - x))
    )


def volume_at_pressure(p, parameters):
    v0 = parameters[0]
    return brentq(lambda v: float(pressure(v, parameters)) - p, v0 * 0.1, v0 * 1.1)


def ruby2020_from_mao1986(p):
    """Rounded Mao pressure -> inferred shift -> Shen (2020) Equation 3.

    Raw R1 wavelengths are absent, so this is a rounding-limited coordinate
    reconstruction, not recovery of original optical measurements.
    """
    shift = np.expm1(np.log1p(np.asarray(p) * 7.665 / 1904) / 7.665)
    return 1870 * shift * (1 + 5.63 * shift)


def fit_one_parameter(v, p, material, residual_kind="pressure"):
    published = PARAMETERS[material]
    v, p = np.asarray(v), np.asarray(p)

    def residual(kp):
        parameters = (published[0], published[1], kp[0], published[3])
        if residual_kind == "pressure":
            return pressure(v, parameters) - p
        return np.array([volume_at_pressure(value, parameters) for value in p]) - v

    fit = least_squares(
        residual, [published[2]], bounds=(1, 10), xtol=1e-12, ftol=1e-12, gtol=1e-12
    )
    return {
        "rows": len(v),
        "assumed_residual": f"equal_{residual_kind}_residuals",
        "k0_prime": float(fit.x[0]),
        "published_k0_prime": published[2],
        "difference": float(fit.x[0] - published[2]),
        "rmse": float(np.sqrt(np.mean(fit.fun**2))),
        "residual_unit": "GPa" if residual_kind == "pressure" else "angstrom^3",
        "optimizer_success": bool(fit.success),
        "author_fit_parity": False,
        "author_covariance_recovered": False,
    }


def reproduce():
    pairs = read_rows("sakai-2025-au-calibration-paired-volumes.csv")
    reductions = []
    native_difference = []

    def add(source, row, au_v, calibrant, cal_v, included=True, note=""):
        params = PARAMETERS[calibrant]
        p = float(pressure(cal_v, params))
        native_difference.append(abs(p - float(RydbergStacey(*params).pressure(cal_v))))
        reductions.append(
            {
                "source": source,
                "source_row": row,
                "volume_a3": au_v,
                "calibrant": calibrant,
                "calibrant_volume_a3": cal_v,
                "calculated_pressure_gpa": p,
                "included_in_diagnostic": included,
                "selection_note": note,
            }
        )

    for row in pairs:
        if row["source_table"] == "S4":
            add(
                "Sakai 2025 Table S4",
                row["source_row"],
                float(row["volume_1_a3"]),
                "tungsten",
                float(row["volume_2_a3"]),
            )
    for index, row in enumerate(read_rows("gold-hirose-2008-table1.csv"), 1):
        if float(row["temperature_k"]) == 300:
            add(
                "Hirose 2008 Table 1",
                index,
                float(row["volume_a3"]),
                "mgo",
                float(row["mgo_volume_a3"]),
                note="300 K rows only",
            )
    excluded_re = []
    for index, row in enumerate(read_rows("gold-dewaele-2018-tables-4-5-au-re.csv"), 1):
        p = float(pressure(float(row["re_volume_a3"]), PARAMETERS["rhenium"]))
        included = p >= 200
        if not included:
            excluded_re.append(index)
        add(
            "Dewaele 2018 Tables 4-5",
            index,
            float(row["au_volume_a3"]),
            "rhenium",
            float(row["re_volume_a3"]),
            included,
            "exclude below 200 GPa on recalculated Re scale (Sakai Fig. S6)",
        )
    for index, row in enumerate(
        read_rows("gold-dewaele-2004-table1-compression.csv"), 1
    ):
        reductions.append(
            {
                "source": "Dewaele 2004 Table I",
                "source_row": index,
                "volume_a3": 4 * float(row["atomic_volume_a3"]),
                "calibrant": "IPPS-Ruby2020",
                "calibrant_volume_a3": None,
                "calculated_pressure_gpa": float(
                    ruby2020_from_mao1986(float(row["ruby_pressure_classical_gpa"]))
                ),
                "included_in_diagnostic": True,
                "selection_note": "inferred shift from rounded Mao1986 pressures",
            }
        )

    for index, row in enumerate(
        read_rows("gold-takemura-2008-table3-compression.csv"), 1
    ):
        # Invert the reported P-prime coordinate using Dorogokupets-Oganov
        # (2007) Eq. 16, P=1884*s*(1+5.5*s). This avoids an assumption about
        # the alternate Zha coordinate; both original columns remain archived.
        source_p = float(row["ruby_pressure_dorogokupets_gpa"])
        shift = 2 * source_p / (1884 * (1 + np.sqrt(1 + 22 * source_p / 1884)))
        reductions.append(
            {
                "source": "Takemura-Dewaele 2008 Table III",
                "source_row": index,
                "volume_a3": float(row["lattice_111_angstrom"]) ** 3,
                "calibrant": "IPPS-Ruby2020",
                "calibrant_volume_a3": None,
                "calculated_pressure_gpa": float(1870 * shift * (1 + 5.63 * shift)),
                "included_in_diagnostic": True,
                "selection_note": "inferred shift from rounded Dorogokupets-Oganov Eq. 16 pressures",
            }
        )

    selected = [r for r in reductions if r["included_in_diagnostic"]]
    v = [r["volume_a3"] for r in selected]
    p = [r["calculated_pressure_gpa"] for r in selected]
    conditional = [
        fit_one_parameter(v, p, "gold", kind) for kind in ("pressure", "volume")
    ]
    for fit in conditional:
        fit["conditional_numerical_parity_achieved"] = bool(
            fit["optimizer_success"]
            and round(fit["k0_prime"], 3) == fit["published_k0_prime"]
        )
        fit["uncertainty_parity_achieved"] = False
        fit["source_selection"] = (
            "all five Table S8 Au input families, 300 K only; Dewaele2018 <200 GPa excluded"
        )
        fit["missing_input"] = (
            "author weights, original ruby shifts, exact upstream regressions and unrounded reductions"
        )
        fit["calibrants"] = (
            "published rounded Sakai Table 2 W, MgO, Re; reconstructed IPPS-Ruby2020"
        )

    parity_fit = conditional[0]
    parity_achieved = parity_fit["conditional_numerical_parity_achieved"]
    assessment = {
        "status": (
            "conditional_numerical_parity_achieved"
            if parity_achieved
            else "conditional_numerical_parity_not_achieved"
        ),
        "conditional_numerical_parity_achieved": parity_achieved,
        "scope": "Au K0_prime central value at published precision",
        "criterion": "Successful independently refitted K0_prime rounds to Table 2 at three decimal places",
        "fitted_parameter": "K0_prime",
        "published_value": parity_fit["published_k0_prime"],
        "reconstructed_value": parity_fit["k0_prime"],
        "difference": parity_fit["difference"],
        "published_decimal_places": 3,
        "rows": parity_fit["rows"],
        "assumed_residual": parity_fit["assumed_residual"],
        "fixed_parameters_from_paper": dict(
            zip(
                ("V0", "K0", "K_infinity_prime"),
                (PARAMETERS["gold"][0], PARAMETERS["gold"][1], PARAMETERS["gold"][3]),
            )
        ),
        "calibrants": parity_fit["calibrants"],
        "uncertainty_parity_achieved": False,
        "full_sequential_calibration_reproduced": False,
        "author_fit_parity": False,
        "author_fit_parity_scope": "exact original fitting procedure and uncertainty",
        "original_weights_recovered": False,
        "original_covariance_recovered": False,
    }

    # A measured-volume-only two-stage branch shows that sequential calibration
    # is numerically possible; it omits the external constraints in Table S8.
    w_rows = [r for r in pairs if r["source_table"] == "S3"]
    au_rows = [r for r in pairs if r["source_table"] == "S4"]
    sequential = []
    for kind in ("pressure", "volume"):
        copper_p = pressure(
            [float(r["volume_1_a3"]) for r in w_rows], PARAMETERS["copper"]
        )
        w_fit = fit_one_parameter(
            [float(r["volume_2_a3"]) for r in w_rows], copper_p, "tungsten", kind
        )
        w_params = (
            *PARAMETERS["tungsten"][:2],
            w_fit["k0_prime"],
            PARAMETERS["tungsten"][3],
        )
        au_p = pressure([float(r["volume_2_a3"]) for r in au_rows], w_params)
        au_fit = fit_one_parameter(
            [float(r["volume_1_a3"]) for r in au_rows], au_p, "gold", kind
        )
        sequential.append(
            {
                "path": "Cu -> W (S3) -> Au (S4)",
                "w_fit": w_fit,
                "au_fit": au_fit,
                "scope": "measured rows only; external Table S8 inputs omitted",
            }
        )

    anchor = read_rows("sakai-2025-table-s14-cu-anchor-calculated.csv")
    fit_anchor = [
        r for r in anchor if float(r["fratanduono_2020_pressure_gpa"]) <= 1200
    ]
    av = np.array([float(r["volume_a3"]) for r in fit_anchor])
    ap = np.array([float(r["fratanduono_2020_pressure_gpa"]) for r in fit_anchor])
    cu_fit = least_squares(
        lambda pars: pressure(av, (*PARAMETERS["copper"][:2], *pars)) - ap,
        [5.2035, 2.0],
        xtol=1e-13,
        ftol=1e-13,
        gtol=1e-13,
    )
    au_p = np.array([r["calculated_pressure_gpa"] for r in reductions[:19]])
    au_v = np.array([r["volume_a3"] for r in reductions[:19]])
    return {
        "source_recovery": {
            "paired_rows": len(pairs),
            "rows_by_table": {
                table: sum(r["source_table"] == table for r in pairs)
                for table in ("S1", "S2", "S3", "S4", "S7")
            },
            "dewaele_2018_rows": 84,
            "dewaele_2018_excluded_below_200_gpa": excluded_re,
        },
        "full_sequential_calibration_reproduced": False,
        "original_weights_recovered": False,
        "original_covariance_recovered": False,
        "author_fit_parity": False,
        "conditional_numerical_parity_achieved": parity_achieved,
        "uncertainty_parity_achieved": False,
        "reconstruction_assessment": assessment,
        "native_vs_independent_max_difference_gpa": max(native_difference),
        "table_s4_published_curve_residual": {
            "rmse_gpa": float(
                np.sqrt(np.mean((pressure(au_v, PARAMETERS["gold"]) - au_p) ** 2))
            ),
            "max_abs_gpa": float(
                np.max(np.abs(pressure(au_v, PARAMETERS["gold"]) - au_p))
            ),
            "calculated_pressure_range_gpa": [float(min(au_p)), float(max(au_p))],
        },
        "cu_anchor_conversion_diagnostic": {
            "kind": "fit_to_calculated_external_scale_inputs",
            "rows": len(fit_anchor),
            "k0_prime": float(cu_fit.x[0]),
            "k_infinity_prime": float(cu_fit.x[1]),
            "rmse_gpa": float(np.sqrt(np.mean(cu_fit.fun**2))),
            "max_abs_gpa": float(np.max(np.abs(cu_fit.fun))),
            "selection": "Table S14 E/F; <=1200 GPa; 1209.215 GPa last row retained in archive but not fit",
            "optimizer_success": bool(cu_fit.success),
            "author_fit_parity": False,
        },
        "conditional_au_fits": conditional,
        "measurement_only_sequential_branch": sequential,
        "pressure_reductions": reductions,
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    report = reproduce()
    encoded = json.dumps(report, indent=2) + "\n"
    if args.output:
        args.output.write_text(encoded, encoding="utf-8")
    else:
        print(encoded, end="")
