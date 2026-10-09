"""Validate Au pressure equations, source curves and qualified reconstructions.

Run ``python -m scripts.audit_gold_pvt --check`` to verify the archived audit.
Equation parity, agreement with calculated source outputs, measured-data
comparisons and reproduction of an author's original fit are separate results.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad

from peritheos import get_eos_record

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
CURVES = ROOT / "docs/data/fei-2007-gold-source-curves.json"
OUTPUT = ROOT / "docs/data/gold-pvt-validation.json"
RECORD = "gold_fei_2007_vinet_2"
PARAMETERS = dict(
    V0=67.85, K0=167.0, K0_prime=6.0, theta0=170.0, gamma0=2.97, q=0.6, Tr=300.0, n=1
)


def fei_pressure(volume, temperature):
    """Fei Table 1/Eqs. 2-3 in SI, independent of the library evaluator.

    Input is four-atom fcc cell Å³. The printed theta law uses gamma(V) as
    the exponent; it is not the integrated constant-q theta law. Electronic
    pressure is not an additional term in this published parameterization.
    """
    volume, temperature = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    if np.any(~np.isfinite(volume) | (volume <= 0)) or np.any(
        ~np.isfinite(temperature) | (temperature <= 0)
    ):
        raise ValueError("Positive finite cell volumes and temperatures required")
    p = PARAMETERS
    ratio = volume / p["V0"]
    x = ratio ** (1 / 3)
    cold = 3 * p["K0"] * (1 - x) / x**2 * np.exp(1.5 * (p["K0_prime"] - 1) * (1 - x))
    gamma = p["gamma0"] * ratio ** p["q"]
    theta = p["theta0"] * ratio ** (-gamma)

    def energy(t, th):
        y = th / t
        integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11, epsrel=1e-11)[
            0
        ]
        return 9 * R * t * integral / y**3

    delta_energy = np.array(
        [
            energy(t, th) - energy(p["Tr"], th)
            for t, th in zip(temperature.flat, theta.flat)
        ]
    ).reshape(volume.shape)
    molar_volume_m3 = volume * N_A * 1e-30 / 4
    pressure = cold + gamma * delta_energy / molar_volume_m3 / 1e9
    return float(pressure) if pressure.ndim == 0 else pressure


def metrics(residual):
    residual = np.asarray(residual)
    return {
        "states": int(residual.size),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_difference_gpa": float(np.max(abs(residual))),
    }


def measured_comparison(filename, *, cold=False):
    with (DATA / filename).open() as stream:
        rows = list(csv.DictReader(stream))
    volume = np.array(
        [
            4 * float(row["atomic_volume_a3"]) if cold else float(row["volume_a3"])
            for row in rows
        ]
    )
    temperature = (
        np.full(len(rows), 300.0)
        if cold
        else np.array([float(row["temperature_k"]) for row in rows])
    )
    pressure = np.array(
        [
            float(row["ruby_pressure_revised_gpa"] if cold else row["pressure_gpa"])
            for row in rows
        ]
    )
    return {
        **metrics(fei_pressure(volume, temperature) - pressure),
        "calibration": "Dewaele revised ruby" if cold else "Speziale (2001) MgO",
        "temperature_k": [float(temperature.min()), float(temperature.max())],
        "pressure_gpa": [float(pressure.min()), float(pressure.max())],
        "qualification": (
            "Comparison with source fit inputs, not withheld physical validation; "
            "original uncertainty columns remain in the CSV. No refit is made."
            + (
                " Original Dewaele 298 K rows are evaluated at the paper's 300 K "
                "reference isotherm, an explicit comparison convention."
                if cold
                else ""
            )
        ),
        "points": [
            {
                "source_row_index": i,
                "cell_volume_a3": float(v),
                "temperature_k": float(t),
                "source_pressure_gpa": float(p),
                "model_pressure_gpa": float(m),
            }
            for i, (v, t, p, m) in enumerate(
                zip(volume, temperature, pressure, fei_pressure(volume, temperature)),
                start=1,
            )
        ],
    }


def source_curve_comparison():
    source = json.loads(CURVES.read_text(encoding="utf-8"))
    axes = source["axis_calibration"]
    comparisons = []
    for curve in source["curves"]:
        xy = np.array(curve["points_xy_pt"])
        volume = np.polyval(axes["volume_from_y"], xy[:, 1])
        pressure = np.polyval(axes["pressure_from_x"], xy[:, 0])
        model = fei_pressure(volume, curve["temperature_k"])
        # One full printed stroke in the pressure direction is a graphical
        # comparison scale, not a confidence interval or measurement error.
        stroke_gpa = curve["linewidth_pt"] * axes["pressure_from_x"][0]
        comparisons.append(
            {
                "temperature_k": curve["temperature_k"],
                **metrics(model - pressure),
                "one_stroke_pressure_width_gpa": stroke_gpa,
                "within_one_stroke_pressure_width": bool(
                    np.all(abs(model - pressure) < stroke_gpa)
                ),
                "points": [
                    {
                        "cell_volume_a3": float(v),
                        "figure_pressure_gpa": float(p),
                        "model_minus_figure_gpa": float(d),
                    }
                    for v, p, d in zip(volume, pressure, model - pressure)
                ],
            }
        )
    return {
        "kind": "calculated_source_curves_not_observations",
        "source_pdf_sha256": source["source_pdf_sha256"],
        "location": source["location"],
        "qualification": "Vector path vertices and axis ticks from the primary PDF. "
        "Agreement at graph precision does not establish exact author "
        "numerical output parity or the original fitting objective.",
        "curves": comparisons,
    }


def audit():
    from scripts.fit_yokoo_2009_gold_thermal import (
        OUTPUT as YOKOO_OUTPUT,
    )
    from scripts.fit_yokoo_2009_gold_thermal import (
        check_reconstruction,
        reconstruct,
    )
    from scripts.reproduce_yokoo_2009_gold import reproduce

    record = get_eos_record(RECORD)
    volume = PARAMETERS["V0"] * np.linspace(0.72, 1.04, 17)[:, None]
    temperature = np.array([300.0, 500.0, 1000.0, 1473.0, 2000.0, 2173.0, 2330.0])[
        None, :
    ]
    source = fei_pressure(volume, temperature)
    library = record.pressure(volume, temperature)
    yokoo = reconstruct()
    check_reconstruction(json.loads(YOKOO_OUTPUT.read_text(encoding="utf-8")), yokoo)
    printed = reproduce()
    return {
        "scope": "Au pressure-volume-temperature only; caloric validation is not an acceptance criterion.",
        "audit_date": "2026-10-09",
        "input_sha256": {
            str(path.relative_to(ROOT)): hashlib.sha256(path.read_bytes()).hexdigest()
            for path in [
                CURVES,
                DATA / "gold-fei-2004-table1.csv",
                DATA / "gold-dewaele-2004-table1-compression.csv",
                DATA / "gold-yokoo-2009-source.json",
                DATA / "gold-yokoo-2009-table3-isochores.csv",
                DATA / "tsuchiya-kawamura-2002-table1-electronic-pressure.csv",
            ]
        },
        "fei_2007": {
            "record_identifier": RECORD,
            "status": "published_pvt_equation_verified",
            "published_parameters": PARAMETERS,
            "doi": "10.1073/pnas.0609013104",
            "source_locations": ["Table 1", "Equations 2-3", "Figure 1"],
            "independent_equation": metrics(library - source),
            "verification_grid": {
                "cell_volume_a3": volume[:, 0].tolist(),
                "temperature_k": temperature[0].tolist(),
                "qualification": "Numerical test domain, not a measured validity rectangle.",
            },
            "source_curve_comparison": source_curve_comparison(),
            "measured_data_comparisons": {
                "fei_2004_hot": measured_comparison("gold-fei-2004-table1.csv"),
                "dewaele_2004_cold": measured_comparison(
                    "gold-dewaele-2004-table1-compression.csv", cold=True
                ),
            },
            "author_fit_reproduced": False,
            "limitations": "Unrounded source coefficients, covariance, original weights "
            "and complete fitting selections are not recovered. These checks "
            "verify the published PVT equation and available comparison evidence.",
        },
        "yokoo_2009": {
            "doi": yokoo["doi"],
            "published_300k_records": list(printed["branches"]),
            "full_published_pvt_status": "not_reproduced",
            "reconstruction_status": "numerically_verified_derived_output_reconstruction",
            "experimental_observations": yokoo["experimental_observations"],
            "author_fit_reproduced": yokoo["author_fit_reproduced"],
            "primary_fit": yokoo["fits"]["primary"]["fit_states"],
            "first_liquid_holdout": yokoo["fits"]["primary"]["first_liquid_states"],
            "volume_holdout": yokoo["volume_interpolation_holdout"]["held_out_states"],
            "independent_quad_max_difference_gpa": yokoo[
                "adaptive_quad_max_difference_gpa"
            ],
            "published_coefficients_changed": False,
            "limits": "Derived reconstruction coefficients differ from the printed phonon "
            "normalization; cold Vc is inferred. Electronic interpolation is an "
            "explicit linear choice. Source phase markers and omitted cells do "
            "not certify the full rectangular domain as solid Au. The supplied "
            "electronic pressure table suffices for this bounded PVT reconstruction; "
            "missing electronic energy is not a PVT blocker.",
        },
    }


def check_saved(saved, current):
    if isinstance(current, dict):
        assert saved.keys() == current.keys()
        for key in current:
            check_saved(saved[key], current[key])
    elif isinstance(current, list):
        assert len(saved) == len(current)
        for left, right in zip(saved, current):
            check_saved(left, right)
    elif isinstance(current, float):
        np.testing.assert_allclose(saved, current, rtol=2e-7, atol=2e-7)
    else:
        assert saved == current


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = audit()
    if args.check:
        check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), result)
    else:
        OUTPUT.write_text(
            json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8"
        )
    print(
        json.dumps(
            {
                "fei_2007": result["fei_2007"]["status"],
                "yokoo_published_full_pvt": result["yokoo_2009"][
                    "full_published_pvt_status"
                ],
                "yokoo_reconstruction": result["yokoo_2009"]["reconstruction_status"],
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
