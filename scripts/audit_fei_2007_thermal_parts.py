"""Compare Fei (2007) thermal diagnostics across Au, Pt, NaCl-B2 and Ne.

Run with ``python -m scripts.audit_fei_2007_thermal_parts``. This writes only
the diagnostic report; published EOS records and the catalog ledger are intact.
"""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos import get_eos_record
from peritheos.materials import NACL_B2_FEI_2007
from scripts.reproduce_fei_2007_platinum import observations as platinum_observations

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
INPUT = ROOT / "docs/data/fei-2007-thermal-inputs.json"
OUTPUT = ROOT / "docs/data/fei-2007-thermal-comparison.json"
# Conventional cell A^3, GPa, dimensionless, dimensionless, dimensionless, K.
PARAMETERS = {
    "gold": (67.85, 167.0, 6.0, 2.97, 0.6, 170.0),
    "platinum": (60.38, 277.0, 5.08, 2.72, 0.5, 230.0),
    "nacl_b2": (41.35, 26.86, 5.25, 1.70, 0.5, 290.0),
    "neon": (88.967, 1.16, 8.23, 2.05, 0.6, 75.1),
}
NAMES = ("V0", "K0", "K0_prime", "gamma0", "q", "theta0")


def pressure(volume, temperature, material, coefficients=None, law="printed"):
    """Independent Vinet-MGD evaluation with Debye quadrature and SI units."""
    v0, k0, kp, gamma0, q, theta0 = (
        PARAMETERS[material] if coefficients is None else coefficients
    )
    volume, temperature = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    logarithm = np.log(volume / v0)
    gamma = gamma0 * np.exp(q * logarithm)
    exponent = (
        -gamma * logarithm
        if law == "printed"
        else -gamma0 * logarithm
        if abs(q) < 1e-10
        else -gamma0 * np.expm1(q * logarithm) / q
    )
    theta = theta0 * np.exp(exponent)

    def energy(t, th):
        y = th / t
        integral = quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-11, epsrel=1e-11)[
            0
        ]
        return 9 * R * t * integral / y**3

    increment = np.array(
        [
            energy(t, th) - energy(300.0, th)
            for t, th in zip(temperature.flat, theta.flat)
        ]
    ).reshape(volume.shape)
    # Energy per mole of formula units: NaCl has 2 atoms and 1 FU/cell;
    # fcc Au/Pt/Ne each have 1 atom per FU and 4 FU/cell.
    atoms, units = (2, 1) if material == "nacl_b2" else (1, 4)
    molar_volume = volume * N_A * 1e-30 / units
    x = np.exp(logarithm / 3)
    cold = 3 * k0 * (1 - x) * np.exp(1.5 * (kp - 1) * (1 - x)) / x**2
    return cold + gamma * atoms * increment / molar_volume / 1e9


def fit(
    volume, temperature, target, material, free=(4,), law="printed", coefficients=None
):
    """Equal pressure weights, source V0/K0/theta0 fixed; no source covariance."""
    base = np.array(PARAMETERS[material] if coefficients is None else coefficients)
    indices = list(free)

    def unpack(x):
        coefficients = base.copy()
        coefficients[indices] = x
        return coefficients

    result = least_squares(
        lambda x: pressure(volume, temperature, material, unpack(x), law) - target,
        base[indices],
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=2000,
    )
    dof = len(volume) - len(indices)
    covariance = np.linalg.inv(result.jac.T @ result.jac) * (
        result.fun @ result.fun / dof
    )
    return {
        "observations": len(volume),
        "law": law,
        "free_parameters": [NAMES[i] for i in indices],
        "fixed_parameters": {
            NAMES[i]: float(base[i]) for i in range(6) if i not in indices
        },
        "parameters": dict(zip([NAMES[i] for i in indices], result.x.tolist())),
        "conditional_standard_errors": dict(
            zip([NAMES[i] for i in indices], np.sqrt(np.diag(covariance)).tolist())
        ),
        "published_rmse_gpa": float(
            np.sqrt(
                np.mean(
                    (pressure(volume, temperature, material, law=law) - target) ** 2
                )
            )
        ),
        "rmse_gpa": float(np.sqrt(np.mean(result.fun**2))),
        "solver_success": bool(result.success),
    }


def read_csv(name):
    with (DATA / name).open() as stream:
        return list(csv.DictReader(stream))


def reproduce():
    inputs = json.loads(INPUT.read_text(encoding="utf-8"))
    hirose = inputs["hirose_2006"]["observations"]
    au_rows = read_csv("gold-fei-2004-table1.csv")
    av, at, ap = [
        np.array([float(row[key]) for row in au_rows])
        for key in ("volume_a3", "temperature_k", "pressure_gpa")
    ]
    hot = [row for row in hirose if row["temperature_k"] > 300]
    av_all = np.r_[av, [row["gold_a_angstrom"] ** 3 for row in hot]]
    at_all = np.r_[at, [row["temperature_k"] for row in hot]]
    ap_all = np.r_[ap, [row["mgo_pressure_gpa"] for row in hot]]
    cold = read_csv("gold-dewaele-2004-table1-compression.csv")
    new_cold = read_csv("gold-fei-2007-figure1-digitized.csv")
    hirose_cold = [row for row in hirose if row["temperature_k"] == 300]
    cold_v = np.r_[
        [4 * float(row["atomic_volume_a3"]) for row in cold],
        [float(row["volume_a3_conventional_cell"]) for row in new_cold],
        [row["gold_a_angstrom"] ** 3 for row in hirose_cold],
    ]
    cold_p = np.r_[
        [float(row["ruby_pressure_revised_gpa"]) for row in cold],
        [float(row["pressure_gpa"]) for row in new_cold],
        [row["mgo_pressure_gpa"] for row in hirose_cold],
    ]
    pv, pt, pp, _, _ = platinum_observations()
    results = {
        "gold": {
            "fei_2004_only_q_diagnostic": fit(av, at, ap, "gold"),
            "q_only_diagnostic": fit(av_all, at_all, ap_all, "gold"),
            "integrated_law_q_only_sensitivity": fit(
                av_all, at_all, ap_all, "gold", law="integrated"
            ),
            "joint_diagnostic": fit(
                np.r_[cold_v, av_all],
                np.r_[np.full(len(cold_v), 300.0), at_all],
                np.r_[cold_p, ap_all],
                "gold",
                free=(2, 3, 4),
            ),
            "scope": "26 Fei (2004) Au-MgO hot rows plus two Hirose (2006) Run 2 hot rows. Use reported Speziale-MgO pressures, not Au-derived pressures. The joint sensitivity additionally includes 37 Dewaele cold rows, six Fei Figure 1 cold rows and one Hirose cold row; original Dewaele 298 K rows are assigned 300 K as a diagnostic choice. Source selection, unrounded inputs and weights are not recovered.",
        },
        "platinum": {
            "q_only_diagnostic": fit(pv, pt, pp, "platinum"),
            "joint_diagnostic": fit(pv, pt, pp, "platinum", free=(2, 3, 4)),
            "q_fixed_at_published_joint_sensitivity": fit(
                pv, pt, pp, "platinum", free=(2, 3)
            ),
            "q_fixed_at_published_hot_sensitivity": fit(
                pv[pt > 300], pt[pt > 300], pp[pt > 300], "platinum", free=(3,)
            ),
            "source_q_fit_status": "Fei (2004), Section 3.2 and Table 3, presents q=0.5(5) among optimized Pt model parameters. Fei (2007) retains the same value/error and describes fitted model parameters but does not explicitly enumerate q as a free parameter in the 2007 optimization. A newly optimized versus retained Pt q cannot be distinguished with certainty from this wording alone.",
            "scope": "Existing 78-row Pt reconstruction: 36 Dewaele cold rows and 42 Fei (2004) paired Au-Pt rows, with pressures re-reduced from measured Au volumes using the published Fei (2007) Au scale. See reproduce_fei_2007_platinum.py. Au diagnostic coefficients are never substituted into this calibration.",
        },
    }
    hot_mask = pt > 300
    for label, cold_mask in (
        ("staged_dewaele_rt_then_fei_hot", np.arange(len(pv)) < 36),
        ("staged_all_rt_then_fei_hot", pt == 300),
    ):
        cold_fit = fit(
            pv[cold_mask], pt[cold_mask], pp[cold_mask], "platinum", free=(2,)
        )
        coefficients = np.array(PARAMETERS["platinum"])
        coefficients[2] = cold_fit["parameters"]["K0_prime"]
        thermal_fit = fit(
            pv[hot_mask],
            pt[hot_mask],
            pp[hot_mask],
            "platinum",
            free=(3, 4),
            coefficients=coefficients,
        )
        coefficients[3:5] = [
            thermal_fit["parameters"]["gamma0"],
            thermal_fit["parameters"]["q"],
        ]
        results["platinum"][label] = {
            "cold_stage": cold_fit,
            "thermal_stage": thermal_fit,
            "all_78_rmse_gpa": float(
                np.sqrt(np.mean((pressure(pv, pt, "platinum", coefficients) - pp) ** 2))
            ),
            "qualification": "Fit K0-prime to RT data first with V0/K0 fixed, then fix that cold EOS and fit gamma0/q to the 35 Fei hot rows with theta0 fixed. Equal pressure weights; thermal-stage errors do not propagate cold-stage or Au-calibration uncertainty. All target pressures use the published Au coefficients, not diagnostic refits.",
        }
    for material, figure in inputs["figures"].items():
        rows = figure["observations"]
        v, t, p = [
            np.array([row[key] for row in rows])
            for key in ("volume_a3", "temperature_k", "pressure_gpa")
        ]
        bounds = {}
        for quantity in ("pressure", "volume"):
            (a, b), (c, d) = figure["axis_calibration_pixels"][quantity]
            bounds[quantity] = figure["coordinate_sensitivity_bound_px"] * abs(
                (d - b) / (c - a)
            )
        extreme_q = [
            fit(
                v + sign * bounds["volume"], t, p + sign * bounds["pressure"], material
            )["parameters"]["q"]
            for sign in (-1, 1)
        ]
        results[material] = {
            "q_only_diagnostic": fit(v, t, p, material),
            "integrated_law_q_only_sensitivity": fit(
                v, t, p, material, law="integrated"
            ),
            "coordinate_offset_q_range": sorted(extreme_q),
            "coordinate_bounds": bounds,
            "scope": "Figure-only partial reconstruction, published Vinet cold coefficients/gamma0/theta0 fixed, q fitted without weights. No measured Pt volumes are available to independently re-reduce the source's Pt-calibrated pressures. Coordinate sensitivity is not a confidence interval. Two touching NaCl markers are excluded; all six separable Ne hot markers are included.",
        }
    records = {
        "gold": get_eos_record("gold_fei_2007_vinet_2"),
        "platinum": get_eos_record("platinum_fei_2007_vinet_mgd"),
        "nacl_b2": NACL_B2_FEI_2007,
        "neon": get_eos_record("neon_fcc_fei_2007_vinet_2"),
    }
    for material, record in records.items():
        ratios = (
            [0.98, 0.85, 0.75]
            if material in ("gold", "platinum")
            else [0.7, 0.6, 0.5]
            if material == "nacl_b2"
            else [0.5, 0.4, 0.3]
        )
        v = PARAMETERS[material][0] * np.array(ratios)[:, None]
        t = np.array([300.0, 1000.0])[None, :]
        delta = float(np.max(abs(record.pressure(v, t) - pressure(v, t, material))))
        assert delta < 1e-8, (material, delta)
        results[material]["native_equation_max_difference_gpa"] = delta
        results[material]["theta_identity_factor_at_half_volume"] = float(
            1 + PARAMETERS[material][4] * np.log(0.5)
        )
    names = [
        "gold-fei-2004-table1.csv",
        "gold-dewaele-2004-table1-compression.csv",
        "gold-fei-2007-figure1-digitized.csv",
        "platinum-dewaele-2004-table1-compression.csv",
        "platinum-fei-2004-table2.csv",
    ]
    hashes = {
        (DATA / name).relative_to(ROOT).as_posix(): hashlib.sha256(
            (DATA / name).read_bytes()
        ).hexdigest()
        for name in names
    }
    hashes[INPUT.relative_to(ROOT).as_posix()] = hashlib.sha256(
        INPUT.read_bytes()
    ).hexdigest()
    for figure in inputs["figures"].values():
        digest = hashlib.sha256((ROOT / figure["image_path"]).read_bytes()).hexdigest()
        assert digest == figure["image_sha256"]
        hashes[figure["image_path"]] = digest
    return {
        "doi": inputs["doi"],
        "primary_source": inputs["primary_source"],
        "input_sha256": hashes,
        "qualification": "Conditional diagnostics, not a recovery of author fitting procedure. All fits use equal pressure weights; errors are residual-scaled linearized errors conditional on the adopted cold scale, pressure standards and data selection. Calibration covariance, source weights and full experimental uncertainty are unavailable. The printed theta law implies -dln(theta)/dln(V)=gamma*(1+q*ln(V/V0)), so its thermodynamic identity issue is common to all four models. Published coefficients remain unchanged.",
        "materials": results,
    }


if __name__ == "__main__":
    result = reproduce()
    OUTPUT.write_text(
        json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    for name, data in result["materials"].items():
        diagnostic = data["q_only_diagnostic"]
        print(name, "q", diagnostic["parameters"]["q"], "RMSE", diagnostic["rmse_gpa"])
