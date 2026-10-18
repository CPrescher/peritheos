#!/usr/bin/env python3
"""Build/evaluate a distinct historical-data NaCl refit, not a corrected Brown EOS.

Run from the repository root. The saved JSON contains the complete cold spline
and thermal settings; load_eos() evaluates it without refitting or reading data.
"""

from __future__ import annotations

import csv
import hashlib
import json
import platform
from dataclasses import asdict, replace
from pathlib import Path

import numpy as np
import scipy
from scipy.interpolate import BSpline
from scipy.optimize import minimize_scalar

try:
    from scripts import brown_nacl_1999_reconstruct as base
except ModuleNotFoundError:  # direct execution from scripts/
    import brown_nacl_1999_reconstruct as base

ID = "nacl_b1_historical_data_2026_refit"
DEST = base.OUT / f"{ID}.json"


class SavedEOS(base.Model):
    """Standalone saved candidate; volumes in cm3/g, pressure in GPa, T in K."""

    def __init__(self, record):
        self.settings = base.Settings(**record["settings"])
        self.thermal = base.Thermal(self.settings)
        cold = record["cold_spline"]
        self.spline = BSpline(cold["knots_x"], cold["coefficients_gpa"], 3)
        self.v = np.array(record["evaluation_domain"]["specific_volume_cm3_g"])
        self.temperature_bounds = record["evaluation_domain"]["temperature_k"]

    def pressure(self, v, T):
        v, T = np.broadcast_arrays(np.asarray(v), np.asarray(T))
        if np.any((v < min(self.v)) | (v > max(self.v))):
            raise ValueError("Volume outside the saved reconstruction domain")
        if np.any((T < self.temperature_bounds[0]) | (T > self.temperature_bounds[1])):
            raise ValueError("Temperature outside the saved reconstruction domain")
        return super().pressure(v, T)


def load_eos(path=DEST):
    record = json.loads(Path(path).read_text(encoding="utf-8"))
    if record["id"] != ID or record["format"] != "peritheos-standalone-refit-v1":
        raise ValueError("Not a supported historical-data NaCl reconstruction")
    expected = dict(
        R=base.R, M=base.M, N=base.N, VREF=base.VREF, TREF=base.TREF, PREF=base.PREF
    )
    if record["constants"] != expected:
        raise ValueError("Saved constants differ from this evaluator")
    return SavedEOS(record)


def validation(model):
    rows = []
    for r in base.read("yamamoto-1987-table2"):
        T = float(r["temperature_k"])
        p = model.properties(base.PREF, T)
        rows.append(
            dict(
                T_k=T,
                model_ks_gpa=p["ks"],
                source_ks_gpa=float(r["ks_kbar"]) * 0.1,
                model_cp_j_g_k=p["cp"],
                source_cp_j_g_k=float(r["cp_1e6_erg_g_k"]) * 0.1,
            )
        )
    cp = []
    for r in base.read("janaf-1971-nacl-cr-cp"):
        T = float(r["temperature_k"])
        if 298 <= T <= 1000:
            observed = float(r["cp_j_mol_k"]) / float(r["historical_molar_mass_g_mol"])
            pred = model.properties(base.PREF, T)["cp"]
            cp.append(dict(T_k=T, source_cp_j_g_k=observed, model_cp_j_g_k=pred))
    v = np.linspace(model.v.min() + 1e-5, model.v.max() - 1e-5, 1001)
    d1 = model.spline(v / base.VREF - 1, nu=1)
    d2 = model.spline(v / base.VREF - 1, nu=2)
    kp = -1 - v / base.VREF * d2 / d1
    return dict(
        yamamoto=rows,
        janaf=cp,
        yamamoto_max_relative_ks_difference=max(
            abs(r["model_ks_gpa"] / r["source_ks_gpa"] - 1) for r in rows
        ),
        janaf_max_relative_cp_difference=max(
            abs(r["model_cp_j_g_k"] / r["source_cp_j_g_k"] - 1) for r in cp
        ),
        cold_kt_range_gpa=[
            float(min(-v / base.VREF * d1)),
            float(max(-v / base.VREF * d1)),
        ],
        cold_kprime_range=[float(min(kp)), float(max(kp))],
        initial_pressure_error_gpa=float(
            model.pressure(base.VREF, base.TREF) - base.PREF
        ),
        independence="Elastic moduli and Cp held out of fitting/model choice; gamma0 adopts Spetzler's ambient value. These sources share historical auxiliary data.",
    )


def build():
    # Deliberately low dimensional: two internal simple knots, six coefficients
    # total. Smoothness/basis choices are analyst assumptions, not Brown settings.
    settings = base.Settings(
        n_basis=6, smoothing=1e-6, join_width=0.025, reference_anchor=True
    )

    def objective(q):
        m = base.Model(replace(settings, q=float(q)))
        return float(np.mean(m.residual[m.groups == "BK"] ** 2))

    result = minimize_scalar(
        objective, bounds=(0.5, 2.5), method="bounded", options={"xatol": 1e-7}
    )
    if not result.success:
        raise RuntimeError("Conditional q fit failed")
    settings = replace(settings, q=float(result.x))
    model = base.Model(settings)
    files = [
        base.DATA / (base.PREFIX + name)
        for name in [
            "boehler-kennedy-1980-table1.csv",
            "fritz-1971-table1.csv",
            "published-coefficients.json",
        ]
    ]
    sources = {
        str(p.relative_to(base.ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
        for p in files
    }
    record = dict(
        format="peritheos-standalone-refit-v1",
        id=ID,
        label="NaCl B1 historical-data reconstruction (2026)",
        record_kind="refit",
        status="staged_independent_reconstruction; not a validated replacement pressure standard",
        published_eos_status="Brown (1999) cannot be reproduced from the available published material; author acknowledges an energy-term error",
        relation_to_brown="Uses the historical experimental sources and quasiharmonic reduction strategy discussed by Brown; neither his Tables I/II nor an inferred correction enters fitting or model selection.",
        material=dict(
            formula="NaCl",
            phase="B1",
            space_group="Fm-3m",
            formula_units_per_conventional_cell=4,
        ),
        determination_method="hybrid",
        constants=dict(
            R=base.R, M=base.M, N=base.N, VREF=base.VREF, TREF=base.TREF, PREF=base.PREF
        ),
        units=dict(v="cm3/g", pressure="GPa", energy="J/g", temperature="K"),
        settings=asdict(settings),
        equations=dict(
            pressure="P(v,T)=P_cold(v)+gamma(v)*e_D(T,theta(v))/(1000*v)",
            cold="Cubic B-spline in x=v/VREF-1; full knots and coefficients below",
            gamma="gamma0*exp(q*(-w*logaddexp(0,-log(v/VREF)/w)))",
            theta="theta_ref*exp(-integral_VREF^v gamma(u)/u du); base.Thermal numerical quadrature",
            e_D="9*n*R*T/M*(T/theta)^3*integral_0^(theta/T) x^3/(exp(x)-1) dx; n=2, all three modes per atom",
            cold_energy="e_cold(VREF)=-e_D(TREF,theta_ref); e_cold(v)=e_cold(VREF)-1000*integral_VREF^v P_cold(u)du",
            energy_reference="Total initial energy is zero. Cold branch represents the fitted effective ground-state contribution, including any zero-point pressure; do not add an extra zero-point term.",
        ),
        cold_spline=dict(
            degree=3,
            knots_x=model.basis_knots.tolist(),
            coefficients_gpa=model.coefficients.tolist(),
        ),
        evaluation_domain=dict(
            specific_volume_cm3_g=[float(min(model.v)), float(max(model.v))],
            temperature_k=[293.0, 1023.15],
            interpretation="Numerical evaluation bounds only, not a rectangular experimental or phase-stability range; no volume/temperature extrapolation in saved evaluator.",
            observations="BK 25-500 C, up to 3.2 GPa; expansion polynomial samples 50-750 C; shock reductions to printed 26.4 GPa with calculated shock temperatures.",
        ),
        fit_provenance=dict(
            script="scripts/brown_nacl_1999_independent_eos.py",
            engine="scripts/brown_nacl_1999_reconstruct.py",
            software=dict(
                python=platform.python_version(),
                numpy=np.__version__,
                scipy=scipy.__version__,
            ),
            dataset=sources,
            source_metadata="peritheos/data/datasets/brown-nacl-1999-source.json",
            row_selection="All 82 BK rows and all 43 Fritz impedance-matched velocity pairs; 15 explicitly synthetic Enck-Dommel polynomial samples at 50 C spacing; no source table pressure/volume overwritten.",
            inner_objective="Iterate RH cold-energy reduction; at each step minimize sum_i(P_cold(x_i)-reduced_target_i)^2 + lambda*integral (d2P_cold/dx2)^2 dx, subject to exact initial P,V,T.",
            outer_objective="Minimize BK pressure mean-square residual versus q in [0.5,2.5], refitting the cold curve each time; conditional staged optimization, not a simultaneous likelihood fit.",
            weights="Equal per reduced-pressure row in GPa; polynomial samples are not independent measurements; source accuracy statements are not Gaussian sigma weights.",
            varied_parameters=[
                "q",
                "six cold-spline coefficients with one exact reference constraint",
            ],
            fixed_parameters="gamma0=1.63 (rounded Spetzler); theta_ref=276 K, rho_ref=2.1645 g/cm3, TREF=293 K (Fritz); gamma join width .025, six basis functions, smoothing 1e-6, and q bounds are analyst choices.",
            model_choice="Low dimensional smooth candidate to avoid the diagnostic natural-spline derivative oscillations; not optimized against Brown tables or held-out thermodynamic values. No assertion that these choices recover Brown's workflow.",
            coefficient_uncertainty=None,
            covariance=None,
            uncertainty_note="Uncertainties are unavailable; structural sensitivity below is not a confidence interval.",
        ),
        fit_statistics=model.stats(),
        validation=validation(model),
        limitations=[
            "Original Enck-Dommel individual data and Brown's private Fritz extended tables are unavailable.",
            "Historical pressure-scale assumptions remain; no modern recalibration or anharmonic correction.",
            "No source covariance or statistically justified coefficient uncertainty; no catalog parent invented.",
            "Independent reconstruction, not a correction to Brown (1999) or a reproduction of his EOS.",
        ],
    )
    record["structural_sensitivity_fixed_q"] = {}
    for name, alt in [
        ("eight_basis", replace(settings, n_basis=8)),
        ("lambda_divided_by_10", replace(settings, smoothing=1e-7)),
        ("lambda_times_10", replace(settings, smoothing=1e-5)),
    ]:
        m = base.Model(alt)
        record["structural_sensitivity_fixed_q"][name] = dict(
            settings=asdict(alt),
            fit_statistics=m.stats(),
            pressure_at_v03709_t300_gpa=float(m.pressure(0.3709, 300)),
        )
    DEST.write_text(json.dumps(record, indent=2, allow_nan=False) + "\n")
    saved = load_eos()
    with DEST.with_suffix(".pvt.csv").open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(["P_GPa", "T_K", "v_cm3_g", "KT_GPa", "KS_GPa", "Cp_J_g_K"])
        for T in [300.0, 500.0, 700.0, 1000.0]:
            for P in [base.PREF, 1.0, 3.0, 5.0, 10.0, 15.0, 20.0]:
                p = saved.properties(P, T)
                writer.writerow([P, T, p["v"], p["kt"], p["ks"], p["cp"]])
    print(
        json.dumps(
            dict(
                id=ID,
                q=settings.q,
                rms_gpa=model.stats()["reduced_pressure_rms_gpa"],
                validation=record["validation"]
                | {"yamamoto": "see JSON", "janaf": "see JSON"},
            ),
            indent=2,
        )
    )
    return record


if __name__ == "__main__":
    build()
