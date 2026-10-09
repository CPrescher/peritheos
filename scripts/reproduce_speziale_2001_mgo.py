"""Validate printed Speziale MgO equations; distinguish later pressure reductions.

python -m scripts.reproduce_speziale_2001_mgo [--check]

No material registration, author-target replacement, or Au/Pt refitting.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
from pathlib import Path

import numpy as np
from scipy.constants import N_A, R
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/speziale-2001-mgo-reproduction.json"
RESIDUALS = ROOT / "docs/data/speziale-2001-mgo-calibration-residuals.csv"
PARAMETERS = dict(
    V0=74.71,
    K0=160.2,
    K0_prime=3.99,
    gamma0=1.524,
    q0=1.65,
    q1=11.8,
    theta0=773.0,
    Tr=300.0,
)
NODES, WEIGHTS = np.polynomial.legendre.leggauss(64)


def read_rows(name):
    with (DATA / name).open() as stream:
        return list(csv.DictReader(stream))


def bm3(v, v0=74.71, k0=160.2, kp=3.99):
    x = (v0 / np.asarray(v)) ** (1 / 3)
    return 1.5 * k0 * (x**7 - x**5) * (1 + 0.75 * (kp - 4) * (x**2 - 1))


def gamma_theta(
    volume, gamma_law="eq11", theta_law="integrated", q0=1.65, q1=11.8, v0=74.71
):
    """Eq. 11 plus integrated theta, or explicitly named diagnostic alternatives.

    q1=0 is the constant-q model underlying Table 3 and Figure 6.
    'local_power' substitutes q(V) into gamma0*x**q; it is NOT Eq. 11.
    """
    v = np.asarray(volume, dtype=float)
    if np.any(~np.isfinite(v) | (v <= 0)):
        raise ValueError("Positive finite conventional-cell volumes required")
    u = np.log(v / v0)

    def gamma_at(logv):
        if gamma_law == "eq11":
            exponent = q0 * logv if abs(q1) < 1e-12 else q0 / q1 * np.expm1(q1 * logv)
        elif gamma_law == "local_power":
            exponent = q0 * np.exp(q1 * logv) * logv
        else:
            raise ValueError(gamma_law)
        return 1.524 * np.exp(exponent)

    gamma = gamma_at(u)
    q = q0 * np.exp(q1 * u)
    if theta_law == "integrated":
        z = u[..., None] * (NODES + 1) / 2
        exponent = -u / 2 * np.sum(WEIGHTS * gamma_at(z), axis=-1)
    elif theta_law == "eq4_local_q":
        exponent = np.divide(
            1.524 - gamma,
            q,
            out=np.array(-1.524 * u, copy=True),
            where=np.abs(q) > 1e-12,
        )
    elif theta_law == "eq4_q0":
        exponent = (1.524 - gamma) / q0 if abs(q0) > 1e-12 else -1.524 * u
    elif theta_law == "power":
        exponent = -gamma * u
    elif theta_law == "fixed":
        exponent = np.zeros_like(u)
    else:
        raise ValueError(theta_law)
    return gamma, 773 * np.exp(exponent)


def energy(theta, temperature, atoms=2):
    theta, t = np.broadcast_arrays(theta, temperature)
    y = theta / t
    z = y[..., None] * (NODES + 1) / 2
    integral = y / 2 * np.sum(WEIGHTS * z**3 / np.expm1(z), axis=-1)
    return 9 * atoms * R * t * integral / y**3


def pressure(
    volume,
    temperature,
    gamma_law="eq11",
    theta_law="integrated",
    tr=300,
    prefactor=1,
    atoms=2,
    formula_units=4,
    q0=1.65,
    q1=11.8,
    v0=74.71,
):
    v, t = np.broadcast_arrays(
        np.asarray(volume, dtype=float), np.asarray(temperature, dtype=float)
    )
    if np.any(~np.isfinite(t) | (t <= 0)) or not np.isfinite(tr) or tr <= 0:
        raise ValueError("Positive finite temperatures required")
    gamma, theta = gamma_theta(v, gamma_law, theta_law, q0, q1, v0)
    delta = prefactor * (energy(theta, t, atoms) - energy(theta, tr, atoms))
    return bm3(v, v0=v0) + gamma * delta / (v * N_A * 1e-30 / formula_units) / 1e9


def adaptive_pressure(volume, temperature):
    """Independent scalar adaptive integration of both thermodynamic integrals."""
    u = np.log(volume / 74.71)
    gamma = 1.524 * np.exp(1.65 / 11.8 * (np.exp(11.8 * u) - 1))
    integral = quad(
        lambda z: 1.524 * np.exp(1.65 / 11.8 * (np.exp(11.8 * z) - 1)),
        0,
        u,
        epsabs=1e-12,
        epsrel=1e-12,
    )[0]
    theta = 773 * np.exp(-integral)

    def e(t):
        y = theta / t
        return (
            18
            * R
            * t
            / y**3
            * quad(lambda z: z**3 / np.expm1(z), 0, y, epsabs=1e-12, epsrel=1e-12)[0]
        )

    f = ((74.71 / volume) ** (2 / 3) - 1) / 2
    cold = 3 * 160.2 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (3.99 - 4) * f)
    return cold + gamma * (e(temperature) - e(300)) * 4 / (volume * N_A * 1e-30) / 1e9


def metrics(residual):
    r = np.asarray(residual)
    return dict(
        rows=len(r),
        rmse_gpa=float(np.sqrt(np.mean(r * r))),
        mean_residual_gpa=float(np.mean(r)),
        max_abs_residual_gpa=float(np.max(np.abs(r))),
    )


def subset_fit(rows, volume_key, ambient=False):
    rows = [r for r in rows if float(r["temperature_k"]) > 300]
    v = np.array([float(r[volume_key]) for r in rows])
    t = np.array([float(r["temperature_k"]) for r in rows])
    p = (
        np.zeros(len(v))
        if ambient
        else np.array([float(r["pressure_gpa"]) for r in rows])
    )
    fit = least_squares(
        lambda z: pressure(v, t, q0=z[0], q1=0) - p,
        [1.3],
        bounds=(0.01, 4),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    residual = pressure(v, t, q0=fit.x[0], q1=0) - p
    out = dict(
        q=float(fit.x[0]),
        free_parameters=["q"],
        objective="equal pressure residuals",
        fixed_parameters={
            k: val for k, val in PARAMETERS.items() if k not in ("q0", "q1")
        },
        constant_q_fit=metrics(residual),
        source_exact_fit_reproduced=False,
    )
    if not ambient:
        pv = pressure(v, t)
        predicted_v = np.array(
            [brentq(lambda vv: pressure(vv, tt) - pp, 45, 80) for tt, pp in zip(t, p)]
        )
        out.update(
            eq11=metrics(pv - p),
            eq11_mean_abs_relative_volume_difference_percent=float(
                np.mean(np.abs(predicted_v / v - 1)) * 100
            ),
        )
    return out


def reproduce():
    from scripts.reproduce_fei_2007_gold import mgo_pressure, observations

    cold = read_rows("mgo-speziale-2001-table1-compression.csv")
    v = np.array([float(r["unit_cell_volume_a3"]) for r in cold])
    p = np.array([float(r["pressure_gpa"]) for r in cold])
    fit = least_squares(
        lambda z: bm3(v, z[0], 160.2, z[1]) - p,
        [74.71, 3.99],
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    rows, d = observations()
    cases = {
        "eq11_integrated_300K": {},
        "eq11_integrated_298K": dict(tr=298),
        "eq11_integrated_298.15K": dict(tr=298.15),
        "eq11_eq4_qV": dict(theta_law="eq4_local_q"),
        "eq11_eq4_q0": dict(theta_law="eq4_q0"),
        "eq11_power_theta": dict(theta_law="power"),
        "eq11_fixed_theta": dict(theta_law="fixed"),
        "constant_q_1.65": dict(q1=0),
        "local_power_integrated_300K": dict(gamma_law="local_power"),
        "local_power_integrated_298K": dict(gamma_law="local_power", tr=298),
        "local_power_integrated_298.15K": dict(gamma_law="local_power", tr=298.15),
        "local_power_eq4_qV_300K": dict(
            gamma_law="local_power", theta_law="eq4_local_q"
        ),
        "local_power_power_theta_300K": dict(
            gamma_law="local_power", theta_law="power"
        ),
        "eq11_n1_wrong_normalization": dict(atoms=1),
        "eq11_Z1_wrong_normalization": dict(formula_units=1),
        "eq11_n8_Z1_cell_normalization": dict(atoms=8, formula_units=1),
        "eq11_R_8.31451": dict(prefactor=8.31451 / R),
        "eq11_51.041_heat_capacity": dict(prefactor=51.041 / (6 * R)),
        "eq11_literal_Table3_specific_heat": dict(prefactor=0.12664 * 40.304 / (6 * R)),
        "eq11_table1_ambient_a0": dict(v0=4.2118**3),
    }
    predictions = {k: pressure(d["vm"], d["t"], **kw) for k, kw in cases.items()}
    comparison = {
        k: dict(
            Fei2004=metrics(pred[:26] - d["p"][:26]),
            Hirose2006=metrics(pred[26:] - d["p"][26:]),
            Hirose2006_pressures_gpa=pred[26:].tolist(),
        )
        for k, pred in predictions.items()
    }
    source_prediction = predictions["eq11_integrated_300K"]
    # Deterministic half-last-digit lattice/pressure limits, not experimental sigmas.
    dp_round = (
        np.abs(
            pressure((d["vm"] ** (1 / 3) - 0.00005) ** 3, d["t"]) - source_prediction
        )
        + np.r_[np.full(26, 0.005), [0.05, 0.05]]
    )
    rounding = dict(
        pressure_plus_lattice_half_digit_bounds_gpa=dp_round.tolist(),
        original_eq11_outside_bounds=int(
            np.sum(np.abs(source_prediction - d["p"]) > dp_round)
        ),
        temperatures_held_at_reported_values=True,
    )
    grid = [
        (74.71 * x, t)
        for x in [0.60, 0.64, 0.70, 0.80, 0.90, 1.0, 1.02]
        for t in [300.0, 500.0, 1100.0, 2000.0, 3663.0]
    ]
    adaptive_difference = max(
        abs(pressure(vv, tt) - adaptive_pressure(vv, tt)) for vv, tt in grid
    )
    gamma, theta = gamma_theta(74.71 * np.array([1.0, 0.85, 0.64, 0.60]))
    shock = read_rows("mgo-svendsen-1987-table6.csv")
    rho = np.array([float(r["initial_density_mg_m3"]) for r in shock])
    us = np.array([float(r["model_shock_velocity_km_s"]) for r in shock])
    ps = np.array([float(r["model_pressure_gpa"]) for r in shock])
    ts = np.array([float(r["greybody_temperature_k"]) for r in shock])
    relative_v = 1 - ps / (rho * us**2)
    measured_v = 4 * 0.040304 / (rho * 1000 * N_A) / 1e-30 * relative_v
    shock_results = {}
    for label, sv in [
        ("per_shot_density_molar_mass_volume", measured_v),
        ("per_shot_ratio_times_adopted_V0", 74.71 * relative_v),
    ]:
        pred_v = np.array(
            [brentq(lambda vv: pressure(vv, tt) - pp, 40, 75) for tt, pp in zip(ts, ps)]
        )
        shock_results[label] = dict(
            cell_volume_a3=sv.tolist(),
            pressure_residual_gpa=(pressure(sv, ts) - ps).tolist(),
            predicted_volume_difference_percent=(100 * (pred_v / sv - 1)).tolist(),
        )
    report = dict(
        format="peritheos.speziale-2001-mgo-validation",
        format_version=1,
        parameters=PARAMETERS,
        atoms_per_formula_unit=2,
        formula_units_per_cell=4,
        printed_equations_validated=True,
        exact_author_thermal_fit_reproduced=False,
        exact_later_pressure_reducer_recovered=False,
        published_coefficients_changed=False,
        cold_table1=dict(
            rows=len(v),
            V0=float(fit.x[0]),
            K0=160.2,
            K0_prime=float(fit.x[1]),
            objective="Independent equal pressure residuals; original objective unknown",
        ),
        source_evidence=dict(
            Speziale2001=dict(
                url="https://duffy.princeton.edu/sites/g/files/toruqf616/files/speziale_et_al-2001-jgrse.pdf",
                pdf_sha256="c0431bdddf8f4bc8269a7e0d1b79d30d92a38b4a778d9edddad9cb4c45b1f94b",
                locations="Eqs. 1-4 and 10-11; Tables 1-3; Sections 4.2 and 5.1-5.4; Figures 6, 9, 10",
                theta_extension="Integrated differential identity; no explicit variable-q theta executable recovered",
            ),
            Fei2004=dict(
                doi="10.1016/j.pepi.2003.09.018",
                pdf_sha256="b3fb00483e702e48d9d4cb5b4fa1dcefbbe176f94b4eca438d84b99a225342fd",
                locations="Original Tables 1 and 3; Section 3.1",
                provenance="Cached primary text and pixels; inherited hash, not a newly fetched PDF",
            ),
            Dorogokupets2010=dict(
                doi="10.1007/s00269-010-0367-2",
                access="Author-uploaded full text on ResearchGate and publisher abstract",
                use="Primary reanalysis corroborates multiple differing historical MgO reductions; does not establish the inferred local-power reducer",
            ),
        ),
        thermodynamic_checkpoints=dict(
            volume_ratios=[1.0, 0.85, 0.64, 0.60],
            gamma=gamma.tolist(),
            theta_k=theta.tolist(),
            gamma_infinite_compression=float(1.524 * np.exp(-1.65 / 11.8)),
        ),
        independent_quadrature=dict(
            states=len(grid), max_difference_gpa=float(adaptive_difference)
        ),
        existing_Au_diagnostic_max_difference_gpa=float(
            np.max(np.abs(source_prediction - mgo_pressure(d["vm"], d["t"])))
        ),
        Dewaele2000_hot41=subset_fit(
            read_rows("mgo-dewaele-2000-table2-pvt.csv"), "mgo_unit_cell_volume_a3"
        ),
        Fiquet1999_COD_hot34=subset_fit(
            read_rows("mgo-fiquet-1999-cod-thermal-expansion.csv"),
            "conventional_cell_volume_a3",
            True,
        ),
        Svendsen1987_shock4=shock_results,
        calibration_conventions=comparison,
        rounding=rounding,
        local_power_is_an_inference=True,
        local_power_q_identity="dln(gamma)/dln(V)=q(V)*(1+q1*ln(V/V0)), not q(V)",
        missing_evidence=[
            "Original author pressure-reduction executable and unrounded inputs",
            "Full selected Fei 1999 and Utsumi 1998 thermal regression inputs",
            "Original thermal objective, weights and covariance",
            "Explicit variable-q theta prescription used by Speziale and later reducers",
            "Exact shock volume normalization used in the source reduction",
        ],
        input_sha256={
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                "mgo-speziale-2001-table1-compression.csv",
                "mgo-dewaele-2000-table2-pvt.csv",
                "mgo-fiquet-1999-cod-thermal-expansion.csv",
                "mgo-svendsen-1987-table6.csv",
                "gold-fei-2004-table1.csv",
                "gold-hirose-2006-table1.csv",
            ]
        },
        source_code_sha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
    )
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(
        [
            "source",
            "row",
            "mgo_volume_a3",
            "temperature_k",
            "reported_pressure_gpa",
            *predictions,
        ]
    )
    for i, row in enumerate(rows):
        writer.writerow(
            [
                row["source"],
                row["row"],
                d["vm"][i],
                d["t"][i],
                d["p"][i],
                *[p[i] for p in predictions.values()],
            ]
        )
    return report, stream.getvalue()


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report, table = reproduce()
    if args.check:
        from scripts.reproduce_fei_2007_gold import check_saved

        check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), report)
        from scripts.check_numerical_archive import check_csv

        calculated = set(table.splitlines()[0].split(",")) - {
            "source",
            "row",
            "mgo_volume_a3",
            "temperature_k",
            "reported_pressure_gpa",
        }
        check_csv(RESIDUALS.read_text(encoding="utf-8"), table, calculated, atol=1e-6)
        print("Speziale MgO validation archive reproduced")
    else:
        OUTPUT.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        RESIDUALS.write_text(table, encoding="utf-8")
        print(OUTPUT.relative_to(ROOT))


if __name__ == "__main__":
    main()
