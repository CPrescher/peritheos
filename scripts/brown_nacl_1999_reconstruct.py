#!/usr/bin/env python3
"""Documented primary-data reconstruction, NOT Brown's recovered author fit.

Brown's calculated Tables I/II enter only comparison(), after fitting.
Units: v cm3/g; pressure GPa; energy J/g; T K. 1 GPa cm3/g = 1000 J/g.
No dependency on the Peritheos runtime or the ferrous-ringwoodite scripts.
"""

from __future__ import annotations

import csv
import json
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, solve_ivp
from scipy.interpolate import BSpline, CubicSpline, make_smoothing_spline
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUT = ROOT / "output/material-candidates/brown-nacl-1999"
PREFIX = "brown-nacl-1999-"
R = 8.31446261815324
M = 58.443
N = 2
VREF = 1 / 2.1645  # Fritz Table2, initial 293 K shock state
TREF = 293.0
PREF = 0.000101325  # GPa, 1 atmosphere


def read(name):
    with (DATA / (PREFIX + name + ".csv")).open() as f:
        return list(csv.DictReader(f))


def mm(a, b):
    """Unoptimized explicit contraction avoids spurious Accelerate BLAS flags."""
    return np.einsum("ij,j...->i...", a, b, optimize=False)


def energy(T, theta):
    """Thermal (zero-point excluded) Debye energy J/g NaCl, two atoms/formula."""
    T, theta = np.broadcast_arrays(
        np.asarray(T, dtype=float), np.asarray(theta, dtype=float)
    )
    # Gauss integration on [0,1]: E=9nRT int_0^1 z^3*y/expm1(y*z) dz / M.
    z, w = np.polynomial.legendre.leggauss(32)
    z = (z + 1) / 2
    w = w / 2
    y = theta / np.maximum(T, 1e-50)
    a = np.minimum(y[..., None] * z, 700)
    result = 9 * N * R * T / M * np.sum(w * z**3 * y[..., None] / np.expm1(a), axis=-1)
    return np.where(T > 0, result, 0.0)


def cv(T, theta):
    e = energy(T, theta)
    return 4 * e / np.asarray(T) - 9 * N * R / M * (np.asarray(theta) / T) / np.expm1(
        np.asarray(theta) / T
    )


LOW = np.array([3.836e-5, 2.979e-8, 1.167e-11])
HIGH = np.array([6.594e-5, -7.499e-8, 1.100e-10])


def length_integral(t):
    def prim(c, x):
        return sum(a * x ** (i + 1) / (i + 1) for i, a in enumerate(c))

    def val(x):
        return prim(LOW, min(x, 550)) + (
            prim(HIGH, x) - prim(HIGH, 550) if x > 550 else 0.0
        )

    return val(t) - val(25.0)


def expansion_volume(T, convention="room_length"):
    delta = np.array([length_integral(float(t) - 273.15) for t in np.atleast_1d(T)])
    initial = length_integral(TREF - 273.15)
    return VREF * (
        ((1 + delta) / (1 + initial)) ** 3
        if convention == "room_length"
        else np.exp(3 * (delta - initial))
    )


@dataclass
class Settings:
    q: float = 1.45
    gamma0: float = 1.63  # rounded Spetzler Table4 gamma(300 K)=1.6275;
    # also consistent with Yamamoto's 1-bar 1.64 +/- .03, no Brown-table fit
    theta_ref: float = 276.0  # Fritz Table2 at initial shock state
    join_width: float = 0.0  # baseline: sharp compressed-power / expanded-constant join
    smoothing: float = 0.0  # unweighted finite-dimensional least-squares baseline
    n_basis: int = 18
    junction_weight: float = 1.0
    shock_mode: str = "observations"
    energy_reference: str = "initial_total_zero"
    expansion_convention: str = "room_length"
    expansion_step_c: float = 50.0
    shock_initial_density: float = 2.1645
    spline_kind: str = "finite_basis"
    reference_anchor: bool = False  # exact measured initial P,V,T in independent refit


class Thermal:
    def __init__(self, s):
        self.s = s
        vv = np.linspace(0.30, 0.55, 8001)
        lg = cumulative_trapezoid(self.gamma(vv) / vv, vv, initial=0)
        integ = CubicSpline(vv, lg)
        self.ltheta = CubicSpline(vv, np.log(s.theta_ref) - (integ(vv) - integ(VREF)))

    def gamma(self, v):
        x = np.log(np.asarray(v) / VREF)
        if self.s.join_width:
            w = self.s.join_width
            x = -w * np.logaddexp(0, -x / w)
        else:
            x = np.minimum(x, 0)
        return self.s.gamma0 * np.exp(self.s.q * x)

    def theta(self, v):
        return np.exp(self.ltheta(v))

    def pressure(self, v, T):
        return self.gamma(v) * energy(T, self.theta(v)) / (1000 * np.asarray(v))


def inputs(s):
    bk = read("boehler-kennedy-1980-table1")
    v25 = expansion_volume([298.15], s.expansion_convention)[0]
    bv = np.array([float(r["volume_ratio_25c_1atm"]) * v25 for r in bk])
    bt = np.array([float(r["temperature_c"]) + 273.15 for r in bk])
    bp = np.array([float(r["pressure_kbar"]) / 10 for r in bk])
    # Polynomial-derived ambient anchors; NOT 77 original observations.
    et = np.arange(50, 751, s.expansion_step_c) + 273.15
    ev = expansion_volume(et, s.expansion_convention)
    rows = read("fritz-1971-table1")
    us = np.array([float(r["shock_velocity_km_s"]) for r in rows])
    up = np.array([float(r["particle_velocity_km_s"]) for r in rows])
    if s.shock_mode != "observations":
        up = np.linspace(0.33, 2.25, 43)
        if s.shock_mode == "linear":
            us = 3.528 + 1.343 * up
        elif s.shock_mode == "quadratic":
            us = 3.403 + 1.5422 * up - 0.07345 * up**2
        elif s.shock_mode == "piecewise":
            uq = 3.403 + 1.5422 * up - 0.07345 * up**2
            # Explicit reconstruction choice: switch at predicted 23 GPa.
            us = np.where(
                s.shock_initial_density * uq * up <= 23, uq, 3.528 + 1.343 * up
            )
        else:
            raise ValueError(s.shock_mode)
    hv = (1 - up / us) / s.shock_initial_density
    hp = s.shock_initial_density * us * up + PREF
    return bv, bt, bp, ev, et, hv, hp


class Model:
    def __init__(self, s):
        self.settings = s
        self.thermal = Thermal(s)
        bv, bt, bp, ev, et, hv, hp = inputs(s)
        v = np.r_[bv, ev, hv]
        self.v = v
        self.groups = np.array(
            ["BK"] * len(bv) + ["ED_polynomial"] * len(ev) + ["shock"] * len(hv)
        )
        self.targets_base = np.r_[
            bp - self.thermal.pressure(bv, bt),
            np.full(len(ev), PREF) - self.thermal.pressure(ev, et),
            hp,
        ]
        self.hv = hv
        self.hp = hp
        self.bv = bv
        self.bt = bt
        self.bp = bp
        x = v / VREF - 1
        xmin = min(x)
        xmax = max(x)
        # Clamped cubic B-spline basis; baseline has no shape/end-derivative constraints.
        interior = np.linspace(xmin, xmax, s.n_basis - 2)[1:-1]
        knots = np.r_[[xmin] * 4, interior, [xmax] * 4]
        n = len(knots) - 4
        if s.spline_kind == "natural_smoothing":
            unique, ix = np.unique(x, return_inverse=True)
            grouped = np.bincount(ix, minlength=len(unique)).astype(float)
            sample = make_smoothing_spline(
                unique, np.zeros(len(unique)), w=grouped, lam=s.smoothing
            )
            knots = sample.t
            n = len(sample.c)
        basis = BSpline(knots, np.eye(n), 3)
        A = basis(x)
        integral = basis.antiderivative()
        # Integral of P0 from vref to each shock v, times -1000 J/g.
        H = -1000 * VREF * (integral(hv / VREF - 1) - integral(0))
        weights = np.ones(len(v))
        if s.junction_weight != 1:
            # Authored sensitivity only; Brown's exact masks/weights unknown.
            mask = ((v > 0.445) & (v < 0.477)) | ((v > 0.412) & (v < 0.427))
            weights[mask] = s.junction_weight
        self.weights = weights
        gx, gw = np.polynomial.legendre.leggauss(180)
        grid = xmin + (gx + 1) / 2 * (xmax - xmin)
        second = basis(grid, nu=2)
        penalty = mm(second.T * (gw * (xmax - xmin) / 2), second)
        if s.spline_kind == "natural_smoothing":
            grouped = np.bincount(ix, weights=weights, minlength=len(unique))
            aggregation = np.zeros((len(unique), len(x)))
            aggregation[ix, np.arange(len(x))] = weights / grouped[ix]
            operator_unique = np.column_stack(
                [
                    make_smoothing_spline(unique, row, w=grouped, lam=s.smoothing).c
                    for row in np.eye(len(unique))
                ]
            )
            operator = mm(operator_unique, aggregation)
        else:
            normal = mm(A.T, weights[:, None] * A) + s.smoothing * penalty
            operator = np.linalg.solve(normal, A.T * weights)
        offset = np.zeros(n)
        if s.reference_anchor:
            if s.spline_kind != "finite_basis":
                raise ValueError("Reference constraint requires finite_basis")
            b = basis(0.0)
            u = np.linalg.solve(normal, b)
            projection = u / np.sum(b * u)
            operator -= projection[:, None] * mm(b[None, :], operator)
            offset = projection * (PREF - float(self.thermal.pressure(VREF, TREF)))
        eh = 0.5 * (hp + PREF) * (1 / s.shock_initial_density - hv) * 1000
        eref = (
            -float(energy(TREF, s.theta_ref))
            if s.energy_reference == "initial_total_zero"
            else 0.0
        )
        c = mm(operator, self.targets_base) + offset
        history = []
        for i in range(200):
            cold = eref + mm(H, c)
            target = self.targets_base.copy()
            target[-len(hv) :] = hp - self.thermal.gamma(hv) * (eh - cold) / (1000 * hv)
            new = mm(operator, target) + offset
            change = float(np.max(abs(mm(A, new - c))))
            history.append(change)
            c = new
            if change < 1e-10:
                break
        if history[-1] >= 1e-10:
            raise RuntimeError("Cold/Hugoniot iteration failed")
        self.spline = BSpline(knots, c, 3)
        self.target = target
        self.residual = mm(A, c) - target
        self.history = history
        self.basis_knots = knots
        self.coefficients = c
        # Verify fixed point independently via the simultaneous linear system.
        L = A.copy()
        L[-len(hv) :] -= self.thermal.gamma(hv)[:, None] * H / (1000 * hv[:, None])
        y = self.targets_base.copy()
        y[-len(hv) :] = hp - self.thermal.gamma(hv) * (eh - eref) / (1000 * hv)
        # Fixed-point smoothing is a reduced-data smoother, not minimizing coupled L residuals.
        self.rh_residual = (
            self.p0(hv)
            - self.thermal.gamma(hv) * mm(H, c) / (1000 * hv)
            - y[-len(hv) :]
        )
        self.cold_energy = eref + mm(H, c)
        self.hugoniot_energy = eh
        self.e_thermal_h = eh - self.cold_energy

    def p0(self, v):
        return self.spline(np.asarray(v) / VREF - 1)

    def pressure(self, v, T):
        return self.p0(v) + self.thermal.pressure(v, T)

    def kt(self, v, T):
        v = np.asarray(v)
        h = 2e-6
        return -v * (self.pressure(v + h, T) - self.pressure(v - h, T)) / (2 * h)

    def volume(self, p, T):
        lo = max(0.305, self.v.min())
        hi = min(0.55, self.v.max())
        return brentq(lambda v: float(self.pressure(v, T)) - p, lo, hi)

    def properties(self, p, T):
        v = self.volume(p, T)
        theta = float(self.thermal.theta(v))
        gamma = float(self.thermal.gamma(v))
        c = float(cv(T, theta))
        kt = float(self.kt(v, T))
        pt = gamma * c / (1000 * v)
        alpha = pt / kt
        cp = c + 1000 * alpha**2 * kt * v * T
        ks = kt * cp / c
        return dict(
            v=v,
            theta=theta,
            gamma=gamma,
            kt=kt,
            ks=ks,
            cv=c,
            cp=cp,
            alpha=alpha,
            bulk_sound_km_s=float(np.sqrt(ks * v)),
            dt_dp_s=T * gamma / ks,
        )

    def stats(self):
        by = {
            g: float(np.sqrt(np.mean(self.residual[self.groups == g] ** 2)))
            for g in np.unique(self.groups)
        }
        return dict(
            settings=asdict(self.settings),
            rows_by_group={
                g: int(sum(self.groups == g)) for g in np.unique(self.groups)
            },
            reduced_pressure_rms_gpa=by,
            iterations=len(self.history),
            last_iteration_change_gpa=self.history[-1],
            min_cold_kt_gpa=float(
                np.min(
                    -np.linspace(self.v.min(), self.v.max(), 501)
                    / VREF
                    * self.spline(
                        np.linspace(self.v.min(), self.v.max(), 501) / VREF - 1, nu=1
                    )
                )
            ),
            rh_residual_equals_fit_residual_max_gpa=float(
                np.max(abs(self.rh_residual - self.residual[-len(self.hv) :]))
            ),
            ambient_300=self.properties(PREF, 300),
            knots_x=self.basis_knots.tolist(),
            coefficients_gpa=self.coefficients.tolist(),
        )


def comparison(model):
    t1 = read("brown-table1-comparison")
    t2 = read("brown-table2-comparison")
    v = np.array([float(r["specific_volume_cm3_g"]) for r in t1])
    p = np.array([float(r["static_pressure_gpa"]) for r in t1])
    b = {
        c: CubicSpline(v, [float(r[c]) for r in t1])
        for c in ["static_pressure_gpa", "gamma", "theta_k"]
    }
    results = []
    full_standard = []
    full_zp = []
    for r in t2:
        vv = float(r["specific_volume_cm3_g"])
        for T in range(300, 1201, 100):
            if not r[f"pressure_{T}k_gpa"]:
                continue
            printed = float(r[f"pressure_{T}k_gpa"])
            # Table arithmetic diagnostic repeats the earlier explicit M=58.44
            # convention; primary reconstruction uses Fritz's M=58.443.
            standard = float(
                b["static_pressure_gpa"](vv)
                + b["gamma"](vv) * energy(T, b["theta_k"](vv)) * M / 58.44 / (1000 * vv)
            )
            zp3 = float(
                b["gamma"](vv)
                * (9 * N * R * b["theta_k"](vv) / 8 / 58.44 / 3)
                / (1000 * vv)
            )
            full_standard.append(standard - printed)
            full_zp.append(standard + zp3 - printed)
            if not model.v.min() <= vv <= model.v.max():
                continue
            results.append(
                [vv, T, printed, float(model.pressure(vv, T)), standard, standard + zp3]
            )
    arr = np.array(results)
    mask = (v >= model.v.min()) & (v <= model.v.max())
    x = v[mask] / VREF - 1
    d1 = model.spline(x, nu=1)
    d2 = model.spline(x, nu=2)
    kt0 = -v[mask] / VREF * d1
    kp0 = -1 - v[mask] / VREF * d2 / d1
    return dict(
        table1_pressure_rms_gpa=float(
            np.sqrt(np.mean((model.p0(v[mask]) - p[mask]) ** 2))
        ),
        table1_kt_rms_gpa=float(
            np.sqrt(
                np.mean((kt0 - np.array([float(r["kt_gpa"]) for r in t1])[mask]) ** 2)
            )
        ),
        table1_kprime_rms=float(
            np.sqrt(
                np.mean((kp0 - np.array([float(r["kprime"]) for r in t1])[mask]) ** 2)
            )
        ),
        table1_gamma_rms=float(
            np.sqrt(
                np.mean(
                    (
                        model.thermal.gamma(v[mask])
                        - np.array([float(r["gamma"]) for r in t1])[mask]
                    )
                    ** 2
                )
            )
        ),
        table1_theta_rms_k=float(
            np.sqrt(
                np.mean(
                    (
                        model.thermal.theta(v[mask])
                        - np.array([float(r["theta_k"]) for r in t1])[mask]
                    )
                    ** 2
                )
            )
        ),
        table2_cells_in_fitted_volume_range=len(results),
        table2_rms_gpa=float(np.sqrt(np.mean((arr[:, 3] - arr[:, 2]) ** 2))),
        table_arithmetic_cells=len(full_standard),
        table1_standard_to_table2_rms_gpa=float(
            np.sqrt(np.mean(np.array(full_standard) ** 2))
        ),
        table1_zp3_signature_to_table2_rms_gpa=float(
            np.sqrt(np.mean(np.array(full_zp) ** 2))
        ),
    ), results


def thermal_mode_factor_check():
    """Temperature increments cancel cold pressure and any T-independent offset.

    This tests a uniform energy multiplier only, not Brown's unspecified buggy
    component or a reconstruction of his original code.
    """
    t1 = read("brown-table1-comparison")
    v = np.array([float(r["specific_volume_cm3_g"]) for r in t1])
    gamma = CubicSpline(v, [float(r["gamma"]) for r in t1])
    theta = CubicSpline(v, [float(r["theta_k"]) for r in t1])
    expected, printed = [], []
    for row in read("brown-table2-comparison"):
        vv = float(row["specific_volume_cm3_g"])
        values = [
            (T, float(row[f"pressure_{T}k_gpa"]))
            for T in range(300, 1201, 100)
            if row[f"pressure_{T}k_gpa"]
        ]
        first_t, first_p = values[0]
        for T, P in values[1:]:
            increment = (
                gamma(vv)
                * (energy(T, theta(vv)) - energy(first_t, theta(vv)))
                * M
                / 58.44
                / (1000 * vv)
            )
            expected.append(float(increment))
            printed.append(P - first_p)
    expected, printed = np.array(expected), np.array(printed)
    scale = float(np.sum(expected * printed) / np.sum(expected**2))
    return {
        "classification": "calculated-table temperature-increment diagnostic; not independent experimental atom-count measurement",
        "increments": len(expected),
        "atom_count_in_standard_formula": N,
        "best_uniform_energy_scale": scale,
        "equivalent_atom_count_conditional_on_printed_gamma_theta": N * scale,
        "increment_rms_gpa_by_scale": {
            label: float(np.sqrt(np.mean((factor * expected - printed) ** 2)))
            for label, factor in [("one_third", 1 / 3), ("standard", 1), ("three", 3)]
        },
        "interpretation": "A uniform factor-three change to the entire thermal energy is inconsistent with the printed temperature increments; the author's omitted factor must be mapped to its actual component/reference treatment before any correction is implemented.",
    }


def high_pressure_shock_check(model):
    """Isolate extrapolated linear/quadratic shock curvature at a common anchor."""
    anchor = 1 / 3.023  # Fritz calculated Table4's 20 GPa Hugoniot density
    integ = model.spline.antiderivative()
    initial = -float(energy(TREF, model.settings.theta_ref)) - 1000 * VREF * (
        integ(anchor / VREF - 1) - integ(0)
    )
    curves = {}

    def shock(v, kind):
        eta = 1 - v / VREF
        if kind == "linear":
            up = 3.528 * eta / (1 - 1.343 * eta)
            us = 3.528 + 1.343 * up
        else:
            up = brentq(
                lambda up: up / (3.403 + 1.5422 * up - 0.07345 * up**2) - eta, 0, 4
            )
            us = 3.403 + 1.5422 * up - 0.07345 * up**2
        ph = 2.1645 * us * up + PREF
        eh = 0.5 * (ph + PREF) * (VREF - v) * 1000
        return ph, eh

    for kind in ["linear", "quadratic"]:

        def rhs(v, e):
            ph, eh = shock(v, kind)
            return [-1000 * ph + model.thermal.gamma(v) / v * (eh - e[0])]

        sol = solve_ivp(
            rhs, [anchor, 0.29], [initial], dense_output=True, rtol=1e-10, atol=1e-9
        )
        assert sol.success

        def press(v, T=300, sol=sol, kind=kind):
            ph, eh = shock(v, kind)
            cold = ph - model.thermal.gamma(v) * (eh - sol.sol(v)[0]) / (1000 * v)
            return float(cold + model.thermal.pressure(v, T))

        curves[kind] = press
    v30 = brentq(lambda v: curves["linear"](v) - 30, 0.29, anchor)
    return {
        "classification": "extrapolation sensitivity using published Us-Up coefficients, not observations",
        "common_cold_energy_anchor_v_cm3_g": anchor,
        "common_cold_energy_J_g": initial,
        "linear_300k_30gpa_v_cm3_g": v30,
        "quadratic_minus_linear_gpa_at_that_volume": curves["quadratic"](v30) - 30,
        "range_note": "both require shock extrapolation beyond the 26.4 GPa measured maximum; common anchor is adopted",
    }


def main():
    OUT.mkdir(exist_ok=True, parents=True)
    cases = [
        ("unweighted_baseline", Settings()),
        ("regularized", Settings(smoothing=1e-6)),
        ("junction_downweight", Settings(smoothing=1e-6, junction_weight=0.1)),
        ("smooth_gamma_join", Settings(smoothing=1e-6, join_width=0.025)),
        (
            "natural_smoothing",
            Settings(smoothing=1e-7, spline_kind="natural_smoothing", join_width=0.025),
        ),
        (
            "published_shock_representations",
            Settings(smoothing=1e-6, shock_mode="piecewise", join_width=0.025),
        ),
        ("q1", Settings(smoothing=1e-6, q=1)),
        ("q2", Settings(smoothing=1e-6, q=2)),
        ("theta285", Settings(smoothing=1e-6, theta_ref=285)),
        (
            "missing_initial_thermal_energy",
            Settings(smoothing=1e-6, energy_reference="cold_at_initial_zero"),
        ),
        (
            "instantaneous_alpha_interpretation",
            Settings(smoothing=1e-6, expansion_convention="instantaneous"),
        ),
        ("shock_linear", Settings(smoothing=1e-6, shock_mode="linear")),
        ("shock_quadratic", Settings(smoothing=1e-6, shock_mode="quadratic")),
        ("shock_piecewise", Settings(smoothing=1e-6, shock_mode="piecewise")),
    ]
    report = {
        "classification": "Brown (1999) EOS cannot be reproduced from the available published material; independent historical-data reconstructions only",
        "constants": dict(
            R_j_mol_k=R,
            M_g_mol=M,
            atoms_per_formula=N,
            initial_v_cm3_g=VREF,
            initial_T_k=TREF,
            initial_p_gpa=PREF,
        ),
        "fit_comparison_separation": "Tables I/II are read only by comparison(), after fitting",
        "cases": {},
    }
    import platform

    import scipy

    report["environment"] = {
        "python": platform.python_version(),
        "numpy": np.__version__,
        "scipy": scipy.__version__,
    }
    correspondence = OUT / "sources/brown-author-reply-2026-10-16.json"
    if correspondence.exists():
        evidence = json.loads(correspondence.read_text())
        report["author_correspondence"] = {
            "source": str(correspondence.relative_to(ROOT)),
            "source_type": evidence["source_type"],
            "author_acknowledges_published_calculation_bug": True,
            "reported_cause": "Omitted factor of three in an energy term for three vibrational modes per ion, across versions of Debye code",
            "exact_affected_term_and_corrected_tables_supplied": False,
            "numerical_reconstruction_changed": False,
        }
    report["thermal_mode_factor_diagnostic"] = thermal_mode_factor_check()
    shock = read("fritz-1971-table1")
    up = np.array([float(a["particle_velocity_km_s"]) for a in shock])
    us = np.array([float(a["shock_velocity_km_s"]) for a in shock])
    c_linear = np.polynomial.polynomial.polyfit(up, us, 1)
    c_quadratic = np.r_[
        3.403, np.linalg.lstsq(np.column_stack([up, up**2]), us - 3.403, rcond=None)[0]
    ]
    report["velocity_fit_verification"] = {}
    for name, c in [
        ("unweighted_linear", c_linear),
        ("unweighted_quadratic_fixed_intercept", c_quadratic),
        ("published_linear", [3.528, 1.343]),
        ("published_quadratic", [3.403, 1.5422, -0.07345]),
    ]:
        report["velocity_fit_verification"][name] = {
            "ascending_coefficients": list(c),
            "rms_us_km_s": float(
                np.sqrt(np.mean((np.polynomial.polynomial.polyval(up, c) - us) ** 2))
            ),
            "objective": "sum (Us_observed - Us(Up))^2, all weights one; Up conditionally exact; printed rounding retained",
        }
    models = {}
    for name, s in cases:
        m = Model(s)
        models[name] = m
        st = m.stats()
        st["comparison"], rows = comparison(m)
        report["cases"][name] = st
        with (OUT / (name + "-brown-comparison.csv")).open("w", newline="") as f:
            w = csv.writer(f)
            w.writerow(
                [
                    "v_cm3_g",
                    "T_k",
                    "Brown_table2_GPa",
                    "reconstructed_GPa",
                    "Brown_table1_standard_GPa",
                    "Brown_table1_plus_zp3_GPa",
                ]
            )
            w.writerows(rows)
    # Conditional q profile: refit cold spline for each q, score BK only.
    qscan = []
    for q in np.linspace(0.7, 2.3, 33):
        m = Model(Settings(smoothing=1e-6, q=float(q)))
        qscan.append([float(q), m.stats()["reduced_pressure_rms_gpa"]["BK"]])
    report["conditional_q_scan"] = qscan
    report["conditional_q_minimum"] = min(qscan, key=lambda a: a[1])
    report["conditional_q_profiles"] = {}
    for kind in ["natural_smoothing", "published_shock_representations"]:
        s0 = models[kind].settings
        from dataclasses import replace

        scan = []
        for q in np.linspace(1, 2, 21):
            m = Model(replace(s0, q=float(q)))
            scan.append([float(q), m.stats()["reduced_pressure_rms_gpa"]["BK"]])
        report["conditional_q_profiles"][kind] = {
            "scan": scan,
            "minimum": min(scan, key=lambda a: a[1]),
        }
    # Smoothing sensitivity with identical data and parameters, not tuned to Brown.
    report["smoothing_scan"] = [
        Model(Settings(smoothing=smoothing)).stats()
        for smoothing in [0, 1e-9, 1e-8, 1e-7, 1e-6, 1e-5, 1e-4]
    ]
    ref = models["natural_smoothing"]
    report["high_pressure_linear_quadratic_check"] = high_pressure_shock_check(ref)
    regularrows = []
    for i, (v, y, r, g, w) in enumerate(
        zip(ref.v, ref.target, ref.residual, ref.groups, ref.weights)
    ):
        regularrows.append([i, g, v, y, ref.p0(v), r, w])
    with (OUT / "natural-smoothing-reduced-data.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "row",
                "group",
                "v_cm3_g",
                "reduced_P0_GPa",
                "spline_P0_GPa",
                "residual_GPa",
                "weight",
            ]
        )
        w.writerows(regularrows)
    with (OUT / "natural-smoothing-hugoniot.csv").open("w", newline="") as f:
        w = csv.writer(f)
        w.writerow(
            [
                "v_cm3_g",
                "Ph_GPa",
                "Eh_J_g",
                "cold_energy_J_g",
                "thermal_energy_J_g",
                "calculated_Th_K",
            ]
        )
        for v, p, eh, ec, eth in zip(
            ref.hv, ref.hp, ref.hugoniot_energy, ref.cold_energy, ref.e_thermal_h
        ):
            th = (
                brentq(lambda T: float(energy(T, ref.thermal.theta(v))) - eth, 1, 6000)
                if eth > 0
                else float("nan")
            )
            w.writerow([v, p, eh, ec, eth, th])
    report["validation"] = {}
    if (DATA / (PREFIX + "yamamoto-1987-table2.csv")).exists():
        vals = []
        for r in read("yamamoto-1987-table2"):
            T = float(r["temperature_k"])
            prop = ref.properties(PREF, T)
            vals.append(
                dict(
                    T=T,
                    model=prop,
                    source={k: float(v) if v else None for k, v in r.items()},
                )
            )
        report["validation"]["yamamoto"] = vals
        # Sound speeds derived from the measured 294 K density, not the
        # inconsistent printed specific-volume column. No shear EOS is fitted.
        first = vals[0]["source"]
        rho_measured = 2.160
        report["validation"]["sound_294k"] = {
            "classification": "derived from resonance elastic constants and independently reported density",
            "source_density_g_cm3": rho_measured,
            "source_bulk_sound_km_s": float(
                np.sqrt(first["ks_kbar"] * 0.1 / rho_measured)
            ),
            "source_vl_100_km_s": float(
                np.sqrt(first["c11_kbar"] * 0.1 / rho_measured)
            ),
            "source_vs_100_km_s": float(
                np.sqrt(first["c44_kbar"] * 0.1 / rho_measured)
            ),
            "model_bulk_sound_km_s": ref.properties(PREF, 294)["bulk_sound_km_s"],
            "model_predicts_shear": False,
        }
    if (DATA / (PREFIX + "spetzler-1972-table4-ks.csv")).exists():
        report["validation"]["spetzler"] = [
            dict(
                source=r,
                model=ref.properties(
                    float(r["pressure_kbar"]) / 10 + PREF, float(r["temperature_k"])
                ),
            )
            for r in read("spetzler-1972-table4-ks")
        ]
    report["validation_model"] = "natural_smoothing"
    report["validation"]["janaf1971"] = []
    for r in read("janaf-1971-nacl-cr-cp"):
        T = float(r["temperature_k"])
        if 298 <= T <= 1000:
            props = ref.properties(PREF, T)
            source = float(r["cp_j_mol_k"]) / float(r["historical_molar_mass_g_mol"])
            report["validation"]["janaf1971"].append(
                {
                    "T_k": T,
                    "model_cp_j_g_k": props["cp"],
                    "source_cp_j_g_k": source,
                    "relative_difference": props["cp"] / source - 1,
                }
            )
    report["validation"]["boehler1981"] = []
    for r in read("boehler-1981-table1"):
        T = float(r["temperature_c"]) + 273.15
        for p in [1, 2, 3, 4, 5]:
            if p == 5 and T > 1000:
                continue  # 800 C data/fit figure stops near 45 kbar
            pred = ref.properties(p + PREF, T)["dt_dp_s"]
            measured = sum(float(r[f"a{i}"]) * (p * 10) ** i for i in range(4)) * 10
            report["validation"]["boehler1981"].append(
                dict(
                    p_gpa=p,
                    T_k=T,
                    model_k_per_gpa=pred,
                    published_polynomial_k_per_gpa=measured,
                    relative_difference=pred / measured - 1,
                )
            )
    # These are source accuracy statements mapped to a pressure scale, not
    # independent Gaussian sigmas. No author-fit confidence interval is implied.
    bounds = []
    for v, T, P, r in zip(ref.bv, ref.bt, ref.bp, read("boehler-kennedy-1980-table1")):
        ratio = float(r["volume_ratio_25c_1atm"])
        dv = 0.007 * abs(1 - ratio) * expansion_volume([298.15])[0]
        kt = float(ref.kt(v, T))
        boundP = 0.004 * P
        boundV = kt / v * dv
        boundT = (
            float(ref.thermal.gamma(v) * cv(T, ref.thermal.theta(v)) / (1000 * v)) * 0.3
        )
        bounds.append(
            dict(
                pressure_gpa=P,
                temperature_k=T,
                pressure_bound_gpa=boundP,
                volume_equivalent_bound_gpa=boundV,
                temperature_equivalent_bound_gpa=boundT,
                linear_sum_bound_gpa=boundP + boundV + boundT,
            )
        )
    report["bk_accuracy_propagation_not_sigmas"] = bounds
    report["checkpoints"] = {}
    g = 1.19
    theta = 375.0
    v = 0.3709
    report["pointwise_published_table_check"] = {
        "v_cm3_g": v,
        "P0_gpa": 8.57,
        "gamma": g,
        "theta_k": theta,
        "M_g_mol": 58.44,
        "standard_P300_gpa": 8.57
        + g * float(energy(300, theta)) * M / 58.44 / (1000 * v),
        "standard_P1200_gpa": 8.57
        + g * float(energy(1200, theta)) * M / 58.44 / (1000 * v),
        "printed_table2_P300_gpa": 9.19,
        "printed_table2_P1200_gpa": 11.60,
        "interpolation_required": False,
        "zp_over_3_adopted": False,
    }
    for name, m in models.items():
        report["checkpoints"][name] = {
            "v03709_p0_gpa": float(m.p0(0.3709)),
            "v03709_p300_gpa": float(m.pressure(0.3709, 300)),
            "v03709_p1200_gpa": float(m.pressure(0.3709, 1200)),
        }
    (OUT / "reconstruction.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    plot(models, report)
    print(
        json.dumps(
            {
                k: {"rms": v["reduced_pressure_rms_gpa"], "comparison": v["comparison"]}
                for k, v in report["cases"].items()
            },
            indent=2,
        )
    )


def plot(models, report):
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    plt.rcParams.update(
        {"font.size": 10, "axes.spines.top": False, "axes.spines.right": False}
    )
    ref = models["natural_smoothing"]
    vv = np.linspace(ref.v.min(), ref.v.max(), 600)
    fig, axs = plt.subplots(2, 2, figsize=(12, 8), layout="constrained")
    ax = axs[0, 0]
    for g, label in [
        ("BK", "Boehler-Kennedy measurements"),
        ("ED_polynomial", "Enck-Dommel polynomial samples"),
        ("shock", "Fritz experimental reductions"),
    ]:
        mask = ref.groups == g
        ax.scatter(ref.v[mask], ref.target[mask], s=15, label=label)
    ax.plot(vv, ref.p0(vv), color="k", label="Reconstructed cold spline")
    ax.set(xlabel="Specific volume (cm³/g)", ylabel="Cold pressure (GPa)")
    ax.legend(fontsize=8)
    ax = axs[0, 1]
    for g in np.unique(ref.groups):
        mask = ref.groups == g
        ax.scatter(ref.v[mask], ref.residual[mask] * 1000, s=15, label=g)
    ax.axhline(0, color="k", lw=0.7)
    ax.set(xlabel="Specific volume (cm³/g)", ylabel="Reconstruction residual (MPa)")
    ax.legend(fontsize=8)
    ax = axs[1, 0]
    for name in [
        "unweighted_baseline",
        "regularized",
        "q1",
        "q2",
        "smooth_gamma_join",
        "published_shock_representations",
    ]:
        m = models[name]
        ax.plot(vv, m.pressure(vv, 300) - ref.pressure(vv, 300), label=name)
    ax.set(
        xlabel="Specific volume (cm³/g)",
        ylabel="300 K pressure difference from natural spline (GPa)",
    )
    ax.legend(fontsize=8)
    ax = axs[1, 1]
    r = read("brown-table1-comparison")
    v = np.array([float(a["specific_volume_cm3_g"]) for a in r])
    p = np.array([float(a["static_pressure_gpa"]) for a in r])
    mask = (v >= ref.v.min()) & (v <= ref.v.max())
    v = v[mask]
    p = p[mask]
    ax.plot(v, ref.p0(v) - p, "o-", label="Reconstructed cold - Brown Table I")
    zp3 = np.array(
        [
            float(a["gamma"])
            * 9
            * N
            * R
            * float(a["theta_k"])
            / 8
            / M
            / 3
            / (1000 * float(a["specific_volume_cm3_g"]))
            for a in r
        ]
    )[mask]
    ax.plot(v, zp3, "--", label="Published-parameter zero-point/3 signature")
    ax.set(xlabel="Specific volume (cm³/g)", ylabel="Pressure difference (GPa)")
    ax.legend(fontsize=8)
    fig.suptitle(
        "Brown NaCl: primary-data reconstruction and independent table comparisons",
        fontsize=13,
    )
    fig.savefig(OUT / "figures/reconstruction.png", dpi=180)
    plt.close(fig)
    fig, axs = plt.subplots(1, 2, figsize=(11, 4), layout="constrained")
    rows = read("fritz-1971-table1")
    u = np.array([float(a["particle_velocity_km_s"]) for a in rows])
    us = np.array([float(a["shock_velocity_km_s"]) for a in rows])
    grid = np.linspace(0, 2.6, 300)
    axs[0].scatter(u, us, s=18, label="Published impedance-match reductions")
    for c, label in [
        ([3.528, 1.343], "Published linear"),
        ([3.403, 1.5422, -0.07345], "Published quadratic"),
    ]:
        axs[0].plot(grid, np.polynomial.polynomial.polyval(grid, c), label=label)
    axs[0].axvspan(
        max(u), 2.6, color="gray", alpha=0.15, label="Beyond measured velocity range"
    )
    axs[0].set(xlabel="Particle velocity (km/s)", ylabel="Shock velocity (km/s)")
    axs[0].legend(fontsize=8)
    for name, a in report["conditional_q_profiles"].items():
        arr = np.array(a["scan"])
        axs[1].plot(arr[:, 0], arr[:, 1] * 1000, label=name.replace("_", " "))
    axs[1].axvline(1.45, color="gray", ls="--", label="Brown reported q ≈ 1.45")
    axs[1].set(
        xlabel="q (cold curve refitted at each value)",
        ylabel="Boehler-Kennedy RMS misfit (MPa)",
    )
    axs[1].legend(fontsize=8)
    fig.savefig(OUT / "figures/shock-and-q.png", dpi=180)
    plt.close(fig)
    if report["validation"].get("yamamoto"):
        vals = report["validation"]["yamamoto"]
        T = [r["T"] for r in vals]
        fig, axs = plt.subplots(1, 3, figsize=(12, 3.6), layout="constrained")
        for ax, key, source, label, scale in zip(
            axs,
            ["ks", "cp", "alpha"],
            ["ks_kbar", "cp_1e6_erg_g_k", "alpha_1e6_per_k"],
            [
                "Adiabatic bulk modulus (GPa)",
                "Cp (J/g/K)",
                "Volumetric expansivity (10⁻⁶/K)",
            ],
            [0.1, 0.1, 1],
        ):
            factor = 1e6 if key == "alpha" else 1
            ax.plot(T, [r["model"][key] * factor for r in vals], label="Reconstruction")
            ax.scatter(
                T,
                [r["source"][source] * scale for r in vals],
                s=20,
                label="Yamamoto Table2",
            )
            ax.set(xlabel="Temperature (K)", ylabel=label)
        axs[0].legend(fontsize=8)
        fig.suptitle("Physical validation (thermal expansion also entered calibration)")
        fig.savefig(OUT / "figures/validation.png", dpi=180)
        plt.close(fig)


if __name__ == "__main__":
    main()
