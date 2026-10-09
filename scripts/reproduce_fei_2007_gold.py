"""Au-only conditional replay of Fei (2007); never changes published EOS records.

python -m scripts.reproduce_fei_2007_gold
python -m scripts.reproduce_fei_2007_gold --check

Fit measured observations only. Calculated Figure 1 paths are separate checks.
No Pt imports/refits, material registration, manifest or ledger regeneration.
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
from scipy.optimize import least_squares

from peritheos import get_eos_record

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/fei-2007-gold-reproduction.json"
RESIDUALS = ROOT / "docs/data/fei-2007-gold-residuals.csv"
CURVES = ROOT / "docs/data/fei-2007-gold-source-curves.json"
NAMES = ("V0", "K0", "K0_prime", "gamma0", "q", "theta0")
PUBLISHED = np.array([67.85, 167.0, 6.0, 2.97, 0.6, 170.0])
NODES, WEIGHTS = np.polynomial.legendre.leggauss(64)


def energy(theta, temperature, atoms=1, prefactor=1.0):
    theta, temperature = np.broadcast_arrays(theta, temperature)
    y = theta / temperature
    z = y[..., None] * (NODES + 1) / 2
    integral = y / 2 * np.sum(WEIGHTS * z**3 / np.expm1(z), axis=-1)
    return prefactor * 9 * atoms * R * temperature * integral / y**3


def pressure(volume, temperature, coefficients=PUBLISHED, law="printed", prefactor=1):
    """Independent four-atom-cell Vinet plus SI Debye pressure evaluation."""
    v, t = np.broadcast_arrays(np.asarray(volume), np.asarray(temperature))
    v0, k0, kp, g0, q, th0 = coefficients
    logv = np.log(v / v0)
    x = np.exp(logv / 3)
    gamma = g0 * np.exp(q * logv)
    if law == "printed":
        exponent = -gamma * logv
    elif law == "integrated":
        exponent = -g0 * (logv if abs(q) < 1e-10 else np.expm1(q * logv) / q)
    else:
        raise ValueError(law)
    theta = th0 * np.exp(exponent)
    cold = 3 * k0 * (1 - x) * np.exp(1.5 * (kp - 1) * (1 - x)) / x**2
    delta = energy(theta, t, prefactor=prefactor) - energy(
        theta, 300, prefactor=prefactor
    )
    return cold + gamma * delta / (v * N_A * 1e-30 / 4) / 1e9


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
    p = dict(zip(NAMES, PUBLISHED))
    p["Tr"] = 300.0
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


def mgo_pressure(volume, temperature):
    """Speziale variable-q reconstruction with integrated dln(theta)=-gamma dln(V).

    The printed gamma equation and SI normalization are independently validated
    by reproduce_speziale_2001_mgo. Exact author/later pressure reducers remain
    unavailable; this is a source-equation calibration sensitivity.
    Constants: Speziale Eqs. 1-4, 10-11 and Fei (2004) Table 3.
    """
    v, t = np.broadcast_arrays(np.asarray(volume), np.asarray(temperature))
    logv = np.log(v / 74.71)
    gamma = 1.524 * np.exp(1.65 / 11.8 * np.expm1(11.8 * logv))
    u = logv[..., None] * (NODES + 1) / 2
    integrand = 1.524 * np.exp(1.65 / 11.8 * np.expm1(11.8 * u))
    theta = 773 * np.exp(-logv / 2 * np.sum(WEIGHTS * integrand, axis=-1))
    x = np.exp(-logv / 3)
    cold = 1.5 * 160.2 * (x**7 - x**5) * (1 + 0.75 * (3.99 - 4) * (x * x - 1))
    delta = energy(theta, t, 2) - energy(theta, 300, 2)
    return cold + gamma * delta / (v * N_A * 1e-30 / 4) / 1e9


def derivative(function, x, step):
    return (function(x + step) - function(x - step)) / (2 * step)


def read_csv(name):
    with (DATA / name).open() as stream:
        return list(csv.DictReader(stream))


def observations():
    rows = []
    for r in read_csv("gold-fei-2004-table1.csv"):
        rows.append(
            dict(
                source="Fei 2004 Table 1",
                row=r["run"],
                **{
                    k: float(r[key])
                    for k, key in (
                        ("v", "volume_a3"),
                        ("t", "temperature_k"),
                        ("p", "pressure_gpa"),
                        ("sv", "volume_uncertainty_a3"),
                        ("sp", "pressure_gpa_uncertainty"),
                        ("vm", "mgo_volume_a3"),
                        ("svm", "mgo_volume_uncertainty_a3"),
                    )
                },
            )
        )
    for r in read_csv("gold-hirose-2006-table1.csv"):
        if r["use_in_fei2007_thermal_replay"] != "true":
            continue
        rows.append(
            dict(
                source="Hirose 2006 Table 1",
                row="run2-cycle" + r["cycle"],
                v=float(r["volume_a3"]),
                t=float(r["temperature_k"]),
                p=float(r["mgo_speziale_pressure_gpa"]),
                sv=np.nan,
                sp=float(r["mgo_speziale_pressure_gpa_uncertainty"]),
                vm=float(r["mgo_volume_a3"]),
                svm=np.nan,
            )
        )
    arrays = {
        k: np.array([r[k] for r in rows])
        for k in ("v", "t", "p", "sv", "sp", "vm", "svm")
    }
    return rows, arrays


def metrics(residual):
    return dict(
        rows=int(len(residual)),
        rmse_gpa=float(np.sqrt(np.mean(residual**2))),
        mean_residual_gpa=float(np.mean(residual)),
        max_abs_residual_gpa=float(np.max(np.abs(residual))),
    )


def effective_sigma(d, coefficients=PUBLISHED):
    """Available terms only; missing Au errors remain missing in the source."""
    pv = derivative(lambda v: pressure(v, d["t"], coefficients), d["v"], 1e-4)
    contribution = np.where(np.isfinite(d["sv"]), (pv * d["sv"]) ** 2, 0.0)
    return np.sqrt(d["sp"] ** 2 + contribution)


def fit(d, free=(4,), mode="equal", law="printed", base=PUBLISHED, prefactor=1):
    base = np.array(base, copy=True)
    sigma = (
        np.ones_like(d["p"])
        if mode == "equal"
        else d["sp"]
        if mode == "pressure_only"
        else effective_sigma(d, base)
    )

    def unpack(values):
        c = base.copy()
        c[list(free)] = values
        return c

    result = least_squares(
        lambda z: (
            (pressure(d["v"], d["t"], unpack(z), law, prefactor) - d["p"]) / sigma
        ),
        base[list(free)],
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
        max_nfev=2000,
    )
    residual = pressure(d["v"], d["t"], unpack(result.x), law, prefactor) - d["p"]
    dof = len(d["p"]) - len(free)
    scale = float(result.fun @ result.fun / dof) if mode == "equal" else 1.0
    covariance = np.linalg.inv(result.jac.T @ result.jac) * scale
    return {
        "free_parameters": [NAMES[i] for i in free],
        "fixed_parameters": {
            NAMES[i]: float(base[i]) for i in range(6) if i not in free
        },
        "parameters": dict(zip([NAMES[i] for i in free], result.x.tolist())),
        "conditional_curvature_widths": dict(
            zip([NAMES[i] for i in free], np.sqrt(np.diag(covariance)).tolist())
        ),
        "width_convention": "Residual-scaled linearized standard errors"
        if mode == "equal"
        else "Absolute inverse-curvature widths assuming the supplied parenthetical measurement widths are 1-sigma; their original coverage and missing errors are not inferred",
        "covariance": covariance.tolist(),
        "covariance_parameter_order": [NAMES[i] for i in free],
        "objective": mode,
        "law": law,
        "energy_prefactor_relative_to_3R": float(prefactor),
        "objective_sum": float(result.fun @ result.fun),
        "nominal_dof": dof,
        "published": metrics(pressure(d["v"], d["t"], base, law, prefactor) - d["p"]),
        "refit": metrics(residual),
        "residuals_gpa": residual.tolist(),
        "pressure_residual_scales_gpa": sigma.tolist(),
        "solver_success": bool(result.success),
        "author_fit_reproduced": False,
    }


def available_error_eiv(d):
    """Latent Au volumes where errors exist; Hirose V is fixed conditionally.

    Minimize pressure residual/sP and measured-volume displacement/sV together.
    Missing predictor errors are not fitted or imputed.
    """
    mask = np.isfinite(d["sv"])

    def state(z):
        c = PUBLISHED.copy()
        c[4] = z[0]
        volume = d["v"].copy()
        volume[mask] += z[1:] * d["sv"][mask]
        return c, volume

    def residual(z):
        c, volume = state(z)
        return np.r_[(pressure(volume, d["t"], c) - d["p"]) / d["sp"], z[1:]]

    result = least_squares(
        residual,
        np.r_[0.6, np.zeros(mask.sum())],
        xtol=1e-11,
        ftol=1e-11,
        gtol=1e-11,
    )
    c, volume = state(result.x)
    return dict(
        q=float(result.x[0]),
        conditional_curvature_width_assuming_1sigma=float(
            np.sqrt(np.linalg.inv(result.jac.T @ result.jac)[0, 0])
        ),
        objective_sum=float(result.fun @ result.fun),
        measured_coordinate_residuals=metrics(pressure(d["v"], d["t"], c) - d["p"]),
        latent_volume_shifts_a3=(volume - d["v"]).tolist(),
        volume_displacements_in_reported_widths=result.x[1:].tolist(),
        solver_success=bool(result.success),
        qualification="Conditional errors-in-variables sensitivity using 26 known Au volume widths and 28 reported pressure widths. Missing Hirose volume and all temperature/calibration covariance remain unavailable. Assuming quoted widths are 1-sigma is only the curvature-width convention, not a recovered source confidence statement. No author weights recovered.",
        author_fit_reproduced=False,
    )


def cold_diagnostics():
    rows = read_csv("gold-dewaele-2004-table1-compression.csv")
    v = np.array([4 * float(r["atomic_volume_a3"]) for r in rows])
    p = np.array([float(r["ruby_pressure_revised_gpa"]) for r in rows])
    d = dict(
        v=v,
        t=np.full(len(v), 300.0),
        p=p,
        sp=np.full(len(v), np.nan),
        sv=np.full(len(v), 0.04),
    )
    result = {"dewaele_only_equal_pressure": fit(d, free=(2,))}

    # Inverse P-V fit: pressure errors unavailable, only published volume width.
    def inverse_residual(z):
        c = PUBLISHED.copy()
        c[2] = z[0]
        solved = least_squares(
            lambda vv: pressure(vv, 300, c) - p, v, xtol=1e-12, ftol=1e-12, gtol=1e-12
        )
        return (solved.x - v) / 0.04

    sol = least_squares(inverse_residual, [6.0], xtol=1e-10, ftol=1e-10, gtol=1e-10)
    result["dewaele_only_volume_error"] = dict(
        K0_prime=float(sol.x[0]),
        conditional_curvature_width=float(1 / np.sqrt((sol.jac.T @ sol.jac)[0, 0])),
        solver_success=bool(sol.success),
        qualification="Minimize squared volume residual/.04 cell A^3 at printed revised-ruby pressures; conditional on exact pressures. Dewaele's global pressure-error envelope is not a row-wise sigma.",
    )
    more = read_csv("gold-fei-2007-figure1-digitized.csv")
    h = next(
        r
        for r in read_csv("gold-hirose-2006-table1.csv")
        if r["temperature_k"] == "300"
    )
    d.update(
        v=np.r_[
            v,
            [float(r["volume_a3_conventional_cell"]) for r in more],
            float(h["volume_a3"]),
        ],
        p=np.r_[
            p,
            [float(r["pressure_gpa"]) for r in more],
            float(h["mgo_speziale_pressure_gpa"]),
        ],
    )
    d["t"] = np.full(len(d["v"]), 300.0)
    d["sv"] = np.r_[d["sv"], np.full(7, np.nan)]
    d["sp"] = np.r_[d["sp"], np.full(7, np.nan)]
    result["all_available_cold_equal_pressure"] = fit(d, free=(2,))
    result["qualification"] = (
        "V0/K0 fixed. Original Dewaele 298 K table evaluated at Fei's 300 K anchor; exact-298 K sensitivity separately reported. Six new Fei rows are digitized measured markers, not calculated curves; no experimental errors are invented. One Hirose cold row lacks errors. Cold-data selection and original author weights remain unknown."
    )
    d298 = dict(v=v, t=np.full(len(v), 298.0), p=p)
    result["dewaele_at_exact_298k_equal_pressure_sensitivity"] = fit(d298, free=(2,))
    return result, d


def curve_checks():
    source = json.loads(CURVES.read_text(encoding="utf-8"))
    out = []
    for c in source["curves"]:
        xy = np.array(c["points_xy_pt"])
        a = source["axis_calibration"]
        v = np.polyval(a["volume_from_y"], xy[:, 1])
        p = np.polyval(a["pressure_from_x"], xy[:, 0])
        width = c["linewidth_pt"] * a["pressure_from_x"][0]
        laws = {}
        for law in ("printed", "integrated"):
            delta = pressure(v, c["temperature_k"], law=law) - p
            laws[law] = dict(
                **metrics(delta),
                within_one_pressure_direction_stroke=bool(np.all(abs(delta) < width)),
            )
        out.append(
            dict(
                temperature_k=c["temperature_k"],
                pressure_direction_stroke_gpa=width,
                laws=laws,
                max_law_difference_gpa=float(
                    np.max(
                        abs(
                            pressure(v, c["temperature_k"])
                            - pressure(v, c["temperature_k"], law="integrated")
                        )
                    )
                ),
            )
        )
    return dict(
        kind="calculated_source_curves_never_fit_targets",
        inherited_source_pdf_sha256=source["source_pdf_sha256"],
        curves=out,
        qualification="Archived vector coordinates from primary Figure 1; source PDF hash inherited from the earlier extraction. Original PDF could not be fetched anew (HTTP 403); original cached page image and full text inspected. Graph agreement alone does not establish fit parity or uniquely identify theta(V).",
    )


def reproduce():
    rows, d = observations()
    small = {k: v[:26] for k, v in d.items()}
    cold, cd = cold_diagnostics()
    fits = {}
    for label, data in (("fei2004_only", small), ("all28", d)):
        for mode in ("equal", "available_errors", "pressure_only"):
            fits[label + "_" + mode] = fit(data, mode=mode)
    fits["all28_integrated_equal"] = fit(d, law="integrated")
    fits["all28_integrated_available_errors"] = fit(
        d, law="integrated", mode="available_errors"
    )
    mass_ratio = 196.96657 * 0.125 / (3 * R)
    fits["all28_fei2004_mass_prefactor_sensitivity"] = fit(d, prefactor=mass_ratio)
    fits["all28_gamma0_q_free_sensitivity"] = fit(d, free=(3, 4))
    # Anchor q fits to independently refitted cold K' values, retaining gamma/theta.
    for label in ("dewaele_only_equal_pressure", "all_available_cold_equal_pressure"):
        base = PUBLISHED.copy()
        base[2] = cold[label]["parameters"]["K0_prime"]
        fits[label + "_then_q"] = fit(d, base=base)
    joint = dict(
        v=np.r_[cd["v"], d["v"]], t=np.r_[cd["t"], d["t"]], p=np.r_[cd["p"], d["p"]]
    )
    fits["cold44_hot28_joint_Kprime_gamma0_q_sensitivity"] = fit(joint, free=(2, 3, 4))
    sensitivity = {}
    for index, width in ((0, 0.004), (2, 0.02), (3, 0.03)):
        variants = []
        for sign in (-1, 1):
            base = PUBLISHED.copy()
            base[index] += sign * width
            variants.append(fit(d, base=base)["parameters"]["q"])
        sensitivity[NAMES[index]] = dict(
            published_width=width,
            q_at_minus_plus_width=variants,
            qualification="Individual fixed-anchor perturbation, not a joint confidence interval",
        )
    vg = 67.85 * np.linspace(0.72, 1.04, 17)[:, None]
    tg = np.array([300.0, 500, 1000, 1473, 2000, 2173, 2330])[None, :]
    parity = get_eos_record("gold_fei_2007_vinet_2").pressure(vg, tg) - pressure(vg, tg)
    quadrature = fei_pressure(vg, tg) - pressure(vg, tg)
    reconstructed = mgo_pressure(d["vm"], d["t"])
    recal = dict(d, p=reconstructed)
    fits["integrated_MgO_reduction_sensitivity"] = fit(recal)
    pv = derivative(lambda v: pressure(v, d["t"]), d["v"], 1e-4)
    au_t = derivative(lambda t: pressure(d["v"], t), d["t"], 0.01)
    mg_t = derivative(lambda t: mgo_pressure(d["vm"], t), d["t"], 0.01)
    # Zero placeholders are not exported as recovered errors.
    points = []
    for i, r in enumerate(rows):
        points.append(
            dict(
                source=r["source"],
                source_row=r["row"],
                volume_a3=r["v"],
                temperature_k=r["t"],
                reported_MgO_pressure_gpa=r["p"],
                reported_pressure_width_gpa=r["sp"],
                Au_volume_width_a3=float(r["sv"]) if np.isfinite(r["sv"]) else None,
                Au_volume_contribution_gpa=abs(float(pv[i] * r["sv"]))
                if np.isfinite(r["sv"])
                else None,
                partial_effective_pressure_width_gpa=float(effective_sigma(d)[i]),
                published_minus_reported_gpa=float(pressure(r["v"], r["t"]) - r["p"]),
                integrated_MgO_pressure_gpa=float(reconstructed[i]),
                integrated_MgO_minus_reported_gpa=float(reconstructed[i] - r["p"]),
                Au_dP_dT_gpa_k=float(au_t[i]),
                MgO_dP_dT_gpa_k=float(mg_t[i]),
                shared_temperature_residual_slope_gpa_k=float(au_t[i] - mg_t[i]),
            )
        )
    paths = [
        DATA / name
        for name in (
            "gold-fei-2004-table1.csv",
            "gold-dewaele-2004-table1-compression.csv",
            "gold-fei-2007-figure1-digitized.csv",
            "gold-hirose-2006-table1.csv",
            "gold-hirose-2006-source.json",
            "gold-fei-2007-replay-source.json",
            "fei-2004-pressure-scales-source.json",
        )
    ] + [CURVES]
    return dict(
        format="peritheos.fei-2007-gold-conditional-replay",
        audit_date="2026-10-09",
        classification=classification(cold, fits),
        primary_source="https://pmc.ncbi.nlm.nih.gov/articles/PMC1890468/",
        doi="10.1073/pnas.0609013104",
        source_locations=[
            "Table 1, Equations 2-3 and Figure 1 on page 9183",
            "Au refit paragraph continuing onto page 9184",
            "Fei 2004 Table 1",
            "Hirose 2006 Table 1",
            "Dewaele 2004 Table I",
        ],
        input_sha256={
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in paths
        },
        published_parameters=dict(zip(NAMES, PUBLISHED.tolist())),
        published_error_widths=dict(
            V0=0.004, K0=None, K0_prime=0.02, gamma0=0.03, q=0.3, theta0=None
        ),
        published_error_confidence=None,
        published_covariance=None,
        independent_equation_max_difference_gpa=float(np.max(abs(parity))),
        independent_quadrature_max_difference_gpa=float(np.max(abs(quadrature))),
        independent_equation_states=int(parity.size),
        source_curve_checks=curve_checks(),
        selection=dict(
            thermal_rows=28,
            Fei2004_rows=26,
            Hirose2006_rows=2,
            cold_rows=44,
            author_exact_machine_readable_selection_recovered=False,
            evidence="Prose explicitly refits Fei 2004 and Hirose 2006 PVT data; Figure 1 identifies the two Run 2 hot Hirose states. Fei 2004 Table 1 is the paired Au-MgO table. Fei 2007 Table 1 Au footnote omits Hirose despite the prose. Full 26-row selection is source-supported but not author-confirmed; Figure 1 plots only selected Fei isotherms.",
        ),
        fits=fits,
        available_error_eiv=available_error_eiv(d),
        cold_diagnostics=cold,
        anchor_sensitivity=sensitivity,
        points=points,
        MgO_reduction=dict(
            model="Speziale variable-q gamma with thermodynamically integrated Debye temperature; fixed constants",
            source="https://duffy.psb-test.princeton.edu/sites/g/files/toruqf616/files/speziale_et_al-2001-jgrse.pdf",
            Fei2004_difference=metrics(reconstructed[:26] - d["p"][:26]),
            Hirose2006_difference=metrics(reconstructed[26:] - d["p"][26:]),
            status="not_exactly_reproduced",
            followup_report="docs/data/speziale-2001-mgo-reproduction.json",
            qualification="The printed MgO gamma equation, differential identities, quadrature and SI normalization have been independently validated. Substitution of q(V) into gamma0*(V/V0)**q, a different gamma law, closely reproduces later reported pressures; that is diagnostic inference, not recovered author code or Speziale Eq. 11. Exact original thermal fitting and later reducer details remain unavailable. Reported targets and Au parity are unchanged.",
        ),
        objectives=dict(
            equal="Minimize sum(P_Au-P_reported_MgO)^2 at measured V,T.",
            available_errors="Minimize sum((P_Au-P_reported_MgO)/s)^2; s^2=sP^2+(dP_Au/dV*sV)^2 using available widths. Scales frozen at adopted anchor and published q=.6; predictor errors treated independent. Hirose sV absent: only its reported sP term is available, not a claim of zero volume error.",
            pressure_only="Uses only source pressure widths; omits known Au predictor errors. Weight-sensitivity control, not the preferred measurement-error diagnostic.",
            shared_temperature="For re-reduced paired measurements, Var(residual) includes (dP_Au/dT-dP_MgO/dT)^2*sT^2, not the sum of separate squared slopes. A full propagation also needs Au-MgO lattice covariance, calibrant parameter covariance and correlations across rows. These are unavailable. Hirose reported pressure widths already include temperature contributions: adding a second temperature term would double count without decomposing them.",
        ),
        limitations=[
            "No author software, exact weights, unrounded inputs, complete source selection file or confidence/covariance supplied.",
            "Source parenthetical widths retained without assigning a confidence level; conditional curvature widths are model assumptions, not recovered author errors.",
            "No row-wise Fei temperature errors; its Type-C thermocouple temperatures were not corrected for pressure effects. Reported spatial gradient ~30 C/500 um over <=250 um is not a standard deviation.",
            "Hirose Run 2 Au/MgO lattice errors absent; spatial variation <+/-10% is not a row-wise temperature sigma. Weighted objectives are incomplete.",
            "Dewaele reports a global pressure-error range but no row-wise pressure sigmas; digitized Fei cold marker errors are graphical, not experimental.",
            "The Fei 2004 mass-specific Au prefactor 0.125 J/g/K differs from modern 3R/M by ~1.29%; included as normalization sensitivity only, not adopted for Fei 2007.",
            "Printed theta law has -dln(theta)/dln(V)=gamma*(1+q*ln(V/V0)); pressure replay does not validate caloric consistency. Integrated law comparison is a sensitivity, not a default change.",
            "Available source inputs are training comparisons, not withheld validation or evidence of absolute model accuracy.",
            "Au pressure calibration is explicitly recorded: RT uses Dewaele revised ruby and hot targets use Speziale MgO. The independent MgO implementation is a sensitivity, not an exact recovered author reduction.",
        ],
        author_fit_reproduced=False,
        published_coefficients_changed=False,
        defaults_changed=False,
    )


def classification(cold, fits):
    """Conditional numerical parity against printed widths, without assigning sigma."""
    cold_fit = cold["dewaele_only_equal_pressure"]
    hot_fit = fits["dewaele_only_equal_pressure_then_q"]
    parameters = []
    for name, published, width, fitted in (
        ("K0_prime", 6.0, 0.02, cold_fit),
        ("q", 0.6, 0.3, hot_fit),
    ):
        value = fitted["parameters"][name]
        difference = value - published
        parameters.append(
            dict(
                parameter=name,
                published=published,
                published_error=width,
                refit=value,
                refit_error=fitted["conditional_curvature_widths"][name],
                difference=difference,
                relative_difference=abs(difference) / published,
                within_combined_2sigma=None,
                within_reported_error=abs(difference) <= width,
                similar=bool(
                    abs(difference) <= (1.0 if name == "K0_prime" else 0.25 * published)
                ),
            )
        )
    return dict(
        status="parity"
        if all(p["within_reported_error"] and p["similar"] for p in parameters)
        else "parity_not_achieved",
        parity_basis="within_reported_parameter_errors",
        parameters=parameters,
        qualification="Conditional numerical parity: staged equal-pressure-weight replay of 37 revised-ruby cold observations followed by 28 reported-MgO hot observations reproduces K0-prime and q within printed parameter error widths. V0/K0/gamma0/theta0 are adopted and fixed. Source weights, exact MgO reduction and uncertainty confidence/covariance remain unresolved; no formal combined-two-sigma or source-exact author-fit claim.",
        author_fit_reproduced=False,
    )


def ledger_outcome(record):
    """Dedicated staged Au result consumed by the shared reproduction ledger."""
    cold, _ = cold_diagnostics()
    _, data = observations()
    base = PUBLISHED.copy()
    base[2] = cold["dewaele_only_equal_pressure"]["parameters"]["K0_prime"]
    hot = fit(data, base=base)
    result = classification(cold, {"dewaele_only_equal_pressure_then_q": hot})
    return dict(
        **result,
        dataset_identifiers=[
            "gold_dewaele_2004_table1_compression",
            "gold_fei_2004_table1",
            "gold_hirose_2006_table1",
        ],
        observations=65,
        cold_observations=37,
        thermal_observations=28,
        fit_kind="staged_vinet_mgd_conditional_replay",
        objective="equal_pressure_residuals",
        free_parameters=["K0_prime", "q"],
        fixed_parameters=["V0", "K0", "gamma0", "theta0"],
        primary_data_status="bundled",
        sigma_absolute=False,
        rmse_gpa=hot["refit"]["rmse_gpa"],
        published_rmse_gpa=metrics(pressure(data["v"], data["t"]) - data["p"])[
            "rmse_gpa"
        ],
        residual_scope="Thermal 28-row RMSE; cold-stage residuals archived separately",
        reproduction_report="docs/data/fei-2007-gold-reproduction.json",
    )


def check_saved(saved, current, *, rtol=1e-5, atol=2e-6):
    if isinstance(current, dict):
        assert saved.keys() == current.keys()
        for k, v in current.items():
            # Latent EIV shifts are expressed in measured line widths. Different
            # numerical Jacobians can move the optimum by 0.0001 of a width.
            if k == "available_error_eiv":
                check_saved(saved[k], v, rtol=1e-4, atol=1e-4)
            else:
                check_saved(saved[k], v, rtol=rtol, atol=atol)
    elif isinstance(current, list):
        assert len(saved) == len(current)
        for a, b in zip(saved, current):
            check_saved(a, b, rtol=rtol, atol=atol)
    elif isinstance(current, float):
        # Optimizer-derived errors vary slightly across SciPy/BLAS builds.
        assert np.isclose(saved, current, rtol=rtol, atol=atol), (saved, current)
    else:
        assert saved == current, (saved, current)


def save_residuals(report):
    fields = list(report["points"][0])
    with RESIDUALS.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, lineterminator="\n")
        writer.writeheader()
        writer.writerows(report["points"])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    if args.check:
        check_saved(json.loads(OUTPUT.read_text(encoding="utf-8")), report)
    else:
        OUTPUT.write_text(
            json.dumps(report, indent=2, sort_keys=True, allow_nan=False) + "\n",
            encoding="utf-8",
        )
        save_residuals(report)
    for label in ("all28_equal", "all28_available_errors", "fei2004_only_equal"):
        r = report["fits"][label]
        print(label, r["parameters"], r["refit"])


if __name__ == "__main__":
    main()
