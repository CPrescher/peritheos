"""Audit Katsura (2009) wadsleyite without changing catalog records.

The independent SI calculation distinguishes source observations, digitized
source curves, conditional refits, and deliberately unphysical normalization
diagnostics. Run as a repository maintenance module, not an installed API.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

import numpy as np
from scipy.integrate import quad
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets/wadsleyite-katsura-2009-table1-pvt.csv"
SOURCE = ROOT / "docs/data/wadsleyite-katsura-2009-source-checks.json"
OUTPUT = ROOT / "docs/data/wadsleyite-katsura-2009-audit.json"
V0 = 538.49
K0 = 169.2
R = 8.31446261815324
NA = 6.02214076e23
Z = 8.0
N = 7.0
NODES, WEIGHTS = np.polynomial.legendre.leggauss(96)
NODES, WEIGHTS = (NODES + 1.0) / 2.0, WEIGHTS / 2.0


def unpack(token: str) -> tuple[float, float]:
    value, digits = token[:-1].split("(")
    return float(value), int(digits) * 10.0 ** (-len(value.split(".")[1]))


def load_data() -> dict[str, np.ndarray]:
    with DATA.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    source = json.loads(SOURCE.read_text())
    tokens = source["table1_tokens"]
    if len(rows) != len(tokens) or len(rows) != 85:
        raise ValueError("Source table row count changed")
    for row, token in zip(rows, tokens):
        p, sp = unpack(token["pressure_gpa_token"])
        x, sx = unpack(token["volume_ratio_token"])
        expected = [float(token["temperature_k_token"]), p, sp, x, sx]
        if not np.allclose(
            list(map(float, row.values())), expected, atol=1e-12, rtol=0
        ):
            raise ValueError(f"Source token mismatch at row {token['row']}")
    return {key: np.array([float(row[key]) for row in rows]) for key in rows[0]}


def verify_pdf(path: Path) -> None:
    """Independently extract the full two-column table from the final PDF."""
    source = json.loads(SOURCE.read_text())
    if hashlib.sha256(path.read_bytes()).hexdigest() != source["article_pdf"]["sha256"]:
        raise ValueError("PDF differs from the visually inspected final article")
    text = subprocess.check_output(
        ["pdftotext", "-raw", "-f", "3", "-l", "3", str(path), "-"], text=True
    )
    extracted = re.findall(
        r"(?m)^\s*(\d{3,4})\s+(\d+\.\d+\(\d+\))\s+(0\.\d+\(\d+\))\s*$", text
    )
    expected = [
        (r["temperature_k_token"], r["pressure_gpa_token"], r["volume_ratio_token"])
        for r in source["table1_tokens"]
    ]
    if extracted != expected:
        raise ValueError("Final PDF and pinned source tokens differ")


def cold_pressure(x, kp=4.1):
    eta = np.asarray(x, dtype=float) ** (-1.0 / 3.0)
    return 1.5 * K0 * (eta**7 - eta**5) * (1 + 0.75 * (kp - 4) * (eta**2 - 1))


def energy(theta, temperature, n=N):
    theta, temperature = np.broadcast_arrays(theta, temperature)
    y = theta / temperature
    argument = y[..., None] * NODES
    integral = y * np.sum(WEIGHTS * argument**3 / np.expm1(argument), axis=-1)
    return 9 * n * R * temperature * integral / y**3


def pressure(
    x,
    temperature,
    *,
    kp=4.1,
    gamma0=1.64,
    q=1.5,
    n=N,
    z=Z,
    law="integrated_gruneisen",
    tr=300.0,
    gamma_ratio=None,
):
    """Return GPa using m3/mol formula units, with no Peritheos evaluation."""
    x, temperature = np.broadcast_arrays(x, temperature)
    base = x if gamma_ratio is None else np.asarray(gamma_ratio)
    gamma = gamma0 * base**q
    if law == "integrated_gruneisen":
        exponent = -gamma0 * np.log(base) if q == 0 else (gamma0 - gamma) / q
        theta = 814.0 * np.exp(exponent)
    elif law == "variable_exponent":
        theta = 814.0 * base ** (-gamma)
    else:
        raise ValueError(law)
    molar_volume_m3 = x * V0 * 1e-30 * NA / z
    thermal = gamma * (energy(theta, temperature, n) - energy(theta, tr, n))
    return cold_pressure(x, kp) + thermal / molar_volume_m3 / 1e9


def volume_ratio(p, t, **options):
    return brentq(lambda x: float(pressure(x, t, **options)) - p, 0.80, 1.18)


def metrics(residual):
    a = np.asarray(residual)
    return dict(
        rows=len(a),
        rms_gpa=float(np.sqrt(np.mean(a * a))),
        mean_gpa=float(np.mean(a)),
        minimum_gpa=float(np.min(a)),
        maximum_gpa=float(np.max(a)),
        maximum_absolute_gpa=float(np.max(abs(a))),
    )


def later_author_code_pressure(x, temperature):
    """Translate the deposited 2022 pressure routine with the 2009 inputs.

    This is a cross-check against a later author implementation, not execution
    of MATLAB or recovery of the original 2009 fitting/plotting program.
    """
    path = (
        ROOT
        / "docs/data/forsterite-katsura-2022-sources/forsterite-katsura-2022-author-geotherm.m"
    )
    payload = path.read_bytes()
    expected = "26375509cb6605523f45a8f127a28cde6be0d22de5353e0cee3134a559373251"
    if hashlib.sha256(payload).hexdigest() != expected:
        raise ValueError("Deposited author code differs from the reviewed version")
    text = payload.decode("utf-8", errors="replace")
    energy_source = text.split("function E_th = ThermalEnergy", 1)[1]
    vector = re.search(r"df\s*=\s*\[([^]]+)\]", energy_source).group(1)
    coefficients = np.array(
        [float(value) for value in vector.replace("...", " ").split()]
    )
    gas_constant = float(re.search(r"R\s*=\s*([0-9.]+);", energy_source).group(1))
    avogadro = float(re.search(r"N_Avogadro\s*=\s*([0-9.E+]+);", text).group(1))
    x, temperature = np.broadcast_arrays(x, temperature)
    gamma = 1.64 * x**1.5
    theta = 814 * np.exp(-(gamma - 1.64) / 1.5)
    vm = avogadro * V0 * x * 1e-30 / Z

    def polynomial_energy(t):
        return (
            9
            * N
            * gas_constant
            * t
            * np.polynomial.polynomial.polyval(theta / t, coefficients)
        )

    # Literal BM3 factors from the author's separate BM3 function.
    cold = (
        (3 / 2)
        * K0
        * (x ** (-7 / 3) - x ** (-5 / 3))
        * (1 + (3 / 4) * (4.1 - 4) * (x ** (-2 / 3) - 1))
    )
    total = (
        cold
        + gamma / vm * (polynomial_energy(temperature) - polynomial_energy(300)) * 1e-9
    )
    return total, dict(
        source_doi="10.5281/zenodo.5903286",
        original_filename="Katsura Geotherm calculation v2.m",
        version="1.1.0, published 25 January 2022",
        sha256=expected,
        source_routines=["MGD_EOS_RJ", "ThermalEnergy", "BM3"],
        evaluation="Literal Python translation of pressure/energy formulas; MATLAB program not executed",
        coefficient_target="2009 wadsleyite coefficients, not the different 2022 fit vector",
        polynomial_coefficients=coefficients.tolist(),
        gas_constant_j_mol_k=gas_constant,
        avogadro_per_mol=avogadro,
        source_2022_wadsleyite_n=7,
        source_2022_wadsleyite_z=8,
        original_2009_program_recovered=False,
    )


def conditional_fit(x, t, p, kp):
    """High-pressure-only diagnostic; never the omitted Inoue joint objective."""
    fit = least_squares(
        lambda a: pressure(x, t, kp=kp, gamma0=a[0], q=a[1]) - p,
        [1.64, 1.5],
        bounds=([0.2, -2.0], [3.0, 8.0]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    dof = len(p) - 2
    covariance = np.linalg.inv(fit.jac.T @ fit.jac) * (fit.fun @ fit.fun) / dof
    return dict(
        gamma0=float(fit.x[0]),
        q=float(fit.x[1]),
        fixed_k0_prime=kp,
        conditional_standard_errors=np.sqrt(np.diag(covariance)).tolist(),
        conditional_covariance=covariance.tolist(),
        covariance_order=["gamma0", "q"],
        degrees_of_freedom=dof,
        objective="equal pressure residuals; high-pressure data only",
        free_parameters=["gamma0", "q"],
        fixed_parameters={
            "V0_cell_a3": V0,
            "K0_gpa": K0,
            "K0_prime": kp,
            "theta0_k": 814.0,
            "Tr_k": 300.0,
            "n": N,
            "Z": Z,
        },
        initial_parameters=[1.64, 1.5],
        parameter_bounds=[[0.2, -2.0], [3.0, 8.0]],
        solver="scipy.optimize.least_squares, default trf, finite-difference Jacobian",
        solver_tolerances={"xtol": 1e-12, "ftol": 1e-12, "gtol": 1e-12},
        original_joint_fit_reproduced=False,
        uncertainty="RSS/dof-scaled local covariance; no calibration, fixed-input, predictor or shared-data covariance",
        success=bool(fit.success),
        active_bounds=fit.active_mask.tolist(),
        residuals=metrics(fit.fun),
    )


def audit() -> dict:
    data = load_data()
    x, t, p = (data[k] for k in ["normalized_volume", "temperature_k", "pressure_gpa"])
    rt = t <= 320
    source = json.loads(SOURCE.read_text())
    baseline = pressure(x, t)
    author_pressure, author_check = later_author_code_pressure(x, t)
    author_check["source_2009_coefficients_residuals"] = metrics(author_pressure - p)
    author_check["maximum_pressure_difference_from_exact_quadrature_gpa"] = float(
        np.max(abs(author_pressure - baseline))
    )
    variants = {}
    for name, options in [
        ("physical_n7_integrated", {}),
        ("physical_n7_variable_exponent", {"law": "variable_exponent"}),
        ("reference_298k_sensitivity", {"tr": 298.0}),
        ("unphysical_n5_diagnostic_only", {"n": 5}),
    ]:
        calculated = pressure(x, t, **options)
        variants[name] = dict(
            all_rows=metrics(calculated - p),
            near_rt_rows=metrics((calculated - p)[rt]),
            heated_rows=metrics((calculated - p)[~rt]),
        )

    # Eq. (4) of the cited 2009 forsterite method uses V(P,300). Test that
    # notation as an implicit pressure equation, without adopting a new EOS.
    def cold_reference_prediction(xx, tt):
        def residual(pp):
            xc = brentq(lambda vv: cold_pressure(vv) - pp, 0.7, 1.2)
            return pp - float(pressure(xx, tt, gamma_ratio=xc))

        return brentq(residual, 0, 45)

    cold_reference = np.array(
        [cold_reference_prediction(xx, tt) for xx, tt in zip(x, t)]
    )
    variants["source_notation_cold_reference_gamma_n7"] = dict(
        all_rows=metrics(cold_reference - p),
        scope="Implicit equation with gamma/theta evaluated at BM3 V(P,300); denominator remains actual molar V. Not adopted.",
    )
    cold_fit = least_squares(lambda a: cold_pressure(x[rt], a[0]) - p[rt], [4.1])
    actual_t_fit = least_squares(
        lambda a: pressure(x[rt], t[rt], kp=a[0]) - p[rt], [4.1]
    )
    thermal = baseline - cold_pressure(x)
    factor = float(thermal @ (p - cold_pressure(x)) / (thermal @ thermal))
    from peritheos.eos.rt import BM3
    from peritheos.eos.thermal import MieGruneisenDebye

    scale = NA * 1e-25 / Z
    native = MieGruneisenDebye(
        BM3(V0 * scale, K0, 4.1), Tr=300, theta0=814, gamma0=1.64, q=1.5, n=N
    )
    native_difference = float(
        np.max(abs(native.pressure(x * V0 * scale, t) - baseline))
    )
    native_available = getattr(native, "_native", None) is not None
    native._native = None
    python_difference = float(
        np.max(
            abs(
                native.thermal_pressure(x * V0 * scale, t) + cold_pressure(x) - baseline
            )
        )
    )
    figure = source["figure2"]
    pa = figure["pressure_axis"]
    va = figure["dilatation_axis"]
    pf = (
        pa["pressure_gpa"]
        + (figure["sampling_pixel_x"] - pa["pixel_x"]) / pa["pixels_per_gpa"]
    )
    source_ratios = np.array(
        [
            va["volume_ratio"]
            + (va["pixel_y"] - r["pixel_y"]) / va["pixels_per_unit_ratio"]
            for r in figure["samples"]
        ]
    )
    temperatures = np.array([r["temperature_k"] for r in figure["samples"]])
    figure_results = {}
    ambient = {}
    for n in [7, 5]:
        vref = volume_ratio(pf, 300, n=n)
        calculated = np.array([volume_ratio(pf, tt, n=n) / vref for tt in temperatures])
        delta = calculated - source_ratios
        figure_results[f"n{n}"] = dict(
            calculated_ratios=calculated.tolist(),
            maximum_absolute_ratio_difference=float(np.max(abs(delta))),
            rms_ratio_difference=float(np.sqrt(np.mean(delta**2))),
            within_graphical_bound=bool(
                np.all(abs(delta) <= figure["reading_bound_volume_ratio"])
            ),
        )
        ambient[f"n{n}"] = dict(
            mean_alpha_293_973_k=float(
                np.log(volume_ratio(0, 973, n=n) / volume_ratio(0, 293, n=n)) / 680
            ),
            source_reported_mean_alpha=2.7e-5,
            note="Exploratory cross-model comparison only: the reported expansivity is associated with the separate HTBM fit, so exact agreement with MGD is not required and this is not a decisive MGD benchmark.",
        )
    per_t = []
    for tt in np.unique(t):
        mask = t == tt
        per_t.append(dict(temperature_k=float(tt), **metrics((baseline - p)[mask])))
    alpha_quad = 814 / 1500
    adaptive_integral = quad(lambda z: z**3 / np.expm1(z), 0, alpha_quad, epsabs=1e-12)[
        0
    ]
    quadrature_difference = abs(
        float(energy(814.0, 1500.0))
        - 9 * N * R * 1500 * adaptive_integral / alpha_quad**3
    )
    return dict(
        format="peritheos.wadsleyite-2009-audit",
        format_version=1,
        audit_date="2026-10-30",
        sources=source,
        dataset_sha256=hashlib.sha256(DATA.read_bytes()).hexdigest(),
        material_file_sha256=hashlib.sha256(
            (ROOT / "peritheos/data/materials/wadsleyite.eosmat").read_bytes()
        ).hexdigest(),
        convention=dict(
            volume="angstrom^3 per conventional Imma cell",
            V0=V0,
            formula_units_per_cell=Z,
            atoms_per_formula=N,
            K0_gpa=K0,
            K0_prime=4.1,
            gamma0=1.64,
            q=1.5,
            theta0_k=814,
            reference_temperature_k=300,
        ),
        variants=variants,
        source_model_required_pressure_changes_by_temperature=per_t,
        cold_stage=dict(
            rows=int(rt.sum()),
            selection="all Table 1 rows at 304-316 K",
            k0_prime_when_assigned_300k=float(cold_fit.x[0]),
            k0_prime_with_actual_temperatures=float(actual_t_fit.x[0]),
        ),
        conditional_high_pressure_fits=dict(
            all_85=conditional_fit(x, t, p, 4.1),
            heated_74=conditional_fit(x[~rt], t[~rt], p[~rt], 4.1),
        ),
        normalization_fingerprint=dict(
            factor=factor,
            effective_n=7 * factor,
            note="One scalar projection at unchanged published coefficients. Not a physical atom count, recovered author setting, correction or registered refit.",
        ),
        implementation_checks=dict(
            native_available=native_available,
            native_maximum_pressure_difference_gpa=native_difference,
            python_maximum_pressure_difference_gpa=python_difference,
            independent_adaptive_energy_difference_j_mol=quadrature_difference,
            whole_cell_vs_formula_maximum_pressure_difference_gpa=float(
                np.max(abs(pressure(x, t, n=56, z=1) - baseline))
            ),
        ),
        later_author_code_crosscheck=author_check,
        figure2=dict(
            pressure_gpa=pf,
            temperatures_k=temperatures.tolist(),
            source_ratios=source_ratios.tolist(),
            graphical_bound=figure["reading_bound_volume_ratio"],
            comparisons=figure_results,
        ),
        ambient_source_calculation=ambient,
        pressure_calibration=dict(
            marker="MgO",
            source="Matsui, Parker and Leslie (2000), 10.2138/am-2000-2-308",
            basis="Inherited from the explicitly cited Katsura (2004) procedure, corroborated by the cited 2009 forsterite method; not independently stated in the short 2009 wadsleyite text.",
            paired_measured_mgo_cells_available=False,
            actual_recalibration_performed=False,
            limitation="No original MgO cells, detailed reducer settings or run identities recovered. Model-inverted cells would not be independent calibration observations.",
        ),
        missing_original_fit_inputs=[
            "Inoue (2004) ambient expansion observations and their selection/normalization",
            "Original weights, objective, covariance and fit session",
            "Numeric author atom/mole normalization and precise evaluation prescription",
            "Measured paired MgO cells and exact calibration reducer",
        ],
        disposition="Completed audit: source thermal-model/normalization consistency discrepancy is demonstrated; calibration alone cannot reconcile coefficients with the paper's own curves. Exact author-side cause unresolved. Keep the existing deferred hold and do not adopt n=5 or replace source coefficients.",
        row_residuals=[
            dict(
                row=i + 1,
                temperature_k=float(tt),
                pressure_gpa=float(pp),
                normalized_volume=float(xx),
                physical_n7_pressure_gpa=float(pred),
                residual_gpa=float(pred - pp),
                unphysical_n5_residual_gpa=float(pressure(xx, tt, n=5) - pp),
            )
            for i, (xx, tt, pp, pred) in enumerate(zip(x, t, p, baseline))
        ],
    )


def plot(result: dict, path: Path) -> None:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    d = load_data()
    p = d["pressure_gpa"]
    t = d["temperature_k"]
    x = d["normalized_volume"]
    fig, ax = plt.subplots(1, 3, figsize=(13.8, 4.2), layout="constrained")
    for n, color in [(7, "#b84a36"), (5, "#286a88")]:
        ax[0].scatter(
            t,
            pressure(x, t, n=n) - p,
            s=20,
            alpha=0.7,
            color=color,
            label="Physical n=7" if n == 7 else "Unphysical n=5 diagnostic",
        )
    ax[0].axhline(0, color="black", lw=0.6)
    ax[0].set(
        xlabel="Temperature (K)",
        ylabel="P(model) − P(source) (GPa)",
        title="Unchanged published coefficients",
    )
    ax[0].legend(fontsize=8)
    grid = np.linspace(10, 21, 50)
    for tt, col in [(500, "#666666"), (1300, "#436c92"), (2100, "#c37d29")]:
        for n, style in [(7, "--"), (5, "-")]:
            ratios = [
                volume_ratio(pp, tt, n=n) / volume_ratio(pp, 300, n=n) for pp in grid
            ]
            ax[1].plot(grid, ratios, style, color=col, label=f"{tt} K, n={n}")
        f = result["figure2"]
        index = f["temperatures_k"].index(tt)
        ax[1].errorbar(
            f["pressure_gpa"],
            f["source_ratios"][index],
            yerr=f["graphical_bound"],
            fmt="o",
            color=col,
            ms=4,
        )
    ax[1].set(
        xlabel="Pressure (GPa)",
        ylabel="V(P,T) / V(P,300 K)",
        title="Circles: digitized source curves",
    )
    ax[1].legend(fontsize=7, ncol=2)
    rt = t <= 320
    ax[2].errorbar(
        p[rt],
        x[rt] * V0,
        xerr=d["pressure_uncertainty_gpa"][rt],
        yerr=d["normalized_volume_uncertainty"][rt] * V0,
        fmt="o",
        ms=4,
        color="#333333",
        label="304–316 K source rows",
    )
    ax[2].plot(
        grid,
        [volume_ratio(pp, 300) * V0 for pp in grid],
        color="#286a88",
        label="Published 300 K BM3",
    )
    ax[2].set(
        xlabel="Pressure (GPa)", ylabel="Cell volume (Å³)", title="Cold-curve control"
    )
    ax[2].legend(fontsize=8)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(path, dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-pdf", type=Path)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    if args.source_pdf:
        verify_pdf(args.source_pdf)
    result = audit()
    if args.write:
        OUTPUT.write_text(json.dumps(result, indent=2) + "\n")
    if args.plot:
        plot(result, args.plot)
    print(
        json.dumps(
            {
                k: result[k]
                for k in [
                    "variants",
                    "cold_stage",
                    "normalization_fingerprint",
                    "implementation_checks",
                    "figure2",
                    "ambient_source_calculation",
                    "disposition",
                ]
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
