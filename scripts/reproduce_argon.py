"""Independent primary-data diagnostics for fcc Ar and the legacy hcp fit.

Run with ``python -m scripts.reproduce_argon``. The equations below use
conventional-cell volumes and independent Gauss-Legendre Debye quadrature.
"""

from __future__ import annotations

import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.constants import Boltzmann
from scipy.optimize import brentq, least_squares

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/argon-reproduction.json"
_NODES, _WEIGHTS = np.polynomial.legendre.leggauss(96)
_NODES = (_NODES + 1) / 2
_WEIGHTS = _WEIGHTS / 2


def rows(filename):
    with (DATA / filename).open() as stream:
        return list(csv.DictReader(stream))


def vinet(volume, v0, k0, kp):
    x = (np.asarray(volume) / v0) ** (1 / 3)
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def dewaele(volume, temperature, v0=152.0, kp=7.423):
    v, t = np.broadcast_arrays(np.asarray(volume), np.asarray(temperature))
    x = v / v0
    gamma = 0.5 + 2.2 * x
    theta = 93.3 / np.sqrt(x) * np.exp(2.2 * (1 - x))
    y = theta / t
    # D3(y) = 3/y^3 integral_0^y z^3/(exp(z)-1) dz.
    z = y[..., None] * _NODES
    d3 = 3 * y * np.sum(_WEIGHTS * _NODES**3 / np.expm1(z), axis=-1)
    return vinet(v, v0, 2.65, kp) + 3 * gamma * 4 * Boltzmann * t * d3 / v * 1e21


def bm2(volume, v0=78.0, k0=6.5):
    ratio = v0 / np.asarray(volume)
    return 1.5 * k0 * (ratio ** (7 / 3) - ratio ** (5 / 3))


def rmse(residual):
    return float(np.sqrt(np.mean(np.asarray(residual) ** 2)))


def diagnostic_fit(model, volume, pressure, initial, bounds, sigma=None):
    weights = 1 if sigma is None else sigma
    result = least_squares(
        lambda pars: (model(volume, *pars) - pressure) / weights,
        initial,
        bounds=bounds,
        x_scale="jac",
        ftol=1e-12,
        xtol=1e-12,
        gtol=1e-12,
        max_nfev=3000,
    )
    return dict(
        parameters=result.x.tolist(),
        pressure_rmse_gpa=rmse(model(volume, *result.x) - pressure),
        solver_success=bool(result.success),
    )


@lru_cache(maxsize=1)
def reproduce():
    from peritheos import Material, get_material_document

    material = Material.from_eosmat(get_material_document("argon_fcc"))
    dr = rows("argon-fcc-dewaele-2021-supplement.csv")
    selected = [r for r in dr if r["selected_for_room_temperature_fit"] == "1"]
    dv = np.array([float(r["argon_a_angstrom"]) ** 3 for r in selected])
    dp = np.array([float(r["pressure_gpa"]) for r in selected])
    dmodel = material.get_eos_record("argon_fcc_dewaele_2021_vinet_mgd")
    dchecks = []
    for r in (
        selected[0],
        selected[30],
        max(selected, key=lambda r: float(r["pressure_gpa"])),
    ):
        v, p = float(r["argon_a_angstrom"]) ** 3, float(r["pressure_gpa"])
        pred = float(dewaele(v, 296))
        dchecks.append(
            dict(
                sample=r["sample"],
                volume_a3=v,
                temperature_k=296,
                source_pressure_gpa=p,
                calculated_pressure_gpa=pred,
                tolerance_gpa=max(0.15, 0.025 * p),
            )
        )
    low = [r for r in dr if r["run"] == "5" and r["pressure_gpa"]]
    lv = np.array([float(r["argon_a_angstrom"]) ** 3 for r in low])
    lt = np.array([float(r["temperature_k"]) for r in low])
    lp = np.array([float(r["pressure_gpa"]) for r in low])
    of = rows("argon-fcc-ono-2020-table1.csv")
    ov = np.array([float(r["volume_a3"]) for r in of])
    op = np.array([float(r["pressure_gpa"]) for r in of])
    os = np.array([float(r["pressure_sigma_gpa"]) for r in of])
    ovs = np.array([float(r["volume_sigma_a3"]) for r in of])
    opar = [184.5, 1.07, 8.02]
    deriv = (vinet(ov + 0.0001, *opar) - vinet(ov - 0.0001, *opar)) / 0.0002
    sigma = np.sqrt(os**2 + (deriv * ovs) ** 2)
    orc = material.get_eos_record("argon_fcc_ono_2020_vinet")
    wf = rows("argon-hcp-wittlinger-1997-figure3-digitized.csv")
    wv = np.array([float(r["volume_a3_conventional_cell"]) for r in wf])
    wp = np.array([float(r["pressure_gpa"]) for r in wf])
    ws = np.array(
        [
            (
                float(r["volume_error_minus_a3_conventional_cell"])
                + float(r["volume_error_plus_a3_conventional_cell"])
            )
            / 2
            for r in wf
        ]
    )
    wps = np.array([float(r["pressure_digitization_uncertainty_gpa"]) for r in wf])
    wvs = np.hypot(
        ws,
        78 * np.array([float(r["volume_ratio_digitization_uncertainty"]) for r in wf]),
    )
    wderiv = (bm2(wv + 0.0001) - bm2(wv - 0.0001)) / 0.0002
    weffective = np.hypot(wps, wderiv * wvs)
    # An errors-in-variables fit can move observed volumes. Report both raw
    # and adjusted-coordinate residuals so these cannot be confused.
    eiv = least_squares(
        lambda pars: np.r_[
            (bm2(pars[2:], *pars[:2]) - wp) / wps, (pars[2:] - wv) / wvs
        ],
        np.r_[78.0, 6.5, wv],
        bounds=(
            np.r_[50.0, 0.01, np.full(len(wv), 1)],
            np.r_[150.0, 100.0, np.full(len(wv), 150)],
        ),
        x_scale="jac",
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    comparison = []
    for pressure in [2.0, 5.0, 8.5, 20.0, 50.0, 100.0]:
        va = brentq(lambda v: float(dewaele(4 * v, 296)) - pressure, 3, 38)
        pw = float(bm2(2 * va))
        comparison.append(
            dict(
                fcc_pressure_gpa=pressure,
                atomic_volume_a3=va,
                wittlinger_pressure_gpa=pw,
                relative_pressure_difference=(pw - pressure) / pressure,
                within_wittlinger_pressure_range=1.2 <= pw <= 8.5,
            )
        )
    return dict(
        dewaele_2021=dict(
            total_rows=len(dr),
            selected_rows=len(selected),
            run_counts={
                str(n): sum(r["run"] == str(n) for r in dr) for n in range(1, 6)
            },
            published_rmse_gpa=rmse(dewaele(dv, 296) - dp),
            native_max_difference_gpa=float(
                np.max(np.abs(dmodel.pressure(dv, 296) - dewaele(dv, 296)))
            ),
            refits={
                str(t): diagnostic_fit(
                    lambda v, v0, kp: dewaele(v, t, v0, kp),
                    dv,
                    dp,
                    [152, 7.423],
                    ([100, 1], [220, 15]),
                )
                for t in (296, 300)
            },
            benchmarks=dchecks,
            low_temperature_rows=len(low),
            low_temperature_rmse_gpa=rmse(dewaele(lv, lt) - lp),
            low_temperature_native_max_difference_gpa=float(
                np.max(np.abs(dmodel.pressure(lv, lt) - dewaele(lv, lt)))
            ),
            qualification="95 selected rows; run 5 is independent validation. Equal-pressure weights are a diagnostic assumption: the publication supplies no objective, weights, coordinate uncertainties or covariance. Both 296 K (Methods) and 300 K (fit prose) evaluated. Missing pressures and duplicate rows retained in bundle.",
        ),
        ono_2020=dict(
            total_rows=len(of),
            published_rmse_gpa=rmse(vinet(ov, *opar) - op),
            native_max_difference_gpa=float(
                np.max(np.abs(orc.pressure(ov, 300) - vinet(ov, *opar)))
            ),
            refits={
                key: diagnostic_fit(
                    vinet, ov, op, opar, ([100, 0.001, 1], [400, 30, 20]), s
                )
                for key, s in [("equal_pressure", None), ("effective_errors", sigma)]
            },
            benchmarks=[
                dict(
                    volume_a3=float(ov[i]),
                    source_pressure_gpa=float(op[i]),
                    calculated_pressure_gpa=float(vinet(ov[i], *opar)),
                    tolerance_gpa=2.0,
                )
                for i in (0, 9, 18)
            ],
            qualification="All 19 rows. Effective-error weights combine printed 1-sigma P,V errors using a fixed derivative at published coefficients. Source weighting and covariance are unspecified; no strict parameter parity claimed.",
        ),
        wittlinger_1997=dict(
            reproduction_status="not_reproduced",
            total_rows=len(wf),
            published_parameters=dict(V0=78.0, K0=6.5),
            published_rmse_gpa=rmse(bm2(wv) - wp),
            published_max_abs_residual_gpa=float(np.max(np.abs(bm2(wv) - wp))),
            published_reduced_chi_square=float(
                np.sum(((bm2(wv) - wp) / weffective) ** 2) / (len(wv) - 2)
            ),
            equal_pressure_refit=diagnostic_fit(
                bm2, wv, wp, [78, 6.5], ([50, 0.01], [150, 100])
            ),
            errors_in_variables=dict(
                parameters=eiv.x[:2].tolist(),
                original_coordinate_rmse_gpa=rmse(bm2(wv, *eiv.x[:2]) - wp),
                adjusted_coordinate_rmse_gpa=rmse(bm2(eiv.x[2:], *eiv.x[:2]) - wp),
                reduced_chi_square=float(np.sum(eiv.fun**2) / (len(wv) - 2)),
            ),
            fcc_comparison=comparison,
            qualification="Comparison uses atomic volume: hcp cell/2 and fcc cell/4. Different phases are compared, not identified. Nine digitized points have large volume errors and share the fitted 78 A^3 normalization; unknown covariance limits coefficient inference. No replacement EOS is promoted from this diagnostic. Original article could not be reopened without institutional access on 2026-09-24; source transcription relies on the existing primary-source audit.",
        ),
    )


def ledger_outcome(record):
    if "wittlinger" in record["identifier"]:
        result = reproduce()["wittlinger_1997"]
        fit = result["errors_in_variables"]
        comparisons = [
            dict(
                parameter=name,
                published=published,
                refit=fitted,
                published_error=error,
                refit_error=None,
                relative_difference=abs(fitted - published) / published,
                within_combined_2sigma=None,
                similar=abs(fitted - published) <= 2 * error,
            )
            for name, published, fitted, error in zip(
                ["V0", "K0"], [78.0, 6.5], fit["parameters"], [3.0, 1.3]
            )
        ]
        return dict(
            status="not_refittable",
            reproduction_status="not_reproduced",
            dataset_identifiers=["argon_hcp_wittlinger_1997_figure3_digitized"],
            observations=9,
            selection="All nine experimental markers; fitted ambient marker excluded",
            fit_kind="digitized_bm2_errors_in_variables",
            objective="P and V residuals with digitization uncertainty combined with plotted volume error",
            free_parameters=["V0", "K0"],
            parameters=comparisons,
            rmse_gpa=fit["original_coordinate_rmse_gpa"],
            adjusted_coordinate_rmse_gpa=fit["adjusted_coordinate_rmse_gpa"],
            published_rmse_gpa=result["published_rmse_gpa"],
            reduced_chi_square=fit["reduced_chi_square"],
            solver_success=True,
            reason="NOT REPRODUCED: the available digitization does not establish a source-faithful reproduction of the published fit. Original observations, regression settings and covariance are unavailable; the shared fitted V0 normalization limits independent coefficient recovery. The numerical fits below are diagnostics only. Broad error-bar overlap is not successful reproduction. "
            + result["qualification"],
            reproduction=result,
        )
    is_dewaele = "dewaele" in record["identifier"]
    result = reproduce()["dewaele_2021" if is_dewaele else "ono_2020"]
    fitted = result["refits"]["296" if is_dewaele else "effective_errors"]
    names = ["V0", "K0_prime"] if is_dewaele else ["V0", "K0", "K0_prime"]
    comparisons = []
    for name, value in zip(names, fitted["parameters"]):
        pub = record["eos"]["parameters"][name]
        error = record["parameter_errors"].get(name)
        comparisons.append(
            dict(
                parameter=name,
                published=pub,
                refit=value,
                published_error=error,
                refit_error=None,
                relative_difference=abs(value - pub) / abs(pub),
                within_combined_2sigma=None,
                similar=abs(value - pub) <= max(error or 0, 0.05 * abs(pub)),
            )
        )
    return dict(
        status="similar"
        if all(p["similar"] for p in comparisons)
        else "parity_not_achieved",
        dataset_identifiers=record["fit_datasets"],
        observations=result.get("selected_rows", result["total_rows"]),
        selection="Runs 2-3 with measured pressure, plus runs 1 and 4 at P<=5 GPa"
        if is_dewaele
        else "All 19 Table 1 observations",
        fit_kind="primary_table_diagnostic",
        objective="equal pressure residuals"
        if is_dewaele
        else "fixed effective P,V errors",
        free_parameters=names,
        parameters=comparisons,
        rmse_gpa=fitted["pressure_rmse_gpa"],
        published_rmse_gpa=result["published_rmse_gpa"],
        solver_success=fitted["solver_success"],
        reason=result["qualification"],
        reproduction=result,
    )


if __name__ == "__main__":
    output = reproduce()
    REPORT.write_text(json.dumps(output, indent=2, sort_keys=True) + "\n")
    print(json.dumps(output, indent=2))
