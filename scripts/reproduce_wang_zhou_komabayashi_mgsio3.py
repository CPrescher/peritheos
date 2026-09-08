"""Independent fits to the original published MgSiO3 observations and scales.

No Dorogokupets pressure rescaling, point weighting, or fitted synthetic data.
The source papers do not specify a fully reproducible residual/weighting rule;
unweighted pressure residuals are an explicit reconstruction choice, not a claim
about their original optimizer. Published parameters remain the catalog values.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro
from scipy.optimize import least_squares

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
REPORT = ROOT / "docs/data/wang-zhou-komabayashi-mgsio3-refit.json"
WANG = "akimotoite-wang-2004-tables1-2-pvt.csv"
ZHOU = "akimotoite-zhou-2014-table1-elasticity.csv"
KOMABAYASHI = "mgsio3-post-perovskite-komabayashi-2008-table1-pvt.csv"
IDS = (
    "akimotoite_wang_2004_t0133_bm3",
    "akimotoite_wang_2004_t0150_bm3",
    "akimotoite_zhou_2014_300k_bm3",
    "mgsio3_post_perovskite_komabayashi_2008_thermal_bm3",
)


def rows(filename):
    with (DATA / filename).open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def observations():
    """Explicit source masks; retain every selected row, including ambient rows."""
    wang = rows(WANG)
    result = {}
    for identifier, run in zip(IDS[:2], ("T0133", "T0150")):
        selected = [
            r for r in wang if r["run"] == run and float(r["temperature_k"]) == 298
        ]
        result[identifier] = [
            (
                int(r["source_order"]),
                float(r["unit_cell_volume_a3"]),
                float(r["pressure_original_gpa"]),
                298.0,
            )
            for r in selected
        ]
    result[IDS[2]] = [
        (
            int(r["source_order"]),
            6 * 100.387 * 1e24 / Avogadro / float(r["density_g_cm3"]),
            float(r["pressure_gpa"]),
            300.0,
        )
        for r in rows(ZHOU)
        if float(r["temperature_k"]) == 300
    ]
    result[IDS[3]] = [
        (
            int(r["source_order"]),
            float(r["ppv_volume_a3"]),
            float(r["pressure_mgo_gpa"]),
            float(r["temperature_k"]),
        )
        for r in rows(KOMABAYASHI)
        if r["ppv_fit_included"] == "true"
    ]
    return result


def solve(residual, initial):
    fit = least_squares(
        residual, initial, xtol=1e-12, ftol=1e-12, gtol=1e-12, max_nfev=10000
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    return fit


def thermal_model(k0, alpha1, dkdt):
    # Eq. 4: alpha(T)=alpha0+alpha1*T; alpha(300) is constrained, not alpha0.
    return ThermalReferenceStateEOS(
        BM3(163.813, k0, 4),
        300,
        1.7e-5 - 300 * alpha1,
        dkdt,
        alpha1,
        thermal_expansion_law="linear_temperature",
    )


def comparison(published, errors, fitted):
    return [
        dict(
            parameter=name,
            published=value,
            published_error=errors[name],
            refit=fitted[name],
            refit_error=None,
            difference=fitted[name] - value,
            relative_difference=abs(fitted[name] - value) / abs(value),
            within_reported_error=bool(abs(fitted[name] - value) <= errors[name]),
            within_combined_2sigma=None,
        )
        for name, value in published.items()
    ]


def reproduce():
    results = {}
    for identifier, data in observations().items():
        array = np.array(data)
        v, p, t = array[:, 1:].T
        stages = []
        sensitivity = None
        if identifier in IDS[:2]:
            published = (
                dict(V0=264.2, K0_prime=4.8)
                if identifier == IDS[0]
                else dict(V0=263.9, K0_prime=5.6)
            )
            errors = dict(V0=0.2, K0_prime=0.5 if identifier == IDS[0] else 0.8)
            fit = solve(lambda x: BM3(x[0], 210, x[1]).pressure(v) - p, [264, 5])
            fitted = dict(zip(published, map(float, fit.x)))
            published_model = BM3(published["V0"], 210, published["K0_prime"])
            predicted = published_model.pressure(v)
            fixed = dict(K0=210)
            selection = "All 298 K rows in " + (
                "Table 1, run T0133 (Au)"
                if identifier == IDS[0]
                else "Table 2, run T0150 (NaCl)"
            )
        elif identifier == IDS[2]:
            published, errors = dict(V0=262.45, K0=207), dict(V0=0.26, K0=3)
            fit = solve(lambda x: BM3(x[0], x[1], 4.6).pressure(v) - p, [262, 210])
            fitted = dict(zip(published, map(float, fit.x)))
            predicted = BM3(262.45, 207, 4.6).pressure(v)
            fixed = dict(K0_prime=4.6)
            selection = "All 17 akimotoite 300 K density rows, including three ambient density-only rows; V=6*M/(N_A*rho), M=100.387 g/mol"
        else:
            rt = t == 300
            room = solve(lambda x: BM3(163.813, x[0], 4).pressure(v[rt]) - p[rt], [220])
            k0 = float(room.x[0])

            def thermal_fit(k):
                return solve(
                    lambda x: (
                        thermal_model(k, x[0] * 1e-8, x[1] * 0.01).pressure(v, t) - p
                    ),
                    [0.1, -0.8],
                )

            fit = thermal_fit(k0)
            alpha1, dkdt = float(fit.x[0] * 1e-8), float(fit.x[1] * 0.01)
            published = dict(K0=223.2, alpha1=1.13e-9, dK_dT=-0.0085)
            errors = dict(K0=0.2, alpha1=1.14e-8, dK_dT=0.0011)
            fitted = dict(K0=k0, alpha1=alpha1, dK_dT=dkdt)
            fixed = dict(V0=163.813, K0_prime=4, Tr=300, alpha_at_300k=1.7e-5)
            stages = [
                dict(
                    stage="300 K BM3",
                    observations=int(rt.sum()),
                    free_parameters=["K0"],
                    fixed_parameters=["V0", "K0_prime"],
                    rmse_gpa=float(np.sqrt(np.mean(room.fun**2))),
                ),
                dict(
                    stage="P-V-T",
                    observations=len(p),
                    free_parameters=["alpha1", "dK_dT"],
                    fixed_parameters=["V0", "K0", "K0_prime", "Tr", "alpha_at_300k"],
                    dependent_parameters={"alpha0": "1.7e-5 - 300*alpha1"},
                ),
            ]
            rounded = thermal_fit(223.2)
            sensitivity = dict(
                description="Second stage with the published rounded K0=223.2 GPa fixed",
                alpha1=float(rounded.x[0] * 1e-8),
                dK_dT=float(rounded.x[1] * 0.01),
            )
            # Evaluate the actual printed coefficients, including rounded alpha0.
            predicted = ThermalReferenceStateEOS(
                BM3(163.813, 223.2, 4),
                300,
                1.667e-5,
                -0.0085,
                1.13e-9,
                thermal_expansion_law="linear_temperature",
            ).pressure(v, t)
            selection = "All 21 Table 1 rows with PPv volume, using original Speziale-MgO pressures; the Pv-only row is excluded"
        results[identifier] = dict(
            status="similar",
            observations=len(p),
            selection=selection,
            source_orders=[int(r[0]) for r in data],
            observed_pressure_range_gpa=[float(p.min()), float(p.max())],
            observed_temperature_range_k=[float(t.min()), float(t.max())],
            objective="unweighted pressure residuals",
            fixed_values=fixed,
            free_parameters=list(published),
            parameters=comparison(published, errors, fitted),
            rmse_gpa=float(np.sqrt(np.mean(fit.fun**2))),
            published_rmse_gpa=float(np.sqrt(np.mean((predicted - p) ** 2))),
            solver_success=bool(fit.success),
            stages=stages,
            qualification="All free coefficients agree within the printed source errors. Exact optimizer, weighting, covariance and confidence convention are not recovered; common-ledger classification is similar, not strict statistical parity.",
        )
        if sensitivity is not None:
            results[identifier].update(
                rounded_K0_sensitivity=sensitivity,
                derived_alpha0=1.7e-5 - 300 * fitted["alpha1"],
            )
    return dict(
        format="peritheos.wang-zhou-komabayashi-mgsio3-refit",
        format_version=1,
        generated_with="scripts/reproduce_wang_zhou_komabayashi_mgsio3.py",
        source_sha256={
            f: hashlib.sha256((DATA / f).read_bytes()).hexdigest()
            for f in (WANG, ZHOU, KOMABAYASHI)
        },
        density_conversion=dict(
            molar_mass_g_mol=100.387, formula_units=6, avogadro_mol_inverse=Avogadro
        ),
        fits=results,
    )


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    report = reproduce()
    if not all(
        p["within_reported_error"]
        for f in report["fits"].values()
        for p in f["parameters"]
    ):
        raise AssertionError("A recovered coefficient is outside its source error")
    # Retain full precision in the machine-readable diagnostics.
    text = json.dumps(report, indent=2, allow_nan=False) + "\n"
    if args.check:
        if json.loads(REPORT.read_text()) != report:
            raise SystemExit(f"Stale report: {REPORT}")
    else:
        REPORT.write_text(text, encoding="utf-8")
    print(json.dumps({k: v["status"] for k, v in report["fits"].items()}))


if __name__ == "__main__":
    main()
