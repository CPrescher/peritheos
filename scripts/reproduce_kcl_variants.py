"""Independent primary-table reproduction of nondefault Walker/Tateno/Chidester KCl."""

from __future__ import annotations

import argparse
import csv
import json
from functools import lru_cache
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares
from scipy.special import roots_legendre

from peritheos import get_eos_record, get_eos_record_document

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/kcl-variants-reproduction.json"
A3_PER_MOLAR = 1e24 / 6.02214076e23
GAS_CONSTANT = 0.00831446261815324  # GPa cm3 mol-1 K-1
IDS = (
    "kcl_b1_walker_2002_bm3_linear_thermal",
    "kcl_b2_tateno_2019_holmes_vinet_mgd",
    "kcl_b2_tateno_2019_holmes_vinet_linear_thermal",
    "kcl_b2_tateno_2019_sokolova_vinet_linear_thermal",
    "kcl_b2_chidester_2021_vinet_mgd",
)


def rows(filename):
    with (ROOT / "peritheos/data/datasets" / filename).open() as stream:
        return list(csv.DictReader(stream))


def column(data, name):
    return np.array([float(row[name]) for row in data])


def vinet(volume, v0, k0, kp):
    """Primary Vinet expression; no Peritheos evaluator."""
    x = (volume / v0) ** (1 / 3)
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def bm3(volume, v0, k0, kp):
    """Compression-positive Eulerian finite strain."""
    f = ((v0 / volume) ** (2 / 3) - 1) / 2
    return 3 * k0 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * f * (kp - 4))


@lru_cache(maxsize=2)
def quadrature(order):
    x, weights = roots_legendre(order)
    return (x + 1) / 2, weights / 2


def mgd(volume, temperature, parameters, order=64):
    """Vinet+MGD in molar cm3/mol; directly integrate the Debye energy."""
    v0, k0, kp, gamma0, q = parameters
    gamma = gamma0 * (volume / v0) ** q
    theta = 235 * np.exp((gamma0 - gamma) / q)
    x, weights = quadrature(order)

    def energy(t):
        y = theta / t
        u = y[:, None] * x
        with np.errstate(over="ignore"):
            debye = 3 * np.sum(weights * u**3 / np.expm1(u), axis=1) / y**2
        return 6 * GAS_CONSTANT * t * debye  # two atoms per KCl formula unit

    return (
        vinet(volume, v0, k0, kp)
        + gamma
        * (energy(temperature) - energy(np.full_like(temperature, 300)))
        / volume
    )


def metrics(prediction, observed):
    residual = prediction - observed
    return {
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_absolute_residual_gpa": float(np.max(np.abs(residual))),
    }


@lru_cache(maxsize=1)
def reproduce():
    output = {}
    walker_all = rows("kcl-walker-2002-table1-pvt.csv")
    walker = [r for r in walker_all if r["included_in_fit"] == "1"]
    v = column(walker, "b1_kcl_cell_volume_a3")
    p = column(walker, "pressure_kbar") * 0.1
    t = column(walker, "temperature_celsius") + 273.15
    shape = bm3(v, 249.53, 1, 5)
    # Walker specifies ordinary pressure least squares. Joint K0 and alpha*K0
    # is the reproducible complete-data diagnostic; B1 staging is not explicit.
    fitted = np.linalg.lstsq(np.column_stack((shape, t - 296.15)), p, rcond=None)[0]
    prediction = shape * 17.7 + 0.00195 * (t - 296.15)
    native = get_eos_record(IDS[0]).pressure(v, temperature=t, check_validity=False)
    output[IDS[0]] = {
        "total_source_rows": len(walker_all),
        "observations": len(walker),
        "selection": "23 measured Table 1 KCl states; exclude six NaCl-only checks and one derived V0 anchor",
        "free_parameters": ["K0", "alpha_KT"],
        "fitted_parameters": dict(zip(("K0", "alpha_KT"), fitted.tolist())),
        "published": metrics(prediction, p),
        "refit": metrics(shape * fitted[0] + fitted[1] * (t - 296.15), p),
        "native_max_difference_gpa": float(np.max(np.abs(native - prediction))),
        "qualification": "Fixed V0=249.53 conventional-cell A3 (Z=4), K0 prime=5, Tr=296.15 K. Complete unweighted joint pressure fit; original B1 staging is not fully explicit. K0 differs by 2.94%, alpha_KT by 3.99%; strict coefficient/uncertainty parity is not established. Individual elastic errors and covariance are not published. The printed BE1 signs are inconsistent with its positive strain definition; the compression-positive BM3 convention used by the existing B2 audit is retained explicitly.",
    }
    official = rows("kcl-tateno-2019-official-table-s1.csv")
    v = column(official, "kcl_unit_cell_volume_a3") / A3_PER_MOLAR
    vp = column(official, "platinum_unit_cell_volume_a3")
    t = column(official, "temperature_k")
    sokolova_p = column(official, "pressure_gpa")
    # Holmes (1989) Equations 11-12: exact theoretical curve mapping and the
    # rounded thermal B_T. These are not fitted against KCl observations.
    holmes_p = vinet(vp, 60.4000884, 798.31 / 3, 1 + 7.2119 / 1.5) + 0.0069426 * (
        t - 300
    )
    for identifier, pressure in zip(IDS[1:4], (holmes_p, holmes_p, sokolova_p)):
        record = get_eos_record_document(identifier)
        static = record["eos"]["parameters"]
        thermal = record["thermal"]["parameters"]
        is_mgd = "gamma0" in thermal
        initial = [static["K0"], static["K0_prime"]] + (
            [thermal["gamma0"], thermal["q"]] if is_mgd else [thermal["alpha_KT"]]
        )
        names = ["K0", "K0_prime"] + (["gamma0", "q"] if is_mgd else ["alpha_KT"])

        def model(z, order=64):
            if is_mgd:
                return mgd(v, t, [54.5 / A3_PER_MOLAR, *z], order)
            return vinet(v, 54.5 / A3_PER_MOLAR, z[0], z[1]) + z[2] * (t - 300)

        fit = least_squares(
            lambda z: model(z) - pressure,
            initial,
            bounds=([1, 2, 0.1, 0.01], [50, 10, 8, 4]) if is_mgd else (-np.inf, np.inf),
            x_scale="jac",
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        pred = model(initial)
        native = get_eos_record(identifier).pressure(
            v * A3_PER_MOLAR, temperature=t, check_validity=False
        )
        output[identifier] = {
            "observations": len(official),
            "free_parameters": names,
            "fitted_parameters": dict(zip(names, fit.x.tolist())),
            "selection": "all 39 correctly paired official Table S1 rows; no exclusions",
            "published": metrics(pred, pressure),
            "refit": metrics(model(fit.x), pressure),
            "solver_success": bool(fit.success),
            "native_max_difference_gpa": float(np.max(np.abs(native - pred))),
            "quadrature_max_difference_gpa": float(
                np.max(np.abs(model(initial, 96) - pred))
            )
            if is_mgd
            else 0,
            "pressure_coordinate": "conditional Holmes Equations 11-12 derived from Pt volumes"
            if identifier != IDS[3]
            else "author-reported Sokolova Pt pressures",
            "qualification": (
                "Equal pressure weights; covariance and original objective details unavailable. "
                "Holmes-coordinate results are conditional: exact author rounding/thermal reduction "
                "is unspecified and Holmes's approximation is stated for T<2000 K. "
                "No pressure-coordinate or statistical parity claimed; source parameter errors "
                "are widths of unknown confidence. Published coefficients remain unchanged."
                if identifier != IDS[3]
                else "All 39 author-reported Sokolova-pressure observations. Equal pressure weights; "
                "covariance and original objective details unavailable. Coefficients agree within "
                "printed error widths of unknown confidence; no statistical parity claimed. "
                "Published coefficients remain unchanged."
            ),
        }
    output["holmes_pressure_coordinate"] = {
        "source": "Holmes (1989) Equations 11-12; paired official Tateno Pt volumes",
        "reference_eos_record": "platinum_holmes_1989_vinet_1",
        "calibration_parameter_errors": None,
        "rows": [
            {
                "source_excel_row": int(row["source_excel_row"]),
                "pressure_gpa": float(pressure),
                "outside_holmes_thermal_approximation": bool(temp >= 2000),
            }
            for row, pressure, temp in zip(official, holmes_p, t)
        ],
        "qualification": "Derived model coordinate, not published Holmes pressures or independent observations. Raw author Sokolova pressures, Pt volumes and coordinate errors remain in the official dataset. Neither Sokolova pressure errors nor missing calibration covariance are relabeled as Holmes errors.",
    }
    cold = rows("kcl-dewaele-2012-table1-compression.csv")
    hot = rows("chidester_2021_kcl_pvt.csv")
    v = np.r_[
        column(cold, "volume_a3_per_formula_unit") / A3_PER_MOLAR,
        column(hot, "V (cc/mole)"),
    ]
    p = np.r_[column(cold, "pressure_gpa"), column(hot, "P (GPa)")]
    t = np.r_[np.full(len(cold), 300), column(hot, "T KCl (K)")]
    initial = [34.3, 13, 6.2, 3.4, 1.0]
    fit = least_squares(
        lambda z: mgd(v, t, z) - p,
        initial,
        bounds=([25, 1, 2, 0.1, 0.01], [45, 50, 10, 8, 4]),
        x_scale="jac",
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    names = ["V0", "K0", "K0_prime", "gamma0", "q"]
    fitted = fit.x.copy()
    fitted[0] *= A3_PER_MOLAR
    pred = mgd(v, t, initial)
    native = get_eos_record(IDS[4]).pressure(
        v * A3_PER_MOLAR, temperature=t, check_validity=False
    )
    output[IDS[4]] = {
        "observations": len(p),
        "room_temperature_rows": len(cold),
        "effective_temperature_rows": len(hot),
        "free_parameters": names,
        "fitted_parameters": dict(zip(names, fitted.tolist())),
        "selection": "all 123 Dewaele room-temperature rows as 300 K source reference plus all 155 Chidester effective-temperature rows",
        "published": metrics(pred, p),
        "refit": metrics(mgd(v, t, fit.x), p),
        "high_temperature_refit": metrics(
            mgd(v, t, fit.x)[len(cold) :], p[len(cold) :]
        ),
        "solver_success": bool(fit.success),
        "native_max_difference_gpa": float(np.max(np.abs(native - pred))),
        "quadrature_max_difference_gpa": float(
            np.max(np.abs(mgd(v, t, initial, 96) - pred))
        ),
        "qualification": "All five Table I coefficients are fitted, including q. Fixed theta0=235 K, Tr=300 K and n=2. Equal pressure weights; no original covariance or explicit row weights. Integrated-Gruneisen theta law is inferred by complete-data numerical reproduction; Equations 3-4 do not explicitly define theta(V). Dewaele 298 K measurements enter the source regression as its 300 K reference, preserving raw temperature provenance. High-temperature KCl temperatures are geometry-conditioned modeled effective temperatures. Cold ruby and high-T Pt calibration ancestry remain distinct.",
    }
    return output


def ledger_outcome(record):
    result = reproduce()[record["identifier"]]
    comparisons = []
    for name in result["free_parameters"]:
        block = record if name in record["eos"]["parameters"] else record["thermal"]
        parameters = (
            record["eos"]["parameters"] if block is record else block["parameters"]
        )
        published, fitted = parameters[name], result["fitted_parameters"][name]
        error = block["parameter_errors"].get(name)
        # Missing elastic errors cannot prove uncertainty parity. Five percent
        # is a declared numerical similarity criterion for the Walker diagnostic.
        similar = (error is not None and abs(fitted - published) <= error) or (
            abs(fitted / published - 1) <= 0.05
        )
        comparisons.append(
            {
                "parameter": name,
                "published": published,
                "published_error": error,
                "refit": fitted,
                "refit_error": None,
                "difference": fitted - published,
                "relative_difference": (fitted - published) / published,
                "within_combined_2sigma": None,
                "similar": similar,
            }
        )
    return {
        "status": "similar"
        if all(c["similar"] for c in comparisons)
        else "parity_not_achieved",
        "dataset_identifiers": record["fit_datasets"],
        "observations": result["observations"],
        "selection": result["selection"],
        "fit_kind": "independent_primary_equation_pressure_fit",
        "objective": "equal pressure weights",
        "free_parameters": result["free_parameters"],
        "parameters": comparisons,
        "rmse_gpa": result["refit"]["rmse_gpa"],
        "published_rmse_gpa": result["published"]["rmse_gpa"],
        "solver_success": result.get("solver_success", True),
        "reason": result["qualification"],
        "reproduction": result,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    assert all(
        result[identifier]["native_max_difference_gpa"] < 1e-9 for identifier in IDS
    )
    if args.check:
        saved = json.loads(OUTPUT.read_text())

        def check(actual, expected, path="report"):
            if isinstance(actual, dict):
                assert actual.keys() == expected.keys(), path
                for key in actual:
                    check(actual[key], expected[key], f"{path}.{key}")
            elif isinstance(actual, list):
                assert len(actual) == len(expected), path
                for index, (a, e) in enumerate(zip(actual, expected)):
                    check(a, e, f"{path}[{index}]")
            elif isinstance(actual, float):
                assert np.isclose(actual, expected, rtol=1e-6, atol=1e-9), path
            else:
                assert actual == expected, path

        check(result, saved)
    else:
        OUTPUT.write_text(
            json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        )
    print(json.dumps({i: result[i]["refit"] for i in IDS}, indent=2))


if __name__ == "__main__":
    main()
