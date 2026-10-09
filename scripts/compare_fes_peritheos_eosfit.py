"""Compare native Peritheos pressures and a joint refit with retained EosFit runs.

This uses the public Peritheos EOS and fitting APIs, not the independent audit
formula. No source coefficients, observations or catalog records are changed.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.constants import Avogadro

import peritheos
from peritheos import Material, get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_joint_eos

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/fes-peritheos-eosfit-comparison.json"
DIRECT = ROOT / "docs/data/fes-eosfit7c"
STAGED = ROOT / "docs/data/fes-eosfit7c-staged"
PUBLISHED = {"V0": 15.4, "K0": 115.5, "Kp": 4.99, "Gamm0": 2.42}


def model(parameters, *, gas_constant=None):
    """Peritheos uses J/bar/mol: divide EosFit's cm3/mol volumes by ten."""
    config = {} if gas_constant is None else {"Cvmax": 6 * gas_constant}
    return MieGruneisenDebye(
        BM3(parameters["V0"] / 10, parameters["K0"], parameters["Kp"]),
        Tr=300,
        theta0=417,
        gamma0=parameters["Gamm0"],
        q=1,
        n=2,
        debye_temperature_law="integrated_gruneisen",
        thermal_pressure_reference="reference_temperature",
        **config,
    )


def input_coordinates(path, *, cell_units=False):
    data = np.loadtxt(path, skiprows=5)
    factor = Avogadro / 4e24 if cell_units else 1
    return data[:, 0], data[:, 2] * factor / 10, data[:, 4]


def statistics(residual):
    return {
        "rows": len(residual),
        "rmse_gpa": float(np.sqrt(np.mean(residual**2))),
        "max_abs_gpa": float(np.max(np.abs(residual))),
        "rows_exceeding_3_gpa": (np.flatnonzero(np.abs(residual) > 3) + 1).tolist(),
    }


def comparison(parameters, coordinates, reference):
    observed, volumes, temperatures = coordinates
    eos = model(parameters)
    if not hasattr(eos, "_native"):
        raise RuntimeError("This audit requires the actual Peritheos native evaluator")
    predicted = np.asarray(eos.pressure(volumes, temperatures))
    residual = predicted - observed
    eosfit_residual = np.array([row["residual_gpa"] for row in reference["rows"]])
    if len(residual) != len(eosfit_residual):
        raise ValueError("Observation count differs from the EosFit output")
    # Cross-check Peritheos's thermal Python fallback at the same coefficients.
    del eos._native
    python_pressure = np.asarray(eos.pressure(volumes, temperatures))
    difference = residual - eosfit_residual
    # Diagnostic only: retain the default-constant result as the actual replay.
    rounded_r_pressure = model(parameters, gas_constant=8.314).pressure(
        volumes, temperatures
    )
    result = {
        "parameters_cm3_mol_gpa": parameters,
        "peritheos": statistics(residual),
        "eosfit": statistics(eosfit_residual),
        "max_pressure_difference_gpa": float(np.max(np.abs(difference))),
        "max_native_python_difference_gpa": float(
            np.max(np.abs(predicted - python_pressure))
        ),
        "max_difference_with_R_8_314_diagnostic_gpa": float(
            np.max(np.abs(rounded_r_pressure - observed - eosfit_residual))
        ),
        "rows": [
            {
                "row": i + 1,
                "observed_pressure_gpa": float(observed[i]),
                "peritheos_pressure_gpa": float(predicted[i]),
                "eosfit_pressure_from_printed_residual_gpa": float(
                    observed[i] + eosfit_residual[i]
                ),
                "difference_gpa": float(difference[i]),
            }
            for i in range(len(observed))
        ],
    }
    if len(observed) >= 146:
        result["peritheos_thermal"] = statistics(residual[:146])
        result["eosfit_thermal"] = statistics(eosfit_residual[:146])
    return result


def reproduce():
    direct = json.loads((DIRECT / "manifest.json").read_text(encoding="utf-8"))
    staged = json.loads((STAGED / "manifest.json").read_text(encoding="utf-8"))
    cases = {}
    inputs = [DIRECT / "manifest.json", STAGED / "manifest.json"]
    for name, case in direct["cases"].items():
        path = DIRECT / name / "input.dat"
        inputs.append(path)
        cold = name.startswith("cold")
        coordinates = input_coordinates(path, cell_units=cold)
        initial = (
            {"V0": 98.96 * Avogadro / 4e24, "K0": 148, "Kp": 4.53, "Gamm0": 0}
            if cold
            else PUBLISHED.copy()
        )
        fitted = initial.copy()
        fitted.update(
            {
                key: item["value"]
                for key, item in case["refined_parameters_final_cycle"].items()
            }
        )
        if cold:
            fitted["V0"] *= Avogadro / 4e24
        cases[name] = {
            "initial": comparison(initial, coordinates, case["initial_source_replay"]),
            "fitted": comparison(fitted, coordinates, case["fitted_replay"]),
        }
    for name, case in staged["cases"].items():
        path = STAGED / name / "all.dat"
        inputs.append(path)
        coordinates = input_coordinates(path)
        initial = {
            key: item["value"]
            for key, item in case["cold_stage"][
                "refined_parameters_final_cycle"
            ].items()
        }
        initial["Gamm0"] = 2.42
        fitted = initial.copy()
        fitted["Gamm0"] = case["thermal_stage"]["refined_parameters_final_cycle"][
            "Gamm0"
        ]["value"]
        cases[f"staged-{name}"] = {
            "initial": comparison(
                initial, coordinates, case["thermal_stage"]["initial_source_replay"]
            ),
            "fitted": comparison(
                fitted, coordinates, case["thermal_stage"]["fitted_replay"]
            ),
        }

    # Independent four-parameter joint fit through the public fitting API.
    observed, volume, temperature = input_coordinates(
        DIRECT / "combined-unit/input.dat"
    )
    fit = fit_joint_eos(
        MieGruneisenDebye,
        BM3,
        volume,
        temperature,
        observed,
        initial={
            "rt_eos.V0": 1.54,
            "rt_eos.K0": 115.5,
            "rt_eos.K0_prime": 4.99,
            "gamma0": 2.42,
        },
        fixed={"Tr": 300.0, "theta0": 417.0, "q": 1.0, "n": 2.0},
        bounds={
            "rt_eos.V0": (0.8, 2.5),
            "rt_eos.K0": (1, 500),
            "rt_eos.K0_prime": (0, 15),
            "gamma0": (0.01, 10),
        },
        max_nfev=1000,
    )
    if not fit.success or not np.all(np.isfinite(fit.covariance)):
        raise RuntimeError(
            "Peritheos joint fit did not converge with finite covariance"
        )
    mapping = {
        "V0": "rt_eos.V0",
        "K0": "rt_eos.K0",
        "Kp": "rt_eos.K0_prime",
        "Gamm0": "gamma0",
    }
    if tuple(mapping.values()) != fit.free_parameters:
        raise RuntimeError("Unexpected fit covariance parameter order")
    units = np.array([10.0, 1.0, 1.0, 1.0])
    parameters = {
        key: {
            "value": fit.parameters[name] * scale,
            "esd": fit.standard_errors[name] * scale,
        }
        for (key, name), scale in zip(mapping.items(), units)
    }
    refit_residual = fit.model.pressure(volume, temperature) - observed
    # Also exercise the catalog adapter's conventional-cell conversion.
    record = Material.from_eosmat(
        get_material_document("fes_vi"), require_primary_validation=False
    ).get_eos_record("fes_vi_morard_2026_bm3_mgd")
    catalog_pressure = record.pressure(
        volume[:146] * 10 * 4e24 / Avogadro, temperature[:146], check_validity=False
    )
    api_pressure = model(PUBLISHED).pressure(volume[:146], temperature[:146])
    from peritheos import _rust

    return {
        "scope": "Native Peritheos forward evaluation of eight actual EosFit cases, plus an independent equal-weight four-parameter joint refit. Weighted EosFit solutions are evaluated, not refitted with a different objective.",
        "qualification": "Same retained input coordinates; 13 cold re-reported VI rows, not confirmed complete author inputs. EosFit coefficients and residuals have limited printed precision. Default Peritheos gas constant retained; R=8.314 is only a separate diagnostic. No calibration covariance or physical validation.",
        "peritheos_version": peritheos.__version__,
        "native_extension_sha256": hashlib.sha256(
            Path(_rust.__file__).read_bytes()
        ).hexdigest(),
        "input_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in inputs
        },
        "volume_conversion": "EosFit cm3/mol divided by 10 gives Peritheos J/bar/mol; conventional cell Z=4 converted using exact SI Avogadro constant.",
        "gas_constant_diagnostic": "Peritheos defaults to R=8.31446261815324 J/mol/K. A separate Cvmax=6*8.314 diagnostic tests the rounded constant found in the audited CrysFML source. It does not change default evaluations or prove the exact constants used in the distributed executable. Remaining differences also involve printed coefficients/residuals and implementation precision.",
        "cases": cases,
        "catalog_adapter_max_pressure_difference_gpa": float(
            np.max(np.abs(catalog_pressure - api_pressure))
        ),
        "independent_joint_equal_weight_refit": {
            "parameters_cm3_mol_gpa": parameters,
            "eosfit_parameters_cm3_mol_gpa": direct["cases"]["combined-unit"][
                "refined_parameters_final_cycle"
            ],
            "covariance_parameter_order": list(mapping),
            "covariance": (fit.covariance * np.outer(units, units)).tolist(),
            "all_rows": statistics(refit_residual),
            "thermal_rows": statistics(refit_residual[:146]),
            "solver": fit.to_dict()["solver"],
            "degrees_of_freedom": fit.degrees_of_freedom,
            "qualification": "Conditional residual-scaled parameter errors; same unweighted pressure objective. q, theta0, n and Tr fixed. Peritheos weighted errors-in-variables fitting has a different objective from EosFit iterative weighting and is not claimed equivalent here.",
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    serialized = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.check:
        if OUTPUT.read_text(encoding="utf-8") != serialized:
            raise SystemExit("Peritheos/EosFit comparison is stale")
    else:
        OUTPUT.write_text(serialized, encoding="utf-8")
    print(
        "Maximum pressure difference (GPa):",
        max(
            c[s]["max_pressure_difference_gpa"]
            for c in result["cases"].values()
            for s in ["initial", "fitted"]
        ),
    )
    print(
        "Joint Peritheos fit:",
        result["independent_joint_equal_weight_refit"]["parameters_cm3_mol_gpa"],
    )


if __name__ == "__main__":
    main()
