"""Cross-check actual EosFit q-compromise runs using Peritheos components.

The audit-local model composes native BM3 and native Debye energy at the
reference volume. It is not added as a new catalog or core library model.
"""

from __future__ import annotations

import argparse
import hashlib
import json

import numpy as np

from peritheos.eos import ThermalEOS
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_joint_eos
from scripts.compare_fes_peritheos_eosfit import (
    PUBLISHED,
    ROOT,
    input_coordinates,
    statistics,
)

DIRECT = ROOT / "docs/data/fes-eosfit7c-q-compromise"
STAGED = ROOT / "docs/data/fes-eosfit7c-staged-q-compromise"
OUTPUT = ROOT / "docs/data/fes-q-compromise-peritheos-comparison.json"


class QCompromiseDebye(ThermalEOS):
    """Audit composition: Pth=gamma0/V0 * [E(theta0,T)-E(theta0,Tr)]."""

    def __init__(self, rt_eos, Tr, theta0, gamma0, n):
        super().__init__(rt_eos)
        # Energy evaluated at V0 has theta=theta0. q=0 is used solely inside
        # this energy helper; the composed q-compromise model has no q.
        self.energy_model = MieGruneisenDebye(rt_eos, Tr, theta0, gamma0, 0, n)
        self.Tr, self.theta0 = self.energy_model.Tr, self.energy_model.theta0
        self.gamma0, self.n = self.energy_model.gamma0, self.energy_model.n

    def thermal_pressure(self, V, T):
        volumes, temperatures = self._broadcast_state(V, T)
        v0 = self.rt_eos.V0
        energy = self.energy_model.thermal_energy(v0, temperatures)
        reference = self.energy_model.thermal_energy(v0, self.Tr)
        result = self.gamma0 * (energy - reference) / v0 / 10000
        return self._scalar_or_array(np.broadcast_to(result, volumes.shape))


def model(parameters):
    return QCompromiseDebye(
        BM3(parameters["V0"] / 10, parameters["K0"], parameters["Kp"]),
        300,
        417,
        parameters["Gamm0"],
        2,
    )


def compare(parameters, coordinates, reference):
    observed, volume, temperature = coordinates
    eos = model(parameters)
    if not hasattr(eos.rt_eos, "_native") or not hasattr(eos.energy_model, "_native"):
        raise RuntimeError("Peritheos native components required for the audit")
    predicted = eos.pressure(volume, temperature)
    residual = predicted - observed
    old = np.array([row["residual_gpa"] for row in reference["rows"]])
    if len(old) != len(residual):
        raise ValueError("EosFit and Peritheos row counts differ")
    del eos.energy_model._native
    fallback = eos.pressure(volume, temperature)
    return {
        "parameters_cm3_mol_gpa": parameters,
        "eosfit_all": statistics(old),
        "peritheos_all": statistics(residual),
        "eosfit_thermal": statistics(old[:146]),
        "peritheos_thermal": statistics(residual[:146]),
        "max_pressure_difference_gpa": float(np.max(np.abs(residual - old))),
        "max_native_python_energy_difference_in_pressure_gpa": float(
            np.max(np.abs(predicted - fallback))
        ),
        "rows": [
            {
                "row": i + 1,
                "peritheos_pressure_gpa": float(predicted[i]),
                "eosfit_pressure_from_printed_residual_gpa": float(
                    observed[i] + old[i]
                ),
                "difference_gpa": float(residual[i] - old[i]),
            }
            for i in range(len(observed))
        ],
    }


def reproduce():
    direct = json.loads((DIRECT / "manifest.json").read_text())
    staged = json.loads((STAGED / "manifest.json").read_text())
    cases, files = {}, [DIRECT / "manifest.json", STAGED / "manifest.json"]
    for name, case in direct["cases"].items():
        if name.startswith("cold"):
            continue
        folder = DIRECT / name
        files.append(folder / "input.dat")
        coordinates = input_coordinates(folder / "input.dat")
        fitted = PUBLISHED.copy()
        fitted.update(
            {
                key: p["value"]
                for key, p in case["refined_parameters_final_cycle"].items()
            }
        )
        cases[name] = {
            "initial": compare(
                PUBLISHED.copy(), coordinates, case["initial_source_replay"]
            ),
            "fitted": compare(fitted, coordinates, case["fitted_replay"]),
        }
    for name, case in staged["cases"].items():
        path = STAGED / name / "all.dat"
        files.append(path)
        coordinates = input_coordinates(path)
        initial = {
            key: p["value"]
            for key, p in case["cold_stage"]["refined_parameters_final_cycle"].items()
        }
        initial["Gamm0"] = 2.42
        fitted = initial.copy()
        fitted["Gamm0"] = case["thermal_stage"]["refined_parameters_final_cycle"][
            "Gamm0"
        ]["value"]
        cases[f"staged-{name}"] = {
            "initial": compare(
                initial, coordinates, case["thermal_stage"]["initial_source_replay"]
            ),
            "fitted": compare(
                fitted, coordinates, case["thermal_stage"]["fitted_replay"]
            ),
        }

    observed, volume, temperature = input_coordinates(
        DIRECT / "combined-unit/input.dat"
    )
    fit = fit_joint_eos(
        QCompromiseDebye,
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
        fixed={"Tr": 300.0, "theta0": 417.0, "n": 2.0},
        bounds={
            "rt_eos.V0": (0.8, 2.5),
            "rt_eos.K0": (1, 500),
            "rt_eos.K0_prime": (0, 15),
            "gamma0": (0.01, 10),
        },
        max_nfev=1000,
    )
    mapping = {
        "V0": "rt_eos.V0",
        "K0": "rt_eos.K0",
        "Kp": "rt_eos.K0_prime",
        "Gamm0": "gamma0",
    }
    if (
        not fit.success
        or not np.all(np.isfinite(fit.covariance))
        or tuple(mapping.values()) != fit.free_parameters
    ):
        raise RuntimeError("Independent q-compromise joint fit failed")
    scales = np.array([10.0, 1.0, 1.0, 1.0])
    parameters = {
        key: {
            "value": float(fit.parameters[name] * scale),
            "esd": float(fit.standard_errors[name] * scale),
        }
        for (key, name), scale in zip(mapping.items(), scales)
    }
    residual = fit.model.pressure(volume, temperature) - observed
    return {
        "scope": "Actual EosFit7c q-compromise replay and refits, forward-checked with native Peritheos BM3 and Debye-energy components; independent Peritheos API equal-weight joint refit using an audit-local composition.",
        "qualification": "GUI estimation uses q-compromise, but GUI usage alone does not prove the final author setting. Same 146 hot and 13 re-reported cold VI points; incomplete original input recovery. Theta0=417 K, n=2, Tr=300 K fixed; q undefined. Weighted EosFit fits are evaluated, not independently refitted with a different objective. No catalog coefficients changed or physical validation claimed.",
        "input_sha256": {
            str(p.relative_to(ROOT)): hashlib.sha256(p.read_bytes()).hexdigest()
            for p in files
        },
        "cases": cases,
        "independent_joint_equal_weight_refit": {
            "parameters_cm3_mol_gpa": parameters,
            "covariance_parameter_order": list(mapping),
            "covariance": (fit.covariance * np.outer(scales, scales)).tolist(),
            "solver": fit.to_dict()["solver"],
            "degrees_of_freedom": fit.degrees_of_freedom,
            "all_rows": statistics(residual),
            "thermal_rows": statistics(residual[:146]),
        },
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    text = json.dumps(result, indent=2, allow_nan=False) + "\n"
    if args.check:
        if OUTPUT.read_text() != text:
            raise SystemExit("Q-compromise comparison is stale")
    else:
        OUTPUT.write_text(text)
    print(
        "Max difference GPa:",
        max(
            s["max_pressure_difference_gpa"]
            for c in result["cases"].values()
            for s in c.values()
        ),
    )
    print(
        "Independent joint parameters:",
        result["independent_joint_equal_weight_refit"]["parameters_cm3_mol_gpa"],
    )


if __name__ == "__main__":
    main()
