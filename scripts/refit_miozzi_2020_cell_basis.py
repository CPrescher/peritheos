"""Fit Miozzi's data with Peritheos on per-Fe and whole-hcp-cell bases.

Run: python -m scripts.refit_miozzi_2020_cell_basis
This diagnostic does not register a new EOS or alter source observations.
Both normalizations use the user-relayed three-stage fitting sequence, with
theta0=420 K, Tr=300 K and equal pressure-residual weights. Original author
weights and pressure inputs remain unresolved.
"""

from __future__ import annotations

import hashlib
import json
import platform
from pathlib import Path

import numpy as np

from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_joint_eos, fit_rt_eos, fit_thermal_eos
from scripts.reconstruct_miozzi_2020_iron import (
    NA,
    ROOT,
    R,
    bm3,
    constant_q_gamma_theta,
    debye_energy,
    pressure_column,
    read_data,
)
from scripts.refit_miozzi_2020_tange import target_pressures

OUTPUT = ROOT / "docs/data/miozzi-2020-cell-basis"
CASES = ("printed", "tange_vinet", "speziale_variable_q_debye")


def direct_cell_pressure(volume_a3, temperature, parameters):
    """Independent whole-cell calculation: n=2, k_B and m^3 per cell.

    The quadrature helper normally takes R; substituting k_B makes its energy
    J/cell here. No per-Fe molar-volume conversion enters this expression.
    """
    v0 = parameters["V0_cell_a3"]
    gamma, theta = constant_q_gamma_theta(
        np.asarray(volume_a3) / v0, parameters["gamma0"], parameters["q"], 420
    )
    delta_u_cell = debye_energy(
        theta, temperature, 2, gas_constant=R / NA
    ) - debye_energy(theta, 300, 2, gas_constant=R / NA)
    return (
        bm3(volume_a3, v0, parameters["K0"], parameters["K0_prime"])
        + gamma * delta_u_cell / (np.asarray(volume_a3) * 1e-30) / 1e9
    )


def staged_fit(data, pressures, n):
    # J/bar/mol of Fe for n=1; J/bar/mol of whole hcp cells for n=2.
    scale = NA * 1e-25 / (2 if n == 1 else 1)
    volume = data["volume_a3"] * scale
    temperature = data["temperature_k"]
    cold = temperature <= 300
    bounds = {"V0": (20 * scale, 25 * scale), "K0": (20, 400), "K0_prime": (1, 10)}
    rt = fit_rt_eos(
        BM3,
        volume[cold],
        pressures[cold],
        initial={"V0": 22.81 * scale, "K0": 129, "K0_prime": 6.24},
        bounds=bounds,
        max_nfev=3000,
    )
    thermal_bounds = {"gamma0": (0.1, 5), "q": (-5, 10)}
    thermal = fit_thermal_eos(
        MieGruneisenDebye,
        rt.model,
        volume,
        temperature,
        pressures,
        initial={"gamma0": 1.11, "q": 0.3},
        fixed={"theta0": 420, "Tr": 300, "n": n},
        bounds=thermal_bounds,
        max_nfev=3000,
    )
    joint = fit_joint_eos(
        MieGruneisenDebye,
        BM3,
        volume,
        temperature,
        pressures,
        initial={
            **{f"rt_eos.{key}": rt.parameters[key] for key in bounds},
            **{key: thermal.parameters[key] for key in thermal_bounds},
        },
        fixed={"theta0": 420, "Tr": 300, "n": n},
        bounds={
            **{f"rt_eos.{key}": value for key, value in bounds.items()},
            **thermal_bounds,
        },
        max_nfev=3000,
    )
    stages = (rt, thermal, joint)
    if not all(stage.success for stage in stages):
        raise ValueError("A Peritheos fit stage failed to converge")
    names = ["rt_eos.V0", "rt_eos.K0", "rt_eos.K0_prime", "gamma0", "q"]
    physical_names = ["V0_cell_a3", "K0", "K0_prime", "gamma0", "q"]
    factors = np.array([1 / scale, 1, 1, 1, 1])
    parameters = {
        name: joint.parameters[key] * factor
        for name, key, factor in zip(physical_names, names, factors)
    }
    errors = {
        name: joint.standard_errors[key] * factor
        for name, key, factor in zip(physical_names, names, factors)
    }
    predicted = np.asarray(joint.model.pressure(volume, temperature))
    independent = direct_cell_pressure(data["volume_a3"], temperature, parameters)
    independent_delta = float(np.max(np.abs(independent - predicted)))
    if independent_delta > 1e-8:
        raise ValueError(
            "Native molar model disagrees with direct cell/k_B calculation"
        )
    # Check the *same* parameters under the other normalization separately
    # from independently optimized fits, whose stopping points may differ.
    alternate = MieGruneisenDebye(
        rt_eos=BM3(
            V0=joint.parameters["rt_eos.V0"] * (2 if n == 1 else 0.5),
            K0=parameters["K0"],
            K0_prime=parameters["K0_prime"],
        ),
        theta0=420,
        Tr=300,
        n=3 - n,
        gamma0=parameters["gamma0"],
        q=parameters["q"],
    )
    alternate_delta = float(
        np.max(
            np.abs(
                alternate.pressure(volume * (2 if n == 1 else 0.5), temperature)
                - predicted
            )
        )
    )
    if alternate_delta > 1e-10:
        raise ValueError("Consistent atom-count/volume scaling changed pressure")
    order = [joint.free_parameters.index(name) for name in names]
    covariance = joint.covariance[np.ix_(order, order)] * factors[:, None] * factors
    return {
        "normalization": {
            "n": n,
            "molar_entity": "Fe atom" if n == 1 else "two-atom hcp cell",
            "volume_unit": "J/bar/mol",
            "cell_a3_to_model_volume_scale": scale,
            "initial_V0_model": 22.81 * scale,
            "initial_V0_cm3_mol": 22.81 * scale * 10,
            "volume_error_scale": scale,
        },
        "parameters": parameters,
        "conditional_standard_errors": errors,
        "covariance_parameter_order": physical_names,
        "covariance_in_cell_volume_basis": covariance.tolist(),
        "rmse_gpa": float(np.sqrt(np.mean((predicted - pressures) ** 2))),
        "predicted_pressure_gpa": predicted.tolist(),
        "direct_cell_kB_max_difference_gpa": independent_delta,
        "same_parameters_other_basis_max_difference_gpa": alternate_delta,
        "stages": [stage.to_dict() for stage in stages],
    }


def run():
    from peritheos import get_material_document

    _, data, hashes = read_data()
    cases = {}
    for case in CASES:
        pressures = (
            target_pressures(data, "vinet")
            if case == "tange_vinet"
            else pressure_column(data, case)[0]
        )
        fits = {f"n{n}": staged_fit(data, pressures, n) for n in (1, 2)}
        delta = float(
            np.max(
                np.abs(
                    np.array(fits["n1"]["predicted_pressure_gpa"])
                    - np.array(fits["n2"]["predicted_pressure_gpa"])
                )
            )
        )
        # A separate numerical convergence check, not the exact scaling identity.
        if delta > 1e-4:
            raise ValueError(f"Independently optimized bases diverged: {case}: {delta}")
        cases[case] = {
            "observed_pressure_gpa": pressures.tolist(),
            "fits": fits,
            "independent_fit_max_pressure_difference_gpa": delta,
            "independent_fit_parameter_differences_n2_minus_n1": {
                name: fits["n2"]["parameters"][name] - value
                for name, value in fits["n1"]["parameters"].items()
            },
        }
    return {
        "format": "peritheos.miozzi-2020-cell-basis",
        "format_version": 1,
        "date": "2026-10-03",
        "source_csv_sha256": hashes,
        "script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "dependency_script_sha256": {
            path: hashlib.sha256((ROOT / path).read_bytes()).hexdigest()
            for path in (
                "scripts/reconstruct_miozzi_2020_iron.py",
                "scripts/refit_miozzi_2020_tange.py",
            )
        },
        "tange_calibration_record_sha256": hashlib.sha256(
            json.dumps(
                next(
                    record
                    for record in get_material_document("mgo")["eos_records"]
                    if record["identifier"] == "mgo_b1_tange_2009_vinet"
                ),
                sort_keys=True,
                separators=(",", ":"),
            ).encode()
        ).hexdigest(),
        "environment": {"python": platform.python_version(), "numpy": np.__version__},
        "protocol": "36 RT rows (298/300 K); all 131 rows with RT coefficients fixed; all 131 rows with V0,K0,K0_prime,gamma0,q free. theta0=420 K and Tr=300 K fixed throughout thermal stages.",
        "objective": "Equal pressure-residual weights; central coordinates held fixed; linear least squares. Bounds are identical after conversion to physical cell volumes.",
        "uncertainty": "Conditional local standard errors and covariance scaled by RSS/degrees of freedom. Calibration, fixed-theta, predictor and inter-row covariance uncertainties excluded. Stage-two errors condition on the cold coefficients.",
        "qualification": "Independent Peritheos normalization experiment. Whole-cell n=2 also doubles the molar volumes; this differs from the user-reported author setting n=2 with V0 about 6.87 cm3/mol. It does not resolve original author inputs or physically validate a pressure scale. Library records and defaults are unchanged.",
        "observations": {
            "volume_cell_a3": data["volume_a3"].tolist(),
            "temperature_k": data["temperature_k"].tolist(),
        },
        "cases": cases,
    }


def comparison(report):
    lines = [
        "# Whole-hcp-cell MGD fit with Peritheos",
        "",
        report["protocol"],
        "",
        report["objective"],
        "",
        "The n=1 fit uses volume per mole of Fe; the n=2 fit uses volume per mole of whole two-atom cells. Both start from 22.81 Å³/cell, respectively 6.86825 and 13.73650 cm³/mol (0.686825 and 1.373650 J/bar/mol). The fitted volume below is converted back to the same physical cell basis.",
        "",
        "| Pressure input | n | V0 (Å³/cell) | K0 (GPa) | K0′ | gamma0 | q | RMS (GPa) |",
        "|---|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for case, result in report["cases"].items():
        for key, fit in result["fits"].items():
            values = [
                f"{fit['parameters'][name]:.6f} ± {fit['conditional_standard_errors'][name]:.6f}"
                for name in fit["covariance_parameter_order"]
            ]
            lines.append(
                f"| {case} | {fit['normalization']['n']} | "
                + " | ".join(values)
                + f" | {fit['rmse_gpa']:.8f} |"
            )
    lines += [
        "",
        "## Verification",
        "",
        "| Pressure input | Independently refitted n=1 versus n=2: maximum ΔP (GPa) | n=2 versus direct cell/kB calculation: maximum ΔP (GPa) |",
        "|---|---:|---:|",
    ]
    for case, result in report["cases"].items():
        lines.append(
            f"| {case} | {result['independent_fit_max_pressure_difference_gpa']:.3g} | {result['fits']['n2']['direct_cell_kB_max_difference_gpa']:.3g} |"
        )
    lines += [
        "",
        "All 18 stages converged. The same coefficients under matching atom-count/volume scaling also give identical pressures within numerical precision, separately from optimizer stopping differences.",
        "",
        report["uncertainty"],
        "",
        report["qualification"],
        "",
        "Complete stage results, residuals, covariance, normalization and input fingerprints are retained in [report.json](report.json).",
        "",
        "Reproduce from the repository root:",
        "",
        "```sh",
        ".venv/bin/python -m scripts.refit_miozzi_2020_cell_basis",
        "```",
        "",
    ]
    return "\n".join(lines)


def main():
    report = run()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    (OUTPUT / "report.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    (OUTPUT / "comparison.md").write_text(comparison(report))
    for case, result in report["cases"].items():
        fit = result["fits"]["n2"]
        print(
            case,
            json.dumps(
                {
                    "parameters": fit["parameters"],
                    "rmse_gpa": fit["rmse_gpa"],
                    "basis_pressure_difference_gpa": result[
                        "independent_fit_max_pressure_difference_gpa"
                    ],
                }
            ),
        )


if __name__ == "__main__":
    main()
