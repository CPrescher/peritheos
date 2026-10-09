#!/usr/bin/env python3
"""Diagnose CaSiO3 volume/atom-count normalization with Fu's fixed coefficients.

This does not refit or replace the published EOS. Non-stoichiometric atom counts
are thermal-amplitude diagnostics only. Run with the bundled Sun/Greaux rows.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import audit_fu_2023_casio3_refit as fu
import numpy as np
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
RECORD = "ca_perovskite_fu_2023_bm3_mgd_refit"
BOLTZMANN = 1.380649e-23


def predictions(volume, temperature, n=5.0, formula_units=1.0):
    """Evaluate the full audited model on an explicitly represented volume basis."""
    parameters = fu.PUBLISHED.copy()
    parameters[0] *= formula_units
    original = fu.ATOM_COUNT
    try:
        fu.ATOM_COUNT = n
        return fu._predictions(parameters, volume * formula_units, temperature)
    finally:
        fu.ATOM_COUNT = original


def atomic_pressure(volume, temperature):
    """Independent SI evaluation: energy per CaSiO3 formula divided by its m^3."""
    v0, k0, _, _, gamma0, q, _ = fu.PUBLISHED
    f = 0.5 * ((v0 / volume) ** (2.0 / 3.0) - 1.0)
    cold = 3.0 * k0 * f * (1.0 + 2.0 * f) ** 2.5
    gamma = gamma0 * (volume / v0) ** q
    theta = fu.DEBYE_TEMPERATURE_K * np.exp((gamma0 - gamma) / q)

    def energy(t):
        integral = quad(lambda x: x**3 / np.expm1(x), 0.0, theta / t)[0]
        return 9.0 * 5.0 * BOLTZMANN * t * (t / theta) ** 3 * integral

    return cold + gamma * (energy(temperature) - energy(300.0)) / (volume * 1e-30) / 1e9


def rmse(values):
    return float(np.sqrt(np.mean(np.asarray(values) ** 2)))


def amplitude_diagnostic(observed, cold, thermal):
    scale = float(np.dot(thermal, observed - cold) / np.dot(thermal, thermal))
    return {
        "outputs": len(observed),
        "effective_n": 5.0 * scale,
        "thermal_scale_relative_to_n5": scale,
        "rmse_gpa": rmse(cold + scale * thermal - observed),
        "interpretation": "Unweighted scalar-output amplitude diagnostic; other coefficients fixed. Not a physical atom count or a new EOS.",
    }


def audit():
    sun, greaux = fu.load_bundled_observations()
    x = fu._arrays(sun, greaux)
    sun5 = predictions(x["sun_v"], x["sun_t"])
    sun0 = predictions(x["sun_v"], np.full_like(x["sun_t"], 300.0))
    g5 = predictions(x["g_v"], x["g_t"])
    g0 = predictions(x["g_v"], np.full_like(x["g_t"], 300.0))
    observed = [x["sun_p"], x["g_p"], x["g_k"], x["g_mu"]]
    cold = [sun0[0], g0[0], g0[1], g0[2]]
    thermal = [sun5[0] - sun0[0], *(g5[i] - g0[i] for i in range(3))]
    keys = [
        "sun_pressure",
        "greaux_pressure",
        "greaux_bulk_modulus",
        "greaux_shear_modulus",
    ]
    amplitudes = {
        k: amplitude_diagnostic(o, c, d)
        for k, o, c, d in zip(keys, observed, cold, thermal)
    }
    amplitudes["joint_equal_scalar_gpa_weights"] = amplitude_diagnostic(
        np.concatenate(observed), np.concatenate(cold), np.concatenate(thermal)
    )

    record = Material.from_eosmat(
        get_material_document("ca_perovskite"), record_identifiers=[RECORD]
    ).get_eos_record(RECORD)
    all_v = np.concatenate([x["sun_v"], x["g_v"]])
    all_t = np.concatenate([x["sun_t"], x["g_t"]])
    molar = predictions(all_v, all_t)[0]
    atomic = np.array([atomic_pressure(v, t) for v, t in zip(all_v, all_t)])
    loaded = np.array(
        [record.pressure(float(v), float(t)) for v, t in zip(all_v, all_t)]
    )
    checks = {
        "public_a3_to_model_j_per_bar_per_mol": float(record.volume_scale),
        "expected_public_a3_to_model_scale": fu.AVOGADRO * 1e-25,
        "max_atomic_vs_molar_pressure_difference_gpa": float(
            np.max(np.abs(atomic - molar))
        ),
        "max_loaded_vs_independent_atomic_pressure_difference_gpa": float(
            np.max(np.abs(loaded - atomic))
        ),
        "consistent_basis_changes": [],
    }
    full = predictions(all_v, all_t)
    for z in (2.0, 4.0):
        scaled = predictions(all_v, all_t, n=5.0 * z, formula_units=z)
        checks["consistent_basis_changes"].append(
            {
                "formula_units": z,
                "n": 5.0 * z,
                "V0_a3": fu.PUBLISHED[0] * z,
                "max_output_difference_gpa": max(
                    float(np.max(np.abs(a - b))) for a, b in zip(full, scaled)
                ),
            }
        )
    assert abs(record.volume_scale - fu.AVOGADRO * 1e-25) < 1e-14
    assert checks["max_atomic_vs_molar_pressure_difference_gpa"] < 1e-10
    assert checks["max_loaded_vs_independent_atomic_pressure_difference_gpa"] < 1e-8
    assert all(
        c["max_output_difference_gpa"] < 1e-10
        for c in checks["consistent_basis_changes"]
    )

    cases = []
    mask2200 = x["sun_t"] == 2200.0
    for n in (1.0, 2.0, 5.0, 10.0, 20.0, amplitudes["sun_pressure"]["effective_n"]):
        sp = predictions(x["sun_v"], x["sun_t"], n=n)[0]
        gp, gk, gm = predictions(x["g_v"], x["g_t"], n=n)
        cases.append(
            {
                "n": n,
                "volume_formula_units": 1.0,
                "physical_stoichiometry": n == 5.0,
                "sun_pressure_rmse_gpa": rmse(sp - x["sun_p"]),
                "sun_2200k_pressure_rmse_gpa": rmse((sp - x["sun_p"])[mask2200]),
                "greaux_pressure_rmse_gpa": rmse(gp - x["g_p"]),
                "greaux_bulk_modulus_rmse_gpa": rmse(gk - x["g_k"]),
                "greaux_shear_modulus_rmse_gpa": rmse(gm - x["g_mu"]),
            }
        )

    checkpoints = []
    for n in (5.0, 10.0):
        p, k, mu = predictions(np.array([38.5]), np.array([2200.0]), n=n)
        p0 = predictions(np.array([38.5]), np.array([300.0]), n=n)[0]
        checkpoints.append(
            {
                "volume_a3": 38.5,
                "temperature_k": 2200.0,
                "observed_pressure_gpa": 72.8,
                "n": n,
                "pressure_gpa": float(p[0]),
                "thermal_increment_gpa": float(p[0] - p0[0]),
            }
        )

    # Separate pressure amplitude from the volume dependence of gamma. This is
    # a conditional diagnostic, not Fu's undisclosed joint regression.
    def pressure_residual(thermal_pair):
        parameters = fu.PUBLISHED.copy()
        parameters[4:6] = thermal_pair
        return fu._predictions(parameters, x["sun_v"], x["sun_t"])[0] - x["sun_p"]

    fit = least_squares(
        pressure_residual,
        [1.42, 2.65],
        bounds=([0.1, 0.05], [5.0, 6.0]),
    )
    assert fit.success
    fitted = fu.PUBLISHED.copy()
    fitted[4:6] = fit.x
    gp, gk, gm = fu._predictions(fitted, x["g_v"], x["g_t"])
    thermal_pair_diagnostic = {
        "scope": "Unweighted pressure-only fit to 140 Sun rows, n=5 and all reference/shear/theta parameters fixed; not a replacement EOS or original-fit reconstruction",
        "gamma0": float(fit.x[0]),
        "q": float(fit.x[1]),
        "sun_pressure_rmse_gpa": rmse(fit.fun),
        "greaux_pressure_rmse_gpa": rmse(gp - x["g_p"]),
        "greaux_bulk_modulus_rmse_gpa": rmse(gk - x["g_k"]),
        "greaux_shear_modulus_rmse_gpa": rmse(gm - x["g_mu"]),
        "gamma_at_38_5_a3": {
            "published": float(
                fu.PUBLISHED[4] * (38.5 / fu.PUBLISHED[0]) ** fu.PUBLISHED[5]
            ),
            "pressure_only_diagnostic": float(
                fitted[4] * (38.5 / fitted[0]) ** fitted[5]
            ),
        },
    }

    return {
        "record": RECORD,
        "status": "normalization_consistent_cause_of_published_parameter_mismatch_unresolved",
        "sources": {
            "published_coefficients": "Fu et al. (2023), supplement Table S3",
            "equations": "Fu et al. (2023), main text Equations 23–28; standard dimensionally correct CV in Equation 27",
            "observations": "140 Sun (2016) Table 1 rows at 1200–2200 K and 34 cubic Gréaux (2019) Figure 3b rows; exact Fu selection undisclosed",
            "greaux_pressure": "Measured PNaCl; derived PFS is not treated as an independent measurement",
        },
        "normalization": {
            "n": 5.0,
            "formula_units": 1.0,
            "V0_a3": 45.4,
            "V0_cm3_per_mol": 45.4 * fu.AVOGADRO * 1e-24,
            "heat_capacity_limit_j_per_mol_k": 3.0 * 5.0 * fu.R,
            "conversion": "1 A^3 per CaSiO3 formula = 0.602214076 cm^3 per mol CaSiO3 = 0.0602214076 J/bar per mol CaSiO3",
        },
        "normalization_checks": checks,
        "fixed_coefficient_n_diagnostics": cases,
        "optimal_amplitudes_by_observable": amplitudes,
        "conditional_thermal_pair_diagnostic": thermal_pair_diagnostic,
        "2200k_checkpoint": checkpoints,
        "conclusion": "A consistent change of volume/formula basis leaves all outputs unchanged. Doubling n alone improves Sun pressure but worsens Greaux elastic moduli. No single atom-count/volume-density scaling reconciles the joint observations with the fixed published coefficients. The checks do not determine which normalization or model the authors actually used.",
    }, x


def plot_diagnostics(x, destination):
    import matplotlib.pyplot as plt

    fig, axes = plt.subplots(2, 2, figsize=(11.0, 7.2), constrained_layout=True)
    panels = [
        ("Sun pressure", "sun_v", "sun_t", "sun_p", 0),
        ("Gréaux pressure", "g_v", "g_t", "g_p", 0),
        ("Gréaux bulk modulus", "g_v", "g_t", "g_k", 1),
        ("Gréaux shear modulus", "g_v", "g_t", "g_mu", 2),
    ]
    for ax, (title, vk, tk, ok, index) in zip(axes.flat, panels):
        for n, color, marker in [(5.0, "#b13c2e", "o"), (10.0, "#2468a2", "x")]:
            residual = predictions(x[vk], x[tk], n=n)[index] - x[ok]
            ax.scatter(
                x[tk],
                residual,
                color=color,
                marker=marker,
                s=22,
                alpha=0.75,
                label=f"n = {n:g}; RMSE {rmse(residual):.2f} GPa",
            )
        ax.axhline(0, color="#666666", linewidth=0.8)
        ax.set(
            title=title, xlabel="Temperature (K)", ylabel="Model - observation (GPa)"
        )
        ax.legend(fontsize=9)
        ax.grid(alpha=0.2)
    fig.suptitle(
        "Fu (2023) CaSiO3: changing n alone does not reconcile the joint data",
        fontsize=14,
    )
    fig.savefig(destination, dpi=180)
    plt.close(fig)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--output",
        type=Path,
        default=ROOT / "docs/data/fu-2023-normalization-audit.json",
    )
    parser.add_argument("--plot", type=Path)
    args = parser.parse_args()
    result, arrays = audit()
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(result, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    if args.plot:
        args.plot.parent.mkdir(parents=True, exist_ok=True)
        plot_diagnostics(arrays, args.plot)
    print(
        json.dumps(
            {
                "output": str(args.output),
                "plot": str(args.plot) if args.plot else None,
                "status": result["status"],
                "optimal_effective_n": {
                    k: v["effective_n"]
                    for k, v in result["optimal_amplitudes_by_observable"].items()
                },
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
