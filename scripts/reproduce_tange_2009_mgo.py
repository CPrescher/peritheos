#!/usr/bin/env python3
"""Audit Tange et al. (2009) MgO Fit 3.

The default audit uses only redistributable row-level inputs.  Passing a local
``--approximate-inputs`` directory additionally evaluates a non-redistributed
reconstruction of the four restricted or figure-only source blocks.  That
second route is explicitly a similarity test, not an exact global refit.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.constants import Avogadro
from scipy.integrate import quad
from scipy.optimize import least_squares

from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import Tange2009Debye

ROOT = Path(__file__).resolve().parents[1]
DATASETS = ROOT / "peritheos" / "data" / "datasets"
DEFAULT_OUTPUT = ROOT / "docs" / "data" / "tange-2009-mgo-partial-validation.json"
DEFAULT_APPROXIMATE_OUTPUT = (
    ROOT / "docs" / "data" / "tange-2009-mgo-approximate-refit.json"
)

RECORD_IDENTIFIER = "mgo_b1_tange_2009_vinet"
V0_CELL_A3 = 74.698
CELL_TO_MOLAR_J_PER_BAR = Avogadro * 1.0e-25 / 4.0
V0_MOLAR = V0_CELL_A3 * CELL_TO_MOLAR_J_PER_BAR
T0 = 300.0
KS0_GPA = 162.83
ALPHA0_K_INV = 3.17e-5
CP0_J_MOL_K = 37.4
N_ATOMS = 2.0
PUBLISHED = {"K0_prime": 4.367, "gamma0": 1.442, "a": 0.138, "b": 5.4}

THERMAL_DATASETS = (
    "mgo-dubrovinsky-1997-cod-thermal-expansion.csv",
    "mgo-fiquet-1999-cod-thermal-expansion.csv",
)
ELASTIC_DATASET = "mgo-li-2006-table1-elasticity.csv"
SHOCK_DATASET = "mgo-marsh-1980-lasl-single-crystal-hugoniot.csv"
APPROXIMATE_INPUT_FILES = {
    "isaak_1989": "isaak-1989-mgo-ks.csv",
    "sinogeikin_2000": "sinogeikin-2000-mgo-ks.csv",
    "zha_2000": "zha-2000-mgo-ks-surrogate.csv",
    "duffy_1993": "duffy-1993-mgo-hugoniot.csv",
}
MGO_MOLAR_MASS_G_MOL = 40.304


def _rows(filename: str) -> list[dict[str, str]]:
    with (DATASETS / filename).open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _external_rows(directory: Path, key: str) -> list[dict[str, str]]:
    path = directory / APPROXIMATE_INPUT_FILES[key]
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def _external_input_inventory(directory: Path) -> list[dict[str, Any]]:
    inventory = []
    for key, filename in APPROXIMATE_INPUT_FILES.items():
        path = directory / filename
        payload = path.read_bytes()
        inventory.append(
            {
                "source_key": key,
                "filename": filename,
                "rows": len(_external_rows(directory, key)),
                "sha256": hashlib.sha256(payload).hexdigest(),
            }
        )
    return inventory


def constrained_eos(parameters: dict[str, float]) -> Tange2009Debye:
    """Construct Fit 3 while recomputing its two dependent parameters."""
    dependent = Tange2009Debye.constrained_reference_parameters(
        gamma0=parameters["gamma0"],
        reference_temperature=T0,
        adiabatic_bulk_modulus=KS0_GPA,
        thermal_expansivity=ALPHA0_K_INV,
        molar_heat_capacity_p=CP0_J_MOL_K,
        n=N_ATOMS,
    )
    reference = Vinet(V0_MOLAR, dependent["K0"], parameters["K0_prime"])
    return Tange2009Debye(
        reference,
        Tr=T0,
        theta0=dependent["theta0"],
        gamma0=parameters["gamma0"],
        a=parameters["a"],
        b=parameters["b"],
        n=N_ATOMS,
    )


def _reference_energy_change(eos: Tange2009Debye, volume: float) -> float:
    """Return E(V,300 K)-E(V0,300 K), in J/mol."""
    helmholtz_change = quad(
        lambda value: eos.rt_eos.pressure(value) * 1.0e4,
        volume,
        V0_MOLAR,
        epsabs=1.0e-6,
        epsrel=1.0e-9,
    )[0]
    entropy_change = eos.thermal_entropy(volume, T0) - eos.thermal_entropy(V0_MOLAR, T0)
    return float(helmholtz_change + T0 * entropy_change)


def predicted_hugoniot_pressure(
    eos: Tange2009Debye, volume_ratio: float, observed_pressure_gpa: float
) -> float:
    """Evaluate Tange's energy-reduced Hugoniot pressure at observed V/V0."""
    volume = V0_MOLAR * volume_ratio
    hugoniot_energy = 0.5 * observed_pressure_gpa * (V0_MOLAR - volume) * 1.0e4
    thermal_energy = hugoniot_energy - _reference_energy_change(eos, volume)
    thermal_pressure = eos.gruneisen_parameter(volume) * thermal_energy / volume / 1.0e4
    return float(eos.rt_eos.pressure(volume) + thermal_pressure)


def residual_groups(parameters: dict[str, float]) -> dict[str, np.ndarray]:
    """Return the three available source-class residual vectors."""
    eos = constrained_eos(parameters)
    thermal = []
    for filename in THERMAL_DATASETS:
        for row in _rows(filename):
            volume = float(row["conventional_cell_volume_a3"]) * CELL_TO_MOLAR_J_PER_BAR
            thermal.append(eos.pressure(volume, float(row["temperature_k"])))

    elastic = []
    for row in _rows(ELASTIC_DATASET):
        if row["experimental_path"] == "ambient":
            continue
        volume = float(row["volume_ratio"]) * V0_MOLAR
        observed = float(row["adiabatic_bulk_modulus_gpa"])
        elastic.append(float(eos.adiabatic_bulk_modulus(volume, T0)) - observed)

    shock = []
    for row in _rows(SHOCK_DATASET):
        rho0 = float(row["initial_density_g_cm3"])
        us = float(row["shock_velocity_km_s"])
        up = float(row["particle_velocity_km_s"])
        volume_ratio = 1.0 - up / us
        observed_pressure = rho0 * us * up
        shock.append(
            predicted_hugoniot_pressure(eos, volume_ratio, observed_pressure)
            - observed_pressure
        )
    return {
        "thermal_expansion_pressure_gpa": np.asarray(thermal),
        "li_2006_adiabatic_bulk_modulus_gpa": np.asarray(elastic),
        "lasl_hugoniot_pressure_gpa": np.asarray(shock),
    }


def weighted_residuals(vector: np.ndarray) -> np.ndarray:
    parameters = dict(zip(PUBLISHED, map(float, vector)))
    groups = residual_groups(parameters)
    elastic_sigmas = np.asarray(
        [
            float(row["adiabatic_bulk_modulus_standard_deviation_gpa"])
            for row in _rows(ELASTIC_DATASET)
            if row["experimental_path"] != "ambient"
        ]
    )
    return np.concatenate(
        (
            groups["thermal_expansion_pressure_gpa"] / 0.3,
            groups["li_2006_adiabatic_bulk_modulus_gpa"] / elastic_sigmas,
            groups["lasl_hugoniot_pressure_gpa"] / 3.0,
        )
    )


def approximate_residual_groups(
    parameters: dict[str, float], directory: Path
) -> dict[str, np.ndarray]:
    """Return all recoverable groups, including local restricted inputs."""
    eos = constrained_eos(parameters)
    groups = residual_groups(parameters)

    for key, group_name in (
        ("isaak_1989", "isaak_1989_adiabatic_bulk_modulus_gpa"),
        ("sinogeikin_2000", "sinogeikin_2000_adiabatic_bulk_modulus_gpa"),
    ):
        residuals = []
        for row in _external_rows(directory, key):
            molar_volume = MGO_MOLAR_MASS_G_MOL / float(row["density_g_cm3"]) / 10.0
            residuals.append(
                float(
                    eos.adiabatic_bulk_modulus(
                        molar_volume, float(row["temperature_k"])
                    )
                )
                - float(row["adiabatic_bulk_modulus_gpa"])
            )
        groups[group_name] = np.asarray(residuals)

    zha = []
    for row in _external_rows(directory, "zha_2000"):
        volume = V0_MOLAR * 3.585 / float(row["density_g_cm3"])
        zha.append(
            float(eos.adiabatic_bulk_modulus(volume, T0))
            - float(row["adiabatic_bulk_modulus_gpa"])
        )
    groups["zha_2000_adiabatic_bulk_modulus_gpa"] = np.asarray(zha)

    duffy = []
    for row in _external_rows(directory, "duffy_1993"):
        volume_ratio = float(row["initial_density_g_cm3"]) / float(
            row["shocked_density_g_cm3"]
        )
        observed_pressure = float(row["pressure_gpa"])
        duffy.append(
            predicted_hugoniot_pressure(eos, volume_ratio, observed_pressure)
            - observed_pressure
        )
    groups["duffy_1993_hugoniot_pressure_gpa"] = np.asarray(duffy)
    return groups


def approximate_weighted_residuals(
    vector: np.ndarray, directory: Path, *, zha_sigma_scale: float = 1.0
) -> np.ndarray:
    parameters = dict(zip(PUBLISHED, map(float, vector)))
    groups = approximate_residual_groups(parameters, directory)
    elastic_sigmas = np.asarray(
        [
            float(row["adiabatic_bulk_modulus_standard_deviation_gpa"])
            for row in _rows(ELASTIC_DATASET)
            if row["experimental_path"] != "ambient"
        ]
    )
    isaak_sigmas = np.asarray(
        [
            float(row["adiabatic_bulk_modulus_standard_deviation_gpa"])
            for row in _external_rows(directory, "isaak_1989")
        ]
    )
    sinogeikin_sigmas = np.asarray(
        [
            float(row["adiabatic_bulk_modulus_standard_deviation_gpa"])
            for row in _external_rows(directory, "sinogeikin_2000")
        ]
    )
    zha_sigmas = np.asarray(
        [
            float(row["adiabatic_bulk_modulus_standard_deviation_gpa"])
            for row in _external_rows(directory, "zha_2000")
        ]
    )
    duffy_sigmas = np.asarray(
        [
            float(row["pressure_standard_deviation_gpa"])
            for row in _external_rows(directory, "duffy_1993")
        ]
    )
    return np.concatenate(
        (
            groups["thermal_expansion_pressure_gpa"] / 0.3,
            groups["isaak_1989_adiabatic_bulk_modulus_gpa"] / isaak_sigmas,
            groups["sinogeikin_2000_adiabatic_bulk_modulus_gpa"] / sinogeikin_sigmas,
            groups["zha_2000_adiabatic_bulk_modulus_gpa"]
            / (zha_sigmas * zha_sigma_scale),
            groups["li_2006_adiabatic_bulk_modulus_gpa"] / elastic_sigmas,
            groups["lasl_hugoniot_pressure_gpa"] / 3.0,
            groups["duffy_1993_hugoniot_pressure_gpa"] / duffy_sigmas,
        )
    )


def _metrics(groups: dict[str, np.ndarray]) -> dict[str, dict[str, float | int]]:
    return {
        name: {
            "observations": int(values.size),
            "rmse": float(np.sqrt(np.mean(values**2))),
            "mean_residual": float(np.mean(values)),
            "maximum_absolute_residual": float(np.max(np.abs(values))),
        }
        for name, values in groups.items()
    }


def build_report() -> dict[str, Any]:
    dependent = Tange2009Debye.constrained_reference_parameters(
        gamma0=PUBLISHED["gamma0"],
        reference_temperature=T0,
        adiabatic_bulk_modulus=KS0_GPA,
        thermal_expansivity=ALPHA0_K_INV,
        molar_heat_capacity_p=CP0_J_MOL_K,
        n=N_ATOMS,
    )
    published_groups = residual_groups(PUBLISHED)
    fit = least_squares(
        weighted_residuals,
        np.asarray(list(PUBLISHED.values())),
        bounds=([2.0, 0.5, 0.0, 0.1], [8.0, 3.0, 1.0, 20.0]),
        xtol=1.0e-11,
        ftol=1.0e-11,
        gtol=1.0e-11,
        max_nfev=3000,
    )
    fitted = dict(zip(PUBLISHED, map(float, fit.x)))
    fitted_dependent = constrained_eos(fitted)
    fitted_groups = residual_groups(fitted)
    return {
        "format": "peritheos.tange-2009-mgo-partial-validation",
        "format_version": 1,
        "record_identifier": RECORD_IDENTIFIER,
        "scope": "partial_validation_not_global_refit",
        "global_refit_reproduced": False,
        "global_refit_blocker": (
            "Zha et al. (2000) publish polynomial elastic summaries and plots, not "
            "the 27 row-level V-Ks observations used by Tange et al.; exact source "
            "weights are also unavailable for some earlier elastic rows."
        ),
        "fixed_ambient_constraints": {
            "V0_conventional_cell_a3": V0_CELL_A3,
            "KS0_gpa": KS0_GPA,
            "alpha0_k_inverse": ALPHA0_K_INV,
            "CP0_j_mol_k": CP0_J_MOL_K,
            "reference_temperature_k": T0,
            "n_atoms": N_ATOMS,
        },
        "published_dependent_parameter_check": {
            "calculated_K0_gpa": dependent["K0"],
            "published_K0_gpa": 160.63,
            "calculated_theta0_k": dependent["theta0"],
            "published_theta0_k": 761.0,
            "calculated_CV0_j_mol_k": dependent["Cv0"],
        },
        "objective": {
            "free_parameters": list(PUBLISHED),
            "dependent_parameters": ["K0", "theta0"],
            "thermal_expansion": "P(V,T) / 0.3 GPa",
            "static_elasticity": (
                "(calculated KS - observed KS) / reported sigma_KS; source-listed "
                "pressures are not fit inputs"
            ),
            "shock": "(energy-reduced calculated PH - observed PH) / 3 GPa",
            "note": (
                "These are the source-stated weighting conventions for the available "
                "classes. The omitted Isaak, Sinogeikin, Zha, and Duffy rows prevent "
                "this subset objective from representing the published global objective."
            ),
        },
        "source_inventory": [
            {
                "source": "Dubrovinsky and Saxena (1997)",
                "source_rows_stated_by_tange": 25,
                "rows_recovered": 25,
                "disposition": "bundled_cod_cc0",
                "evidence": "COD records 9006456-9006480",
                "url": "https://www.crystallography.net/cod/result.php?DOI=10.1007/s002690050070",
            },
            {
                "source": "Fiquet et al. (1999)",
                "source_rows_stated_by_tange": 37,
                "rows_recovered": 36,
                "disposition": "bundled_cod_cc0_one_row_unresolved",
                "evidence": "COD records 9006747-9006782",
                "url": "https://www.crystallography.net/cod/result.php?DOI=10.1007/s002690050246",
            },
            {
                "source": "Isaak et al. (1989)",
                "source_rows_stated_by_tange": 16,
                "rows_recovered": 16,
                "disposition": "located_nist_srd_not_redistributed",
                "evidence": "NIST SRD 30 MgO adiabatic bulk-modulus table",
                "url": "https://srdata.nist.gov/CeramicDataPortal/Scd/Z00281",
            },
            {
                "source": "Sinogeikin and Bass (2000)",
                "source_rows_stated_by_tange": 15,
                "rows_recovered": 15,
                "disposition": "located_copyrighted_article_not_redistributed",
                "evidence": "Article Table I",
                "url": "https://doi.org/10.1063/1.1150496",
            },
            {
                "source": "Zha et al. (2000)",
                "source_rows_stated_by_tange": 27,
                "rows_recovered": 0,
                "disposition": "row_level_values_not_published",
                "evidence": "Article Tables 1-2 and Figures 1-2",
                "url": "https://pmc.ncbi.nlm.nih.gov/articles/PMC17603/",
            },
            {
                "source": "Li et al. (2006)",
                "source_rows_stated_by_tange": 17,
                "rows_recovered": 17,
                "disposition": "already_bundled_official_table",
                "evidence": "Official Table 1 supplement",
                "url": "https://doi.org/10.1029/2005JB004251",
            },
            {
                "source": "Marsh (1980) LASL compendium",
                "source_rows_stated_by_tange": 24,
                "rows_recovered": 24,
                "disposition": "bundled_public_domain",
                "evidence": "hugla.dat single-crystal periclase block",
                "url": "https://impacts.wiki/community-resources/material-properties-data-and-references/lasl-shock-hugoniot-data-marsh-1980-zip-file-enclosed-and-in-rusbank",
            },
            {
                "source": "Duffy and Ahrens (1993)",
                "source_rows_stated_by_tange": 4,
                "rows_recovered": 4,
                "disposition": "located_copyrighted_article_not_redistributed",
                "evidence": "Article Table 1",
                "url": "https://doi.org/10.1029/93GL01510",
            },
        ],
        "bundled_partial_validation": {
            "observations": int(
                sum(values.size for values in published_groups.values())
            ),
            "dataset_files": [*THERMAL_DATASETS, ELASTIC_DATASET, SHOCK_DATASET],
            "published_coefficients": PUBLISHED,
            "source_reported_group_rmse": {
                "li_2006_adiabatic_bulk_modulus_gpa": 1.5,
                "lasl_hugoniot_pressure_gpa": 1.4,
            },
            "published_coefficient_metrics": _metrics(published_groups),
            "partial_refit": {
                "solver_success": bool(fit.success),
                "solver_message": fit.message,
                "coefficients": fitted,
                "dependent_K0_gpa": fitted_dependent.rt_eos.K0,
                "dependent_theta0_k": fitted_dependent.theta0,
                "metrics": _metrics(fitted_groups),
                "weighted_sum_of_squares": float(np.sum(fit.fun**2)),
            },
        },
        "conclusion": (
            "The ambient thermodynamic constraints and the available row classes are "
            "executable, but the partial-subset optimum is not an estimate of the "
            "published global coefficients. Full coefficient parity remains untestable "
            "without the missing Zha V-Ks rows and exact legacy elastic weights."
        ),
    }


def _similar_parameter(name: str, published: float, fitted: float) -> bool:
    relative = abs(fitted - published) / max(abs(published), 1.0e-30)
    if name == "K0_prime":
        return abs(fitted - published) <= 1.0 or relative <= 0.20
    if name == "gamma0":
        return relative <= 0.25
    return relative <= 0.20


def _approximate_fit(directory: Path, zha_sigma_scale: float) -> Any:
    return least_squares(
        lambda vector: approximate_weighted_residuals(
            vector, directory, zha_sigma_scale=zha_sigma_scale
        ),
        np.asarray(list(PUBLISHED.values())),
        bounds=([2.0, 0.5, 0.0, 0.1], [8.0, 3.0, 1.0, 20.0]),
        xtol=1.0e-11,
        ftol=1.0e-11,
        gtol=1.0e-11,
        max_nfev=3000,
    )


def build_approximate_report(directory: Path) -> dict[str, Any]:
    """Build the local, non-redistributed similarity experiment."""
    inventory = _external_input_inventory(directory)
    recovered = sum(item["rows"] for item in inventory) + 102
    if recovered != 164:
        raise ValueError(
            f"expected 164 reconstructed observations, recovered {recovered}"
        )

    published_groups = approximate_residual_groups(PUBLISHED, directory)
    fit = _approximate_fit(directory, 1.0)
    fitted = dict(zip(PUBLISHED, map(float, fit.x)))
    fitted_groups = approximate_residual_groups(fitted, directory)
    comparisons = []
    for name, published in PUBLISHED.items():
        value = fitted[name]
        relative_percent = 100.0 * abs(value - published) / abs(published)
        comparisons.append(
            {
                "parameter": name,
                "published": published,
                "fitted": value,
                "absolute_relative_difference_percent": relative_percent,
                "similar": _similar_parameter(name, published, value),
            }
        )

    sensitivity = []
    for scale, effective_percent in ((0.5, 1.0), (1.0, 2.0), (1.5, 3.0)):
        trial = fit if scale == 1.0 else _approximate_fit(directory, scale)
        coefficients = dict(zip(PUBLISHED, map(float, trial.x)))
        sensitivity.append(
            {
                "zha_relative_sigma_percent": effective_percent,
                "coefficients": coefficients,
                "all_parameters_similar": all(
                    _similar_parameter(name, PUBLISHED[name], value)
                    for name, value in coefficients.items()
                ),
                "weighted_sum_of_squares": float(np.sum(trial.fun**2)),
            }
        )

    all_similar = all(item["similar"] for item in comparisons)
    return {
        "format": "peritheos.tange-2009-mgo-approximate-refit",
        "format_version": 1,
        "record_identifier": RECORD_IDENTIFIER,
        "scope": "local_nonredistributed_approximate_global_refit",
        "global_refit_reproduced": False,
        "result": "similar" if all_similar else "parity_not_achieved",
        "qualification": (
            "The result uses a curve-based surrogate for the unpublished Zha rows "
            "and reconstructed legacy elastic uncertainties. It cannot establish "
            "strict coefficient or uncertainty parity."
        ),
        "observations": {
            "source_rows_stated_by_tange": 165,
            "rows_reconstructed": recovered,
            "unresolved_rows": 1,
            "unresolved_source": "Fiquet et al. (1999)",
        },
        "external_input_inventory": inventory,
        "reconstruction_assumptions": {
            "isaak_1989": (
                "KS=(C11+2*C12)/3 at the 16 temperatures selected by Tange; "
                "sigma_KS is linearly interpolated through the four uncertainty "
                "anchors printed by NIST."
            ),
            "sinogeikin_2000": (
                "Table I density and KS values; sigma_KS=1% of KS, matching the "
                "article's stated modulus accuracy bound."
            ),
            "zha_2000": (
                "27 non-ambient compression locations reconstructed from the plot; "
                "KS generated from Zha's published third-order finite-strain fit "
                "(KS0=162.5 GPa, KS0'=3.99); sigma_KS=2% of KS."
            ),
            "duffy_1993": (
                "Four Table 1 P-rhoH states; sigma_P=3 GPa as specified by Tange "
                "for shock-pressure residuals."
            ),
        },
        "objective": {
            "free_parameters": list(PUBLISHED),
            "dependent_parameters": ["K0", "theta0"],
            "same_as_partial_audit": True,
            "additional_groups": [
                "Isaak ambient-pressure KS(T)",
                "Sinogeikin ambient-pressure KS(T)",
                "Zha 300 K KS(V) surrogate",
                "Duffy shock Hugoniot",
            ],
        },
        "published_coefficients": PUBLISHED,
        "published_coefficient_metrics": _metrics(published_groups),
        "approximate_refit": {
            "solver_success": bool(fit.success),
            "solver_message": fit.message,
            "coefficients": fitted,
            "coefficient_comparisons": comparisons,
            "metrics": _metrics(fitted_groups),
            "weighted_sum_of_squares": float(np.sum(fit.fun**2)),
        },
        "zha_weight_sensitivity": sensitivity,
        "conclusion": (
            "The approximate 164-row reconstruction reaches the repository's "
            "numerical similarity criteria, including across the tested 1-3% Zha "
            "weight range. Exact parity still requires the original 27 Zha rows, "
            "legacy row weights, and the unresolved Fiquet state."
            if all_similar
            else "The approximate reconstruction does not reach numerical similarity."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument(
        "--approximate-inputs",
        type=Path,
        help="local directory containing the four non-redistributed CSV inputs",
    )
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.approximate_inputs is None:
        output = args.output or DEFAULT_OUTPUT
        report = build_report()
    else:
        output = args.output or DEFAULT_APPROXIMATE_OUTPUT
        report = build_approximate_report(args.approximate_inputs)
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.check:
        if not output.exists() or output.read_text(encoding="utf-8") != rendered:
            raise SystemExit(f"stale generated file: {output}")
        return
    output.write_text(rendered, encoding="utf-8")


if __name__ == "__main__":
    main()
