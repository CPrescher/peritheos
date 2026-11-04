"""Scoped Fe3C source checks; repository maintenance, not a replacement EOS fit."""

from __future__ import annotations

import copy
import csv
import json
from pathlib import Path
from typing import Any

import numpy as np
from scipy.constants import Avogadro

from peritheos import get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_joint_eos, fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
RESEARCH = ROOT / "peritheos/data/research/cementite"
REPORT = ROOT / "docs/data/cementite-reproduction.json"
CELL_TO_MOLAR = Avogadro * 1e-25 / 4
SCOTT_ID = "cementite_scott_2001_bm3_1"
SCOTT_BASIS = (
    "Qualified similarity of a digitized, equal-weight reconstruction. Both free "
    "central coefficients lie outside their printed source error widths; the "
    "broad combined-error screen does not establish coefficient parity."
)
SCOTT_QUALIFICATION = (
    "Scott et al. (2001), page 1876 and Figure 3, report weighted least squares "
    "on observations from two runs. Only 21 digitized filled-symbol observations "
    "are available here; the open temperature-quenched point is excluded. "
    "Incomplete digitized coordinate error bars cause the generic row reader to "
    "omit both sigma arrays, so this reconstruction uses equal pressure weights. "
    "Its retained errors use residual-scaled covariance because no observation "
    "sigmas are supplied, even though the generic ledger records "
    "absolute_sigma=True. They are not recovered source weights or a source "
    "confidence interval. The original row array, run identities, "
    "weights, covariance and precise objective remain unavailable. The source "
    "coefficients and digitization are preserved; this finding does not show "
    "that the published EOS is incorrect."
)


def qualify_scott(outcome: dict[str, Any]) -> dict[str, Any]:
    """Prevent broad uncertainty overlap from overstating this partial replay.

    The dedicated central-value criterion is one printed source error width for
    every free coefficient. This is a numerical screen, not a confidence test:
    Scott does not specify the confidence convention. Missing widths cannot pass.
    Preserve the generic comparison flags and numerical fit without alteration.
    """
    result = copy.deepcopy(outcome)
    comparisons = result.get("parameters", [])
    agreement = bool(comparisons) and all(
        p.get("published_error") is not None
        and p["published_error"] > 0
        and abs(p["refit"] - p["published"]) <= p["published_error"]
        for p in comparisons
    )
    result["source_parameter_widths_agreement"] = agreement
    if result["status"] == "parity" and not agreement:
        result["status"] = "similar"
    result["parity_basis"] = (
        "Every free central coefficient lies within its printed source error width; "
        "this conditional numerical comparison does not establish the source's "
        "confidence convention or weighted fitting workflow."
        if agreement
        else SCOTT_BASIS
    )
    result["qualification"] = SCOTT_QUALIFICATION
    result["weighting"] = (
        "Equal pressure weights; incomplete plotted errors are not used as fit "
        "sigmas. Refit covariance is scaled by residual variance."
    )
    result["covariance_convention"] = (
        "Residual-scaled inverse normal matrix; no measured sigma arrays supplied. "
        "The fitting API scales covariance in that case regardless of absolute_sigma."
    )
    result["original_author_workflow_parity"] = False
    result["reproduction_report"] = "docs/literature-reproductions/cementite.md"
    return result


def reproduce_scott() -> dict[str, Any]:
    from scripts.validate_primary_eos_refits import _fit_record, _series

    document = get_material_document("cementite")
    record = document["eos_records"][0]
    dataset = document["datasets"][0]
    series = _series(document, record, dataset)
    outcome = _fit_record(document, record, dataset)
    result = fit_rt_eos(
        BM3,
        series.volume,
        series.pressure,
        initial={"K0": 175.4, "K0_prime": 5.1},
        fixed={"V0": 155.26},
        bounds={"K0": (8.77, 877.0), "K0_prime": (0, 20)},
        absolute_sigma=True,
        max_nfev=5000,
    )
    source = BM3(V0=155.26, K0=175.4, K0_prime=5.1)
    pressures = np.linspace(0, 73.2, 1001)
    original = source.volume(pressures)
    reconstructed = result.model.volume(pressures)
    covariance = np.asarray(result.covariance)
    return {
        **outcome,
        "pressure_sigmas_used": series.pressure_sigma is not None,
        "volume_sigmas_used": series.volume_sigma is not None,
        "covariance_order": list(result.free_parameters),
        "covariance": covariance.tolist(),
        "K0_K0_prime_correlation": float(
            covariance[0, 1] / np.sqrt(covariance[0, 0] * covariance[1, 1])
        ),
        "curve_comparison": {
            "pressure_range_gpa": [0, 73.2],
            "grid_states": 1001,
            "max_abs_volume_difference_a3": float(
                np.max(np.abs(reconstructed - original))
            ),
            "max_abs_relative_volume_difference": float(
                np.max(np.abs(reconstructed / original - 1))
            ),
            "basis": "Calculated curves, not independent measured observations.",
        },
    }


def reproduce_takahashi() -> dict[str, Any]:
    """Replay the two source surfaces and one stated five-parameter fit setup.

    The authors name five fitted parameters, without specifying fitting order,
    weights or precise row mask. Simultaneous fitting of all 60 rows with equal
    pressure weights is an explicit diagnostic assumption. See the separate
    audit_takahashi_2019 replay for source-coordinate-width sensitivities.
    """
    with (
        ROOT / "peritheos/data/datasets/cementite-takahashi-2019-table-s1.csv"
    ).open() as stream:
        rows = list(csv.DictReader(stream))
    volume = np.array([float(r["volume_a3"]) for r in rows])
    temperature = np.array([float(r["temperature_k"]) for r in rows])
    results = {}
    for scale, values, errors in (
        ("kbr", [148.8, 311.1, 3.4, 1.06, 1.92], [1, 17.1, 0.10, 0.42, 1.73]),
        ("mgo", [147.3, 323, 3.43, 1.37, 0.98], [1, 16.6, 0.09, 0.33, 1.01]),
    ):
        pressure = np.array(
            [float(r[f"pressure_{scale}_gpa"] or r["pressure_sio2_gpa"]) for r in rows]
        )
        v0, k0, kp, gamma, q = values
        model = MieGruneisenDebye(
            BM3(V0=v0 * CELL_TO_MOLAR, K0=k0, K0_prime=kp),
            Tr=300,
            theta0=314,
            gamma0=gamma,
            q=q,
            n=4,
        )
        initial = dict(
            zip(
                ["rt_eos.V0", "rt_eos.K0", "rt_eos.K0_prime", "gamma0", "q"],
                [v0 * CELL_TO_MOLAR, k0, kp, gamma, q],
            )
        )
        fit = fit_joint_eos(
            MieGruneisenDebye,
            BM3,
            volume=volume * CELL_TO_MOLAR,
            temperature=temperature,
            pressure=pressure,
            initial=initial,
            fixed={"Tr": 300, "theta0": 314, "n": 4},
            bounds={
                "rt_eos.V0": (120 * CELL_TO_MOLAR, 180 * CELL_TO_MOLAR),
                "rt_eos.K0": (100, 600),
                "rt_eos.K0_prime": (1, 8),
                "gamma0": (0.01, 5),
                "q": (-5, 8),
            },
            absolute_sigma=False,
            max_nfev=5000,
        )
        comparisons = []
        for name, source, error in zip(initial, values, errors):
            factor = CELL_TO_MOLAR if name == "rt_eos.V0" else 1
            refit = float(fit.parameters[name]) / factor
            comparisons.append(
                {
                    "parameter": name,
                    "published": source,
                    "published_error": error,
                    "refit": refit,
                    "refit_error": float(fit.standard_errors[name]) / factor,
                    "difference": refit - source,
                    "within_printed_source_error_width": abs(refit - source) <= error,
                }
            )
        residual = model.pressure(volume * CELL_TO_MOLAR, temperature) - pressure
        cell_volume = float(model.volume(329, 5000) / CELL_TO_MOLAR)
        density = 4 * (3 * 55.845 + 12.011) / (Avogadro * cell_volume * 1e-24)
        results[scale] = {
            "coefficient_reconstruction": "similar",
            "source_equation_check": "consistent_with_source_density_benchmark",
            "original_author_workflow_parity": False,
            "selection": "all 60 Table S1 rows, including 2600 K row and ten SiO2-marker rows",
            "weighting": "equal pressure weights; source weights and row mask unknown",
            "observations": len(rows),
            "free_parameters": list(fit.free_parameters),
            "fixed_parameters": {"Tr": 300, "theta0": 314, "n": 4},
            "parameters": comparisons,
            "published_pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
            "published_max_abs_pressure_residual_gpa": float(np.max(np.abs(residual))),
            "refit_pressure_rmse_gpa": float(np.sqrt(np.mean(fit.residuals**2))),
            "covariance_order": list(fit.free_parameters),
            "covariance_model_units": np.asarray(fit.covariance).tolist(),
            "covariance_convention": "Residual-scaled inverse normal matrix; conditional on equal pressure weights. V0 uses molar J/bar.",
            "solver_success": bool(fit.success),
            "density_benchmark": {
                "pressure_gpa": 329,
                "temperature_k": 5000,
                "calculated_density_g_cm3": density,
                "source_value_g_cm3": 12.83 if scale == "kbr" else 12.64,
                "tolerance_g_cm3": 0.67 if scale == "kbr" else 0.04,
                "source": "Section 4, printed calculated density and its reported uncertainty"
                if scale == "kbr"
                else "Figure 3a, blue dotted 5000 K curve; manual graphical reading with 0.04 g/cm3 reading bound",
                "qualification": "Extrapolated source calculation beyond measured coverage; not a core-density measurement or author-fit reproduction.",
            },
        }
    return results


def reproduce() -> dict[str, Any]:
    return {
        "format": "peritheos.cementite-reproduction",
        "source_manifest": "peritheos/data/research/cementite/sources.json",
        "scott_2001": reproduce_scott(),
        "takahashi_2019": reproduce_takahashi(),
    }


def takahashi_equal_pressure_ledger_outcome(record: dict[str, Any]) -> dict[str, Any]:
    """Preserve the earlier equal-pressure diagnostic, separately from parity."""
    identifier = record["identifier"]
    scale = "kbr" if "_kbr_" in identifier else "mgo"
    result = reproduce_takahashi()[scale]
    parameters = []
    for row in result["parameters"]:
        value = row["published"]
        difference = row["difference"]
        name = row["parameter"].split(".")[-1]
        relative = abs(difference / value)
        numerical_similarity = (
            abs(difference) <= 1 or relative <= 0.2
            if name == "K0_prime"
            else relative <= {"V0": 0.05, "K0": 0.15}.get(name, 0.25)
        )
        parameters.append(
            {
                **row,
                "relative_difference": relative,
                "within_combined_2sigma": None,
                "similar": numerical_similarity
                or row["within_printed_source_error_width"],
            }
        )
    all_widths = all(row["within_printed_source_error_width"] for row in parameters)
    with (
        ROOT / "peritheos/data/datasets/cementite-takahashi-2019-table-s1.csv"
    ).open() as stream:
        rows = list(csv.DictReader(stream))
    pressure = [
        float(row[f"pressure_{scale}_gpa"] or row["pressure_sio2_gpa"]) for row in rows
    ]
    scale_vector = np.array([1 / CELL_TO_MOLAR, 1, 1, 1, 1])
    covariance = np.array(result["covariance_model_units"]) * np.outer(
        scale_vector, scale_vector
    )
    return {
        "status": "parity" if all_widths else "similar",
        "parity_basis": "Conditional simultaneous five-parameter comparison: every coefficient meets the broader numerical similarity screen or its printed source width, but parity additionally requires every central value within its printed source width. Missing source confidence and fit weights are not inferred. Source-equation density benchmarks are checked separately.",
        "qualification": "All 60 Table S1 rows and equal pressure weights are an explicit conditional reconstruction, including the 2600 K row and ten SiO2-marker rows. Source weights and exact fit selection remain unknown. No source coefficient is replaced and no author-workflow reproduction is claimed. The independent source density benchmarks concern extrapolated calculated states, not measurements or experimental validation.",
        "dataset_identifiers": ["cementite_takahashi_2019_table_s1"],
        "observations": result["observations"],
        "selection": result["selection"],
        "columns": {
            "pressure": f"pressure_{scale}_gpa|pressure_sio2_gpa",
            "volume": "volume_a3",
            "temperature": "temperature_k",
        },
        "observed_pressure_range_gpa": [min(pressure), max(pressure)],
        "observed_temperature_range_k": [300, 2600],
        "observed_volume_range": [
            min(float(r["volume_a3"]) for r in rows),
            max(float(r["volume_a3"]) for r in rows),
        ],
        "fit_kind": "joint_pvt",
        "objective": "pressure_residuals",
        "weighting": result["weighting"],
        "absolute_sigma": False,
        "free_parameters": result["free_parameters"],
        "fixed_parameters": ["Tr", "theta0", "n"],
        "fixed_values": result["fixed_parameters"],
        "parameters": parameters,
        "parameter_covariance": covariance.tolist(),
        "covariance_convention": "Residual-scaled inverse normal matrix in free_parameters order, with V0 converted to conventional-cell A3; conditional on equal pressure weights.",
        "published_rmse_gpa": result["published_pressure_rmse_gpa"],
        "rmse_gpa": result["refit_pressure_rmse_gpa"],
        "degrees_of_freedom": result["observations"] - len(result["free_parameters"]),
        "solver_success": result["solver_success"],
        "original_author_workflow_parity": False,
        "reproduction_report": "docs/literature-reproductions/cementite.md",
        "source_equation_density_benchmark": result["density_benchmark"],
    }


def takahashi_ledger_outcome(record: dict[str, Any]) -> dict[str, Any]:
    """Replay current coordinate-width parity and retain the earlier objective."""
    from scripts.audit_takahashi_2019 import parity_ledger_outcome

    return parity_ledger_outcome(
        record, takahashi_equal_pressure_ledger_outcome(record)
    )


if __name__ == "__main__":
    report = reproduce()
    REPORT.write_text(json.dumps(report, indent=2, allow_nan=False) + "\n")
    print(
        "Replayed Scott and both Takahashi source models; no source coefficients replaced."
    )
