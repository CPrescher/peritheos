"""Public evidence/report adapters for the archived author-processed pyrope fit.

Repository maintenance only. Reads frozen inputs/results; never refits or
replaces the published EOS/covariance or claims raw-study preprocessing parity.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "docs/data/garnet-angel-2022/pyrope"
RECORD = "pyrope_angel_2022_bm3_q_compromise_debye"
REPORT = "docs/literature-reproductions/pyrope-author-input-replay.md"
QUALIFICATIONS = [
    "All 194 author-processed observations, assigned P/T/value errors and scale groups are recovered. Original upstream reductions, calibration covariance and all exclusion decisions are not reconstructed.",
    "The comparison uses the published EOS snapshot; Ross's later final log and its EosFit replay are separately compared in docs/literature-reproductions/pyrope-author-input-replay.md. Neither snapshot replaces the published coefficients or covariance.",
    "The independent pressure-residual IRLS replay uses precise modulus inversion and local derivative effective variances. EosFit uses finite one-error displacements and approximate modulus inversion; the small parameter/error/chi-square differences are documented, not bitwise identity.",
    "Refit errors are conditional local standard errors rescaled by reduced chi-square with effective coordinate weights frozen at the solution. Printed source errors retain their original convention; raw calibration/systematic uncertainty and shared-data cross-fit covariance remain unknown.",
    "V0, Tr, n and volume scales are fixed; only the 15 high-pressure Sinogeikin–Bass (2000) modulus points share the refined scale. The two ambient points are unscaled. The nuisance scale is not a new EOS coefficient.",
    "Numerical parameter agreement and processed-input replay do not establish complete raw-study author-workflow reproduction, experimental accuracy or a recommendation/default change.",
]


def outcome(record):
    from scripts.audit_garnet_angel_2022 import source

    audit = json.loads(
        (DATA / "author-2026-11-14-audit.json").read_text(encoding="utf-8")
    )
    fit = audit["independent_fits"][0]
    p, covariance = source("pyrope")
    mapping = [
        ("rt_eos.K0", "K0_gpa", p[2], record["parameter_errors"]["K0"]),
        ("rt_eos.K0_prime", "K0_prime", p[3], record["parameter_errors"]["K0_prime"]),
        ("theta0", "theta0_k", p[11], record["thermal"]["parameter_errors"]["theta0"]),
        ("gamma0", "gamma0", p[18], record["thermal"]["parameter_errors"]["gamma0"]),
        (
            "scale_sinogeikin_2000_ks",
            "scale_sinogeikin_2000_ks",
            p[53],
            float(np.sqrt(covariance[52, 52])),
        ),
    ]
    comparisons = []
    for name, key, published, error in mapping:
        refit = fit["parameters"][key]
        width = fit["conditional_residual_scaled_widths"][key]
        comparisons.append(
            dict(
                parameter=name,
                published=published,
                published_error=error,
                refit=refit,
                refit_error=width,
                difference=refit - published,
                relative_difference=abs(refit - published) / abs(published),
                within_combined_2sigma=bool(
                    abs(refit - published) <= 2 * np.hypot(error, width)
                ),
                similar=abs(refit - published) / abs(published) < 0.01,
            )
        )
    if not all(p["within_combined_2sigma"] and p["similar"] for p in comparisons):
        raise ValueError(
            "Archived pyrope replay no longer meets the recorded parameter criterion"
        )
    rows = fit["rows"]
    return dict(
        status="parity",
        primary_data_status="author_processed_inputs_recovered",
        primary_data_finding=QUALIFICATIONS[0],
        assessment_basis="Conditional parameter agreement from the complete 194-observation author-processed input regression. All five free parameters meet the recorded combined-error and numerical similarity criteria; the separate final-log comparison differs by less than 0.01 reported standard errors. Published source coefficients/covariance remain unchanged.",
        assessment_qualifications=QUALIFICATIONS,
        qualification=" ".join(QUALIFICATIONS),
        dataset_identifiers=[
            "pyrope_angel_2022_author_pvt",
            "pyrope_angel_2022_author_pkt",
        ],
        fit_kind="author_processed_joint_volume_modulus_q_compromise",
        reproduction_report=REPORT,
        objective="Pressure residuals with iteratively frozen effective variances from all supplied P/T/value errors; separate from the public latent-coordinate fitter",
        selection="All 194 rows of docs/data/garnet-angel-2022/pyrope/author-2026-11-14-inputs.csv: 115 volumes and 79 adiabatic moduli; no further exclusions",
        weighting="Independent coordinate contributions in quadrature; local derivative propagation in the Python replay; original error tokens and exact scale groups retained",
        observations=194,
        degrees_of_freedom=189,
        free_parameters=[r[0] for r in mapping],
        fixed_parameters=["V0", "Tr", "n", "volume_scales", "unscaled_modulus_group"],
        fixed_values={
            "V0_cm3_mol": 113.13,
            "Tr_k": 298.0,
            "n": 20.0,
            "volume_scales": 1.0,
            "unscaled_modulus_group": 1.0,
        },
        solver_success=fit["numerically_converged"],
        solver_message="Two independent starts converge to the archived solution without reaching bounds",
        parameters=comparisons,
        rmse_gpa=float(
            np.sqrt(
                np.mean(
                    [r["calculated_minus_observed_pressure_gpa"] ** 2 for r in rows]
                )
            )
        ),
        reduced_chi_square=fit["conditional_reduced_chi_square"],
        parameter_covariance=(
            np.array(fit["conditional_unscaled_covariance"])
            * fit["conditional_reduced_chi_square"]
        ).tolist(),
        covariance_convention="Conditional inverse-information covariance with final effective coordinate weights frozen, multiplied by chi-square/189; parameter order follows free_parameters",
    )


def comparison():
    from scripts.audit_garnet_angel_2022 import source
    from scripts.reconstruction_diagnostics import DiagnosticComparison, panel
    from scripts.refit_pyrope_angel_2022_subset import arrays, pressure_at_observations

    audit = json.loads(
        (DATA / "author-2026-11-14-audit.json").read_text(encoding="utf-8")
    )
    fit = audit["independent_fits"][0]
    p, _ = source("pyrope")
    rows = fit["rows"]
    if len(rows) != 194:
        raise ValueError("Pyrope author-input selection must retain all 194 rows")
    d = arrays(rows)
    predictions = [
        pressure_at_observations(c, d)
        for c in (
            np.array([p[2], p[3], p[11], p[18], p[53]]),
            np.array(list(fit["parameters"].values())),
        )
    ]
    panels = []
    for kind, title, xlabel in [
        ("V", "Volume observations (115 rows)", "Author volume input (cm³/mol)"),
        (
            "KS",
            "Adiabatic-modulus observations (79 rows)",
            "Author modulus input (GPa)",
        ),
    ]:
        mask = np.array([r["kind"] == kind for r in rows])
        panels.append(
            panel(
                title,
                d["Y"][mask],
                d["P"][mask],
                *[prediction[mask] for prediction in predictions],
                xlabel,
                "Pressure (GPa)",
                "GPa",
                sigma=d["sP"][mask],
                x_sigma=d["sY"][mask],
                observed_label="Author-processed observations",
            )
        )
    return DiagnosticComparison(
        panels,
        "All 194 processed observations are shown at their individual temperatures, without connecting mixed-temperature points. The source curve uses the published EOS snapshot and source scale; the reconstruction uses the archived independent Python fit. Both reconstruct pressure from the same observed V or KS and T, including the source's operational KS conversion and exact scale groups. Residual = calculated − observed pressure (GPa). Panel RMS is unweighted and differs from the effective-variance chi-square; error bars show the supplied P/value errors, while temperature errors also enter the fit. Full input/error/group tokens and the later author-log/native comparison are linked in the dedicated replay report. Original raw-study reduction/exclusion history remains partial.",
        source_title="Published EOS snapshot",
    )
