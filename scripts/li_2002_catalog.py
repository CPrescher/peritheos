"""Catalog mapping of Li (2002); original coefficients remain authoritative."""

from __future__ import annotations

import copy

from scripts.reproduce_li_2002_fe3c import (
    SOURCE,
    V0,
    WIDTHS,
    arrays,
    curve_metrics,
    fit,
    graphical_checkpoints,
    rows,
)
from scripts.standardize_reconstruction_metadata import assessment_for_record

IDENTIFIER = "cementite_li_2002_bm3"
DATASET = "cementite_li_2002_table1"
REPORT = "docs/literature-reproductions/li-2002-fe3c.md"
DOI = "10.1007/s00269-001-0224-4"
REFERENCE = {
    "authors": [
        "J. Li",
        "H. K. Mao",
        "Y. Fei",
        "E. Gregoryanz",
        "M. Eremets",
        "C. S. Zha",
    ],
    "year": 2002,
    "title": "Compression of Fe3C to 30 GPa at room temperature",
    "source": "Physics and Chemistry of Minerals",
    "volume": "29",
    "issue": "3",
    "locator": "166-169",
    "doi": DOI,
}
QUALIFICATION = (
    "Conditional parameter parity compares K0 and K0_prime within the printed "
    "source error widths using all 12 nonambient Table 1 observations. Measured "
    "V0=155.28(11) A^3 per conventional Pnma Fe3C cell (Z=4) is fixed as in the "
    "source, not independently recovered by the fit. Published errors have "
    "unspecified confidence. The diagnostic effective-variance objective treats "
    "printed pressure and volume widths as if independent standard errors; its "
    "unscaled inverse-JtJ covariance is conditional. Original residual direction, "
    "weights, covariance, confidence convention and fit project are unavailable. "
    "Fixed-V0 uncertainty and shared calibration errors are not propagated. "
    "Original ruby wavelengths and separate ruby/neon pressure arrays are missing, "
    "so the pressure averaging cannot be replayed. Original-scale published "
    "coefficients and observations are preserved; Litasov (2013) recalibrated "
    "comparison coefficients are not substituted. This is not author-workflow "
    "reproduction, magnetic-phase validation or a recommendation for extrapolation."
)


def ledger_outcome(record):
    """Compute the declared conditional comparison from the literal source rows."""
    result = fit("effective_variance")
    selected, (p, v, _, _) = arrays()
    source = curve_metrics(SOURCE, p, v)
    comparison = []
    for i, key in enumerate(("K0", "K0_prime")):
        value = result["parameters"][i]
        delta = value - SOURCE[i]
        comparison.append(
            {
                "parameter": key,
                "published": float(SOURCE[i]),
                "published_error": float(WIDTHS[i]),
                "refit": value,
                "refit_error": result["conditional_absolute_widths"][i],
                "difference": delta,
                "relative_difference": abs(delta / SOURCE[i]),
                "within_published_error_width": bool(abs(delta) <= WIDTHS[i]),
                "within_published_printing_precision": bool(
                    abs(delta) <= (0.5, 0.05)[i]
                ),
                "similar": bool(abs(delta) <= WIDTHS[i]),
                "within_combined_2sigma": None,
                "uncertainty_finding": "conditional_width_compatible_at_printed_precision_source_confidence_unknown",
            }
        )
    if not all(row["within_published_error_width"] for row in comparison):
        raise AssertionError("Li conditional parameter comparison diverged")
    return {
        "status": "parity",
        "parity_basis": "within_reported_parameter_errors",
        "comparison_scope": "parameter_agreement",
        "conditional": True,
        "compared_parameters": ["K0", "K0_prime"],
        "author_workflow_reproduced": False,
        "original_author_workflow_parity": False,
        "dataset_identifiers": [DATASET],
        "fit_kind": "isothermal_pv",
        "objective": "effective_variance_pressure_residuals",
        "selection": "Table 1 rows 3-14; reference-card row 1 excluded; measured ambient row 2 supplies fixed V0 only",
        "source_row_ids": [r["row_id"] for r in selected],
        "columns": {
            "pressure": "pressure_gpa",
            "volume": "volume_a3",
            "temperature": None,
        },
        "weighting": "rP/sqrt(sP^2+(dP/dV*sV)^2), parameter-dependent, with printed widths conditionally treated as independent standard errors; no log-variance term",
        "absolute_sigma": True,
        "free_parameters": ["K0", "K0_prime"],
        "fixed_parameters": ["V0"],
        "fixed_values": {"V0": V0},
        "parameters": comparison,
        "observations": 12,
        "degrees_of_freedom": 10,
        "observed_pressure_range_gpa": [float(min(p)), float(max(p))],
        "observed_volume_range": [float(min(v)), float(max(v))],
        "observed_temperature_range_k": [300, 300],
        "solver_success": result["optimizer_success"],
        "rmse_gpa": result["rms_pressure_gpa"],
        "published_rmse_gpa": source["rms_pressure_gpa"],
        "rmse_volume_a3": result["rms_volume_a3"],
        "published_rmse_volume_a3": source["rms_volume_a3"],
        "parameter_covariance": result["conditional_absolute_covariance"],
        "covariance_convention": "Conditional unscaled inverse JtJ for the declared effective-variance residual; source confidence and shared covariance unknown, fixed-V0 uncertainty unpropagated",
        "reduced_chi_square": None,
        "objective_sum_squares_per_dof": result["objective_sum_squares_per_dof"],
        "qualification": QUALIFICATION,
        "reproduction_report": REPORT,
    }


def catalog_record():
    record = {
        "identifier": IDENTIFIER,
        "label": "Li et al. (2002), 300 K BM3",
        "model": "BM3",
        "temperature_ref": 300,
        "parameters": {
            "V0": {"value": V0, "error": 0.11, "fixed": True},
            "K0": {"value": 174, "error": 6},
            "K0_prime": {"value": 4.8, "error": 0.8},
        },
        "parameter_error_confidence": None,
        "experimental_pressure_range_gpa": [0.0001, 30.5],
        "experimental_temperature_range_k": [299, 301],
        "pressure_range_status": "reported_exactly",
        "record_kind": "published",
        "determination_method": "experimental",
        "fit_datasets": [DATASET],
        "reference": copy.deepcopy(REFERENCE),
        "notes": QUALIFICATION
        + " Eq. (2) approximates the 1-bar reference by zero pressure. The source is a 300 K isotherm; the 299-301 K range records reported room-temperature variation, not a thermal model. No crystallographic transition was observed to 30.5 GPa; the magnetic state was not measured. The material's structural reference cell is independent of this EOS V0.",
        "pressure_calibration": {
            "status": "resolved",
            "audit_date": "2026-10-30",
            "methods": [
                {
                    "kind": "ruby_fluorescence",
                    "source_location": "p. 167, Eq. (1)",
                    "reference": {
                        "authors": ["Zha", "Mao", "Hemley"],
                        "year": 2000,
                        "title": "Elasticity of MgO and a primary pressure scale to 55 GPa",
                        "source": "PNAS 97, 13494-13499",
                    },
                    "parameters": {
                        "A_gpa": 1904,
                        "B": 7.715,
                        "reference_wavelength_nm": 694.3,
                    },
                    "equation": "P=(A/B)*((1+delta_lambda/lambda0)^B-1)",
                    "notes": "Original Li/Zha constants; no approximate bundled ruby-scale substitution.",
                },
                {
                    "kind": "equation_of_state",
                    "source_location": "p. 167, Results and discussion; paired Ne lattice parameters in Table 1",
                    "reference": {
                        "authors": ["Hemley", "Zha", "Jephcoat", "Mao", "Finger"],
                        "year": 1989,
                        "title": "X-ray diffraction and equation of state of solid neon to 110 GPa",
                        "source": "Physical Review B 39, 11820-11827",
                    },
                    "scope": "Neon-bearing rows at and above 7.9 GPa",
                    "notes": "Use original Hemley (1989) calibration in source pressure average; no link to a later Fei-recalibrated or independently refitted neon EOS is asserted.",
                },
            ],
            "recalculation": {
                "status": "missing_calibrant_observations",
                "notes": "Ne lattice parameters retained, but original ruby wavelengths and separate ruby/neon pressure arrays are absent; no pressure recalculation applied.",
            },
            "notes": "Neon medium. Source uses the average of neon and ruby pressures where neon diffraction is available. Table 1 retains that original reduced pressure coordinate. Magnetic stability and hydrostatic accuracy are not independently validated.",
        },
        "evidence": {
            "parameter_provenance": {
                "determination_method": "p. 166-168: polycrystalline XRD compression and least-squares BM3 fit at 300 K; measured ambient V0 constrained.",
                "V0": "Table 1 measured ambient row 2: 155.28(11) A^3; fixed choice stated in abstract and p. 168 beside Eq. (2).",
                "K0": "p. 168, paragraph following Eq. (2), and abstract: 174 +/- 6 GPa.",
                "K0_prime": "p. 168, paragraph following Eq. (2), and abstract: 4.8 +/- 0.8.",
                "reference_state": "1 bar and 300 K; Eq. (2) sets P(V0)=0 as the negligible-1-bar approximation. Conventional Pnma Fe3C cell, Z=4; no conversion of source cell volumes.",
                "uncertainty": "Printed widths retained; confidence and source covariance not reported. V0 error is a measured anchor error despite its fixed fit status.",
            },
            "scientific_validation": {
                "status": "primary_source_validated",
                "audit_date": "2026-10-30",
                "note": "Primary equation, original coefficients/errors, fixed measured V0 and original pressure basis verified. Five source-curve checkpoints pass their declared 0.12 A^3 reading bound; conditional K0/K0_prime parity is separate from author-workflow and physical validation.",
                "verified_fields": [
                    "equation",
                    "parameters",
                    "reported_errors",
                    "fixed_V0",
                    "units",
                    "reference_state",
                    "volume_basis",
                    "source_curve_checkpoints",
                    "primary_table",
                    "pressure_basis",
                ],
                "primary_source_check": {
                    "doi": DOI,
                    "access_url": "https://doi.org/" + DOI,
                    "locations": [
                        "pp. 166-168, abstract, Table 1, Eq. (1)-(2), Figure 2"
                    ],
                    "source_pdf_sha256": "3b4761baa81a3e4bd419f8e93d5e453d7c5d7c2b64b400514d0f468dbcb211af",
                    "finding": "Published BM3 coefficients preserved with measured V0 fixed; direct source-equation and conditional coefficient checks pass. See "
                    + REPORT,
                },
                "primary_data_check": {
                    "status": "primary_table",
                    "dataset_identifiers": [DATASET],
                    "finding": "All 14 Table 1 rows, source tokens and errors retained. Row 1 is a reference card; row 2 is the measured anchor; rows 3-14 are the 12 nonambient fit observations. Ten paired neon lattice parameters retained.",
                },
                "unresolved": [QUALIFICATION],
            },
            "reproduction": {
                "conditional": True,
                "compared_parameters": ["K0", "K0_prime"],
                "author_workflow_reproduced": False,
                "source_equation": {
                    "status": "parity",
                    "scope": "source_equation",
                    "report": REPORT,
                    **graphical_checkpoints(),
                },
                "uncertainty_reproduction": "Conditional effective-variance absolute widths reproduce source printing precision; original confidence, covariance, anchor-error propagation and workflow remain unknown. Residual-scaled widths are smaller and not source errors.",
            },
        },
    }
    outcome = ledger_outcome(record)
    outcome.update(material_identifier="cementite", record_identifier=IDENTIFIER)
    record["evidence"]["reproduction"]["assessment"] = assessment_for_record(
        record, outcome
    )
    return record


def catalog_dataset():
    raw, normalized = rows()
    columns = [
        {
            "name": "source_row",
            "quantity": "index",
            "unit": "dimensionless",
            "role": "value",
        },
        {
            "name": "observation_kind_code",
            "quantity": "source_row_kind",
            "unit": "dimensionless",
            "role": "flag",
        },
    ]
    fields = ["pressure_gpa", "volume_a3", "a_a", "b_a", "c_a", "neon_a_a"]
    for field in fields:
        quantity = (
            "pressure"
            if field == "pressure_gpa"
            else "conventional_unit_cell_volume"
            if field == "volume_a3"
            else "lattice_parameter"
        )
        unit = (
            "GPa"
            if field == "pressure_gpa"
            else "angstrom^3"
            if field == "volume_a3"
            else "angstrom"
        )
        columns.extend(
            [
                {"name": field, "quantity": quantity, "unit": unit, "role": "value"},
                {
                    "name": field + "_error_width",
                    "quantity": quantity,
                    "unit": unit,
                    "role": "uncertainty",
                    "of": field,
                },
            ]
        )
    return {
        "identifier": DATASET,
        "kind": "pressure_volume_lattice_with_reference_card",
        "description": "Complete original-scale Table 1. Row 1 is a reference-card comparison (kind 0), row 2 a measured ambient anchor (kind 1), rows 3-14 measured nonambient observations (kind 2). All printed pressure, lattice, volume and paired neon uncertainties retained; confidence unspecified.",
        "reference": copy.deepcopy(REFERENCE),
        "source_location": "Table 1 and footnotes a-b, printed p. 167, at 300 K",
        "columns": columns,
        "rows": [
            [
                i + 1,
                0 if i == 0 else 1 if i == 1 else 2,
                *[r[key] for f in fields for key in (f, f + "_error_width")],
            ]
            for i, r in enumerate(normalized)
        ],
        "used_by_eos_records": [IDENTIFIER],
        "source_tokens": raw,
        "uncertainty": {
            "reported_as": "parenthetical_last_digit_width_or_explicit_width; confidence_unspecified",
            "covariance": "not_reported; independence_not_established",
            "notes": "Widths reflect calibration variation and reproducibility, not established independent standard deviations. Ambient pressure token 0.000100 (1) retained literally.",
        },
        "source_rights": "Springer-Verlag 2002 article copyright; source numerical data license unspecified. No source text, PDF or figure is relicensed.",
        "contributor_rights": "CC0 only for contributor-held factual transcription, normalization, metadata and arrangement rights.",
        "provenance": {
            "type": "primary_table_transcription",
            "source_pdf_sha256": "3b4761baa81a3e4bd419f8e93d5e453d7c5d7c2b64b400514d0f468dbcb211af",
            "report": REPORT,
            "normalization": "Printed conventional-cell volumes used directly; a*b*c and Ne a^3 are not substituted as observations.",
        },
        "observation_views": [
            {
                "identifier": IDENTIFIER + "_observations",
                "eos_records": [IDENTIFIER],
                "evidence_kind": "experimental",
                "notes": "All 13 measured Table 1 states; reference-card row excluded. The declared 12-point fit excludes the measured ambient fixed-V0 anchor. Pressure is the original source average where Ne is present; no recalibration.",
                "row_selection": [
                    {"column": "observation_kind_code", "values": [1, 2]}
                ],
                "coordinates": {
                    "pressure": {"kind": "column", "column": "pressure_gpa"},
                    "volume": {
                        "kind": "column",
                        "column": "volume_a3",
                        "unit": "angstrom^3/conventional_unit_cell",
                        "reference": "Li (2002), Table 1: printed conventional Pnma Fe3C cell volume; no conversion",
                    },
                },
                "uncertainties": {
                    "pressure": "pressure_gpa_error_width",
                    "volume": "volume_a3_error_width",
                },
                "fit_selection": {
                    "row_selection": [
                        {"column": "observation_kind_code", "values": [2]}
                    ],
                    "notes": "Conditional diagnostic replay; original author weight/objective settings unavailable.",
                },
            }
        ],
    }
