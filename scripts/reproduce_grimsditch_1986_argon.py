"""Audit the printed EOS and explicitly non-independent refit diagnostics."""

from __future__ import annotations

import csv
import json
from pathlib import Path

import numpy as np
from scipy.optimize import brentq

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
RECORD = "argon_fcc_grimsditch_1986_density_polynomial"
MASS_FACTOR = 4 * 39.948 / 6.02214076e23 * 1e24


def pressure(rho):
    """Independent literal transcription of Eq.(6), GPa with rho in g/cm^3."""
    return 12.65 - 11.43 * rho + 1.5 * rho**2 + 0.68 * rho**3


def bulk(rho):
    """Eq.(7), isothermal derivative; no adiabatic correction is applied."""
    return rho * (-11.43 + 3.0 * rho + 2.04 * rho**2)


def equal_weight_diagnostic(rho, reported, source_rows):
    """Minimize equal-weight pressure residuals; this is not independent PV data."""
    coefficients = np.polynomial.polynomial.polyfit(rho, reported, 3)
    fitted = np.polynomial.polynomial.polyval(rho, coefficients)
    grid = np.linspace(min(rho), max(rho), 1001)
    curve_difference = np.polynomial.polynomial.polyval(grid, coefficients) - pressure(
        grid
    )
    return {
        "status": "non_independent_consistency_refit",
        "objective": "sum((P_polynomial(rho_i)-P_reported_i)**2); equal weight per printed row",
        "source_rows": list(map(int, source_rows)),
        "row_count": len(rho),
        "coefficients_c0_to_c3": list(map(float, coefficients)),
        "pressure_rms_gpa": float(np.sqrt(np.mean((fitted - reported) ** 2))),
        "published_pressure_rms_same_rows_gpa": float(
            np.sqrt(np.mean((pressure(rho) - reported) ** 2))
        ),
        "max_curve_difference_from_published_gpa_on_fitted_density_range": float(
            abs(curve_difference).max()
        ),
        "density_range_g_cm3": [float(min(rho)), float(max(rho))],
    }


def reproduce():
    doc = get_material_document("argon_fcc")
    eos = Material.from_eosmat(doc).get_eos_record(RECORD)
    with (
        ROOT / "peritheos/data/datasets/argon-grimsditch-1986-table1.csv"
    ).open() as stream:
        table = list(csv.DictReader(stream))
    solid = [r for r in table if r["phase"] == "fcc"]
    rho = np.array([float(r["density_g_cm3"]) for r in solid])
    reported = np.array([float(r["pressure_gpa"]) for r in solid])
    calc = pressure(rho)
    volumes = MASS_FACTOR / rho
    native_error = np.asarray(eos.pressure(volumes)) - calc
    residual = calc - reported
    source_rows = np.array([int(r["source_row"]) for r in solid])
    retained = source_rows != 98
    refits = {
        "all_75_solid_rows": equal_weight_diagnostic(rho, reported, source_rows),
        "sensitivity_without_inconsistent_source_row_98": equal_weight_diagnostic(
            rho[retained], reported[retained], source_rows[retained]
        ),
        "exclusion_reason": "Only sensitivity analysis: row98 is retained in the source dataset and primary all-row fit. Its printed density disagrees with the printed equation by1.216GPa.",
        "interpretation": "No independent regression reproduction: TableI densities were adopted from the prior EOS, and cover only the acoustic range, not the original diffraction data to77GPa.",
    }
    # Compare printed nv and C against arithmetic, without treating derived C as input.
    shift = np.array([float(r["brillouin_shift_cm_inverse"]) for r in table])
    nv = np.array([float(r["index_times_velocity_km_s"]) for r in table])
    density = np.array([float(r["density_g_cm3"]) for r in table])
    index = np.array([float(r["refractive_index"]) for r in table])
    modulus = np.array([float(r["effective_longitudinal_modulus_gpa"]) for r in table])
    nv_calculated = shift * 514.5e-7 * 299792.458 / 2
    modulus_calculated = density * (nv / index) ** 2
    with (
        ROOT / "peritheos/data/datasets/argon-fcc-grimsditch-1986-table2.csv"
    ).open() as stream:
        elastic = list(csv.DictReader(stream))
    checks = []
    for r in elastic:
        p = float(r["pressure_gpa"])
        density_at_p = brentq(lambda x: pressure(x) - p, 2, 5)
        b = bulk(density_at_p)
        checks.append(
            {
                "pressure_gpa": p,
                "equation_bulk_gpa": b,
                "printed_bulk_gpa": float(r["bulk_modulus_gpa"]),
                "printed_bulk_error_gpa": float(r["bulk_modulus_error_gpa"]),
            }
        )
    return {
        "record_identifier": RECORD,
        "equation_status": "reproduced",
        "fit_status": "not_reproduced",
        "summary_status": "approximate_internal_consistency_reproduced",
        "display_label": "Approximate internal consistency reproduced with equal weighting",
        "summary_scope": "74-row sensitivity excludes inconsistent source row 98; derived densities, not independent diffraction data. All-75-row diagnostic retained.",
        "reason": "Original X-ray fit rows not recovered; equal weighting is feasible but Table I densities and Table II bulk moduli depend on adopted EOS.",
        "equal_weight_diagnostics": refits,
        "table1_count": len(table),
        "liquid_count": len(table) - len(solid),
        "fcc_count": len(solid),
        "table2_count": len(elastic),
        "native_max_abs_pressure_error_gpa": float(abs(native_error).max()),
        "native_max_abs_bulk_error_gpa": float(
            abs(np.asarray(eos.eos.bulk_modulus(volumes)) - bulk(rho)).max()
        ),
        "native_volume_roundtrip_max_abs_a3": float(
            abs(np.asarray(eos.volume(calc)) - volumes).max()
        ),
        "derived_density_consistency_rms_gpa": float(np.sqrt(np.mean(residual**2))),
        "derived_density_consistency_max_abs_gpa": float(abs(residual).max()),
        "largest_density_discrepancy": {
            "source_row": int(solid[int(abs(residual).argmax())]["source_row"]),
            "note": "Retain source values and coefficients; do not repair through a refit.",
        },
        "nv_arithmetic_max_abs_km_s": float(abs(nv_calculated - nv).max()),
        "modulus_arithmetic_max_abs_gpa": float(
            abs(modulus_calculated - modulus).max()
        ),
        "table2_checks": checks,
        "reference_density_g_cm3": 2.0,
        "reference_pressure_gpa": pressure(2.0),
        "reference_volume_a3": MASS_FACTOR / 2,
        "reference_temperature": "room temperature; numeric 298 K is library convention only",
        "independent_observed_pv_count": 0,
    }


if __name__ == "__main__":
    result = reproduce()
    target = ROOT / "docs/data/argon-grimsditch-1986-reproduction.json"
    target.write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))
