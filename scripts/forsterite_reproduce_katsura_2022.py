#!/usr/bin/env python3
"""Audit author-script forsterite pressure, table comparison and ordinary refits.

This does not reconstruct the unpublished bootstrap realizations or fit-input
export. No refit coefficients are installed in the material catalog.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter
from pathlib import Path

import numpy as np
from scipy.optimize import least_squares

from peritheos import get_eos_record
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.units import cell_volume_to_molar_volume

ROOT = Path(__file__).resolve().parents[1]
TABLE_RECORD = "forsterite_katsura_2022_bm3_mgd"
RECORD = "forsterite_katsura_2022_author_script_bm3_mgd"
SCRIPT_PARAMETERS = (4.21116896988574, 0.997258584556838, 2.43854091371459)
DATA = ROOT / "peritheos/data/datasets/forsterite-katsura-2022-olivine-pvt.csv"
FIT_INPUT = DATA.with_name("forsterite-katsura-2022-candidate-fitting-input.csv")
LEDGER_PARITY_POLICY = (
    "Katsura 2022 forsterite accepts conditional numerical parity within printed "
    "Table 2 parameter error widths. The full-precision theta0=760 K author-script "
    "vector is recommended; the theta0=768 K table alternative is retained. "
    "Table errors are a comparison scale only for the script vector; script-specific "
    "errors, covariance and exact bootstrap reproduction remain unavailable."
)


def load_rows() -> list[dict[str, str]]:
    with DATA.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def author_pressure(
    ratio,
    temperature,
    parameters=SCRIPT_PARAMETERS,
    theta0=760.0,
    avogadro=6.02214076e23,
):
    """Independent translation of deposited MGD_EOS_RJ/ThermalEnergy.

    Defaults retain the geotherm script's full-precision vector, theta0,
    Avogadro constant, R and Debye polynomial. Explicit arguments can evaluate
    the table variant or the fitting script's older Avogadro constant.
    """
    ratio = np.asarray(ratio, dtype=float)
    temperature = np.asarray(temperature, dtype=float)
    kp, gamma0, q = parameters
    gamma = gamma0 * ratio**q
    theta = theta0 * np.exp((gamma0 - gamma) / q)
    df = [
        -0.0000013397592,
        0.000035022653,
        -0.00031034588,
        0.00018490174,
        0.016520035,
        -0.12495878,
        0.33333333,
    ]

    def energy(t):
        return 9 * 7 * 8.314472 * t * np.polyval(df, theta / t)

    molar_m3 = avogadro / 4 * 289.58 * ratio * 1e-30
    cold = 1.5 * 127.4 * (ratio ** (-7 / 3) - ratio ** (-5 / 3))
    cold *= 1 + 0.75 * (kp - 4) * (ratio ** (-2 / 3) - 1)
    return cold + gamma / molar_m3 * (energy(temperature) - energy(300.0)) * 1e-9


def author_ordinary_optimizer(ratio, temperature, pressure):
    """Translate NlLSq, including its .1% derivatives and 1e-4 step stop.

    This is the deterministic optimizer inside the bootstrap, applied once to
    unperturbed observations. Original temporary olivine guesses are retained.
    """
    parameters = np.array([4.0, 1.26, 2.9])
    damping = 0.01

    def prediction(values):
        return author_pressure(ratio, temperature, values, avogadro=6.0221367e23)

    chisq = float(np.sum((prediction(parameters) - pressure) ** 2))
    for iteration in range(1, 1001):
        derivatives = []
        for index in range(3):
            plus, minus = parameters.copy(), parameters.copy()
            plus[index] *= 1.001
            minus[index] *= 0.999
            derivatives.append(
                (prediction(plus) - prediction(minus)) / (plus[index] - minus[index])
            )
        jacobian = np.array(derivatives)
        hessian = jacobian @ jacobian.T
        hessian[np.diag_indices(3)] *= 1 + damping
        step = np.linalg.solve(hessian, jacobian @ (pressure - prediction(parameters)))
        proposed = parameters + step
        new_chisq = float(np.sum((prediction(proposed) - pressure) ** 2))
        if new_chisq < chisq:
            parameters, chisq = proposed, new_chisq
            damping /= 2
        else:
            damping *= 2
        if np.linalg.norm(step / parameters) <= 1e-4:
            return parameters, iteration
    raise RuntimeError("Author optimizer exceeded 1000 steps")


def rerun_reconstructed_input() -> dict:
    """Read the actual six-column candidate CSV; do not fill missing sigmaT."""
    with FIT_INPUT.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    with DATA.with_name("forsterite-katsura-2009-table1.csv").open(
        newline="", encoding="utf-8"
    ) as stream:
        table = list(csv.DictReader(stream))
    matched_ids = {r["workbook_run_id"] for r in table if r["workbook_run_id"]}
    source_rows = load_rows()
    assert len(rows) == len(source_rows)
    ratio = np.array([float(r["volume_ratio"]) for r in rows])
    temperature = np.array([float(r["temperature_k"]) for r in rows])
    pressure = np.array([float(r["pressure_gpa"]) for r in rows])
    v0 = float(cell_volume_to_molar_volume(289.58, 4))

    def native_prediction(parameters, mask):
        model = MieGruneisenDebye(
            BM3(v0, 127.4, parameters[0]),
            300,
            760,
            parameters[1],
            parameters[2],
            7,
        )
        return model.pressure(v0 * ratio[mask], temperature[mask])

    report = {
        "input": str(FIT_INPUT.relative_to(ROOT)),
        "input_sha256": hashlib.sha256(FIT_INPUT.read_bytes()).hexdigest(),
        "objective": "Unweighted squared pressure residuals, as in author ChisqPressure; errors are not weights in that optimizer.",
        "fixed": {
            "V0_a3": 289.58,
            "K0_gpa": 127.4,
            "Tr_k": 300,
            "theta0_k": 760,
            "n": 7,
            "Z": 4,
        },
        "author_script_parameters": dict(
            zip(["K0_prime", "gamma0", "q"], SCRIPT_PARAMETERS)
        ),
        "missing_temperature_sigma_count": sum(
            r["temperature_sigma_k"] == "" for r in rows
        ),
        "bootstrap_performed": False,
        "selections": [],
    }
    for name, mask in [
        ("all_51_deposited_rows", np.ones(len(rows), dtype=bool)),
        (
            "50_rows_matched_to_2009_table1",
            np.array([r["run_id"] in matched_ids for r in source_rows]),
        ),
    ]:

        def author_prediction(values):
            return author_pressure(
                ratio[mask], temperature[mask], values, avogadro=6.0221367e23
            )

        def native_for_selection(values):
            return native_prediction(values, mask)

        methods = []
        for mode in [
            "literal_author_equation_scipy_optimizer",
            "exact_debye_scipy_optimizer",
            "literal_author_equation_author_optimizer",
        ]:
            if mode == "literal_author_equation_author_optimizer":
                fitted, iterations = author_ordinary_optimizer(
                    ratio[mask], temperature[mask], pressure[mask]
                )
                evaluate = author_prediction
            else:
                evaluate = (
                    author_prediction
                    if mode.startswith("literal")
                    else native_for_selection
                )
                fit = least_squares(
                    lambda values: evaluate(values) - pressure[mask],
                    SCRIPT_PARAMETERS,
                    xtol=1e-13,
                    ftol=1e-13,
                    gtol=1e-13,
                )
                assert fit.success
                fitted = fit.x
                iterations = fit.nfev
            residual = evaluate(fitted) - pressure[mask]
            methods.append(
                {
                    "method": mode,
                    "parameters": dict(
                        zip(["K0_prime", "gamma0", "q"], fitted.tolist())
                    ),
                    "difference_from_author_script": dict(
                        zip(
                            ["K0_prime", "gamma0", "q"],
                            (fitted - SCRIPT_PARAMETERS).tolist(),
                        )
                    ),
                    "relative_difference_percent": dict(
                        zip(
                            ["K0_prime", "gamma0", "q"],
                            (
                                100 * (fitted - SCRIPT_PARAMETERS) / SCRIPT_PARAMETERS
                            ).tolist(),
                        )
                    ),
                    "pressure_rmse_gpa": float(np.sqrt(np.mean(residual**2))),
                    "author_script_pressure_rmse_gpa": float(
                        np.sqrt(
                            np.mean((evaluate(SCRIPT_PARAMETERS) - pressure[mask]) ** 2)
                        )
                    ),
                    "max_curve_difference_from_author_script_gpa": float(
                        np.max(abs(evaluate(fitted) - evaluate(SCRIPT_PARAMETERS)))
                    ),
                    "iterations_or_function_evaluations": iterations,
                }
            )
        report["selections"].append(
            {
                "selection": name,
                "row_count": int(mask.sum()),
                "excluded_run_ids": [
                    r["run_id"] for r, keep in zip(source_rows, mask) if not keep
                ],
                "fits": methods,
            }
        )
    # Table 2 error widths are a comparison scale, not script-specific errors.
    document = json.loads(
        (ROOT / "peritheos/data/materials/forsterite.eosmat").read_text(
            encoding="utf-8"
        )
    )
    table_record = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == TABLE_RECORD
    )
    widths = {
        "K0_prime": table_record["parameter_errors"]["K0_prime"],
        "gamma0": table_record["thermal"]["parameter_errors"]["gamma0"],
        "q": table_record["thermal"]["parameter_errors"]["q"],
    }
    main_fit = report["selections"][0]["fits"][0]
    comparisons = [
        {
            "parameter": name,
            "author_script": source,
            "refit": main_fit["parameters"][name],
            "reported_error_width": widths[name],
            "difference_over_reported_error_width": abs(
                main_fit["difference_from_author_script"][name]
            )
            / widths[name],
            "within_reported_error": abs(
                main_fit["difference_from_author_script"][name]
            )
            <= widths[name],
            "similar": (
                abs(main_fit["difference_from_author_script"][name]) <= 1.0
                or abs(main_fit["relative_difference_percent"][name]) <= 20.0
            )
            if name == "K0_prime"
            else abs(main_fit["relative_difference_percent"][name]) <= 25.0,
        }
        for name, source in report["author_script_parameters"].items()
    ]
    report["parameter_reproduction"] = {
        "status": "parity"
        if all(p["within_reported_error"] and p["similar"] for p in comparisons)
        else "parity_not_achieved",
        "parity_basis": "within_reported_parameter_errors",
        "scope": "Conditional numerical parity of the three free author-script coefficients using all 51 reconstructed observations and the author equation.",
        "uncertainty_reference_record": TABLE_RECORD,
        "qualification": "Printed Table 2 bootstrap standard deviations are used as comparison widths only. They accompany the rounded theta0=768 K parameterization and are not assigned as uncertainties of the theta0=760 K script vector. No combined-two-sigma, covariance or exact bootstrap reproduction is asserted.",
        "parameters": comparisons,
    }
    return report


def diagnostics() -> dict:
    rows = load_rows()

    def column(key):
        return np.array([float(row[key]) for row in rows])

    ratio = column("forsterite_volume_ratio")
    sigma_ratio = column("forsterite_volume_ratio_sigma")
    temperature = column("temperature_corrected_k")
    pressure = column("pressure_mean_gpa")
    volume = 289.58 * ratio
    record = get_eos_record(RECORD)
    predicted = record.pressure(volume, temperature, check_validity=False)
    source_predicted = author_pressure(ratio, temperature)
    projected_sigma = np.hypot(
        column("pressure_spread_gpa"),
        record.eos.bulk_modulus(volume * record.volume_scale, temperature)
        * sigma_ratio
        / ratio,
    )
    result = {
        "record": RECORD,
        "source_publication": "10.1029/2021JB023562",
        "source_data": "10.5281/zenodo.5644426",
        "source_code": "10.5281/zenodo.5903286",
        "row_count": len(rows),
        "pressure_rmse_gpa": float(np.sqrt(np.mean((predicted - pressure) ** 2))),
        "author_polynomial_max_difference_gpa": float(
            np.max(abs(predicted - source_predicted))
        ),
        "max_residual_over_projected_spread": float(
            np.max(abs(predicted - pressure) / projected_sigma)
        ),
        "checkpoints": [],
        "ordinary_refits": [],
        "qualification": "Author-script coefficients recommended; table alternative retained. Equation validated; conditional numerical parameter parity within printed Table 2 error widths. Exact bootstrap workflow reproduction remains unresolved. Five candidate fitting-input columns are recoverable from the deposited workbook; its missing temperature errors, exact export selection and bootstrap realizations are not invented.",
    }
    for index in [0, len(rows) - 1]:
        result["checkpoints"].append(
            {
                "source_row": int(rows[index]["source_row"]),
                "run_id": rows[index]["run_id"],
                "volume_ratio": float(ratio[index]),
                "temperature_k": float(temperature[index]),
                "observed_pressure_gpa": float(pressure[index]),
                "calculated_pressure_gpa": float(predicted[index]),
                "projected_spread_gpa": float(projected_sigma[index]),
            }
        )
    worst = int(np.argmax(abs(predicted - pressure) / projected_sigma))
    result["largest_normalized_residual"] = {
        "source_row": int(rows[worst]["source_row"]),
        "run_id": rows[worst]["run_id"],
        "temperature_k": float(temperature[worst]),
        "observed_pressure_gpa": float(pressure[worst]),
        "calculated_minus_observed_pressure_gpa": float(
            predicted[worst] - pressure[worst]
        ),
        "projected_spread_gpa": float(projected_sigma[worst]),
        "disposition": "Retained without exclusion or correction in the 51-row main fit. The recovered 2009 Table 1 does not contain this observation; a separate 50-row table-matched sensitivity fit is reported.",
    }
    v0_molar = float(cell_volume_to_molar_volume(289.58, 4))
    for theta0 in [768.0, 760.0]:

        def residual(parameters):
            model = MieGruneisenDebye(
                BM3(v0_molar, 127.4, parameters[0]),
                Tr=300,
                theta0=theta0,
                gamma0=parameters[1],
                q=parameters[2],
                n=7,
            )
            return model.pressure(v0_molar * ratio, temperature) - pressure

        fit = least_squares(
            residual,
            [4.2, 1.0, 2.4],
            bounds=([0, 0.01, -20], [15, 5, 20]),
            xtol=1e-12,
            ftol=1e-12,
            gtol=1e-12,
        )
        result["ordinary_refits"].append(
            {
                "theta0_k": theta0,
                "parameters": dict(zip(["K0_prime", "gamma0", "q"], fit.x.tolist())),
                "rmse_gpa": float(np.sqrt(np.mean(fit.fun**2))),
                "success": bool(fit.success),
                "objective": "Unweighted sum of squared pressure residuals; AD, AN, U; all 51 rows; V0, K0, theta0, Tr, n fixed. No bootstrap or invented temperature errors.",
            }
        )
    table = get_eos_record(TABLE_RECORD)
    table_pressure = table.pressure(volume, temperature)
    result["table_comparison"] = {
        "record": TABLE_RECORD,
        "pressure_rmse_gpa": float(np.sqrt(np.mean((table_pressure - pressure) ** 2))),
        "max_script_minus_table_pressure_gpa": float(
            np.max(abs(predicted - table_pressure))
        ),
        "author_polynomial_max_difference_gpa": float(
            np.max(
                abs(
                    table_pressure
                    - author_pressure(
                        ratio,
                        temperature,
                        parameters=(4.2, 1.0, 2.4),
                        theta0=768,
                        avogadro=6.0221367e23,
                    )
                )
            )
        ),
    }
    result["reconstructed_input_refits"] = rerun_reconstructed_input()
    result["fit_reproduction"] = result["reconstructed_input_refits"][
        "parameter_reproduction"
    ]
    return result


def ledger_outcome(record: dict, report: dict | None = None) -> dict:
    """Return the dedicated conditional result for the global EOS ledger."""
    report = diagnostics() if report is None else report
    script_variant = record["identifier"] == RECORD
    if record["identifier"] not in {RECORD, TABLE_RECORD}:
        raise ValueError("Not a Katsura 2022 forsterite record")
    comparison = report["reconstructed_input_refits"]["parameter_reproduction"]
    if script_variant:
        fit = report["reconstructed_input_refits"]["selections"][0]["fits"][0]
        parameters = fit["parameters"]
        published_rmse = fit["author_script_pressure_rmse_gpa"]
        rmse = fit["pressure_rmse_gpa"]
    else:
        fit = next(r for r in report["ordinary_refits"] if r["theta0_k"] == 768)
        parameters = fit["parameters"]
        published_rmse = report["table_comparison"]["pressure_rmse_gpa"]
        rmse = fit["rmse_gpa"]
    widths = {
        p["parameter"]: p["reported_error_width"] for p in comparison["parameters"]
    }
    source = {
        "K0_prime": record["eos"]["parameters"]["K0_prime"],
        "gamma0": record["thermal"]["parameters"]["gamma0"],
        "q": record["thermal"]["parameters"]["q"],
    }
    comparisons = []
    for name, value in source.items():
        difference = parameters[name] - value
        relative = abs(difference / value)
        comparisons.append(
            {
                "parameter": name,
                "published": value,
                "published_error": None if script_variant else widths[name],
                "comparison_error_width": widths[name],
                "refit": parameters[name],
                "refit_error": None,
                "difference": difference,
                "relative_difference": relative,
                "similar": (abs(difference) <= 1.0 or relative <= 0.20)
                if name == "K0_prime"
                else relative <= 0.25,
                "within_reported_error": abs(difference) <= widths[name],
                "within_combined_2sigma": None,
            }
        )
    return {
        "status": "parity"
        if all(p["within_reported_error"] and p["similar"] for p in comparisons)
        else "parity_not_achieved",
        "parity_basis": "within_reported_parameter_errors",
        "parameters": comparisons,
        "dataset_identifiers": ["forsterite_katsura_2022_olivine_pvt"],
        "observations": 51,
        "fit_kind": "reconstructed_input_conditional_bm3_mgd",
        "objective": "unweighted_pressure_residuals",
        "fixed_parameters": ["V0", "K0", "Tr", "theta0", "n", "Z"],
        "free_parameters": list(source),
        "selection": "All 51 deposited observations; AD pressure, AN volume ratio, U temperature; M375006 retained.",
        "published_rmse_gpa": published_rmse,
        "rmse_gpa": rmse,
        "solver_success": True,
        "original_author_workflow_parity": False,
        "recommended_thermal_record": RECORD,
        "recommended": script_variant,
        "uncertainty_reference_record": TABLE_RECORD,
        "qualification": LEDGER_PARITY_POLICY
        + " Original 2009 combined-fit reproduction is a separate unresolved audit.",
        "reason": (
            "Recommended thermal option: full-precision author-script coefficients, theta0=760 K. "
            if script_variant
            else "Printed Table 2 alternative, theta0=768 K. "
        )
        + "Conditional numerical parity within reported error widths; exact bootstrap reproduction unverified. See [the source comparison](literature-reproductions/forsterite-katsura-2022-refit-comparison.md).",
        "reproduction_report": "docs/data/forsterite-katsura-2022-reproduction.json",
    }


def update_global_ledgers(report: dict) -> None:
    """Upsert only these two records, preserving unrelated ledger results."""
    from scripts.generate_paper_investigation_ledger import render as render_papers
    from scripts.validate_primary_eos_refits import render_markdown

    document = json.loads(
        (ROOT / "peritheos/data/materials/forsterite.eosmat").read_text(
            encoding="utf-8"
        )
    )
    records = [
        r for r in document["eos_records"] if r["identifier"] in {RECORD, TABLE_RECORD}
    ]
    audit_path = ROOT / "peritheos/data/primary-source-audit.json"
    refit_path = ROOT / "docs/data/primary-eos-refits.json"
    audit = json.loads(audit_path.read_text(encoding="utf-8"))
    refits = json.loads(refit_path.read_text(encoding="utf-8"))
    for record in records:
        validation = record["scientific_validation"]
        outcome = ledger_outcome(record, report)
        base = {
            "material_identifier": "forsterite",
            "record_identifier": record["identifier"],
            "label": record["label"],
            "model": record["eos"]["type"],
            "primary_data_status": validation["primary_data_check"]["status"],
            "primary_source": validation["primary_source_check"]["access_url"],
            "source_notes": record.get("notes"),
            "primary_data_finding": validation["primary_data_check"]["finding"],
            "experimental_pressure_range_gpa": record.get(
                "experimental_pressure_range_gpa"
            ),
            "experimental_temperature_range_k": record.get(
                "experimental_temperature_range_k"
            ),
        }
        entry = {
            "material": "forsterite",
            "file": "forsterite.eosmat",
            "record": record["identifier"],
            "label": record["label"],
            "doi": "10.1029/2021jb023562",
            "status": validation["status"],
            "note": validation["note"],
            "primary_source_check": validation["primary_source_check"],
            "primary_data_check": validation["primary_data_check"],
            "reproduction_status": outcome["status"],
            "parity_basis": outcome["parity_basis"],
            "recommended_thermal_record": RECORD,
            "recommended": outcome["recommended"],
            "reproduction_qualification": outcome["qualification"],
        }
        for ledger, key, value in [
            (audit, "record", entry),
            (refits, "record_identifier", {**base, **outcome}),
        ]:
            old_index = next(
                (
                    i
                    for i, row in enumerate(ledger["records"])
                    if row[key] == record["identifier"]
                ),
                None,
            )
            if old_index is None:
                static_index = next(
                    i
                    for i, row in enumerate(ledger["records"])
                    if row[key] == "forsterite_finkelstein_2014_bm3_1"
                )
                ledger["records"].insert(
                    static_index + 1 + records.index(record), value
                )
            else:
                ledger["records"][old_index] = {**ledger["records"][old_index], **value}
    for ledger, total in [(audit, "records"), (refits, "total")]:
        counts = Counter(row["status"] for row in ledger["records"])
        ledger["summary"] = {
            total: len(ledger["records"]),
            **dict(sorted(counts.items())),
        }
    if LEDGER_PARITY_POLICY not in refits["policy"]["parity"]:
        refits["policy"]["parity"] += " " + LEDGER_PARITY_POLICY
    audit_path.write_text(
        json.dumps(audit, indent=2, ensure_ascii=False, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    refit_path.write_text(
        json.dumps(refits, indent=2, sort_keys=True, allow_nan=False) + "\n",
        encoding="utf-8",
    )
    (ROOT / "docs/primary-eos-refits.md").write_text(
        render_markdown(refits), encoding="utf-8"
    )
    (ROOT / "docs/paper-investigation-ledger.md").write_text(
        render_papers(), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-report", action="store_true")
    parser.add_argument(
        "--update-ledgers",
        action="store_true",
        help="Refresh the two forsterite entries and render the global EOS/paper ledgers",
    )
    args = parser.parse_args()
    result = diagnostics()
    serialized = json.dumps(result, indent=2) + "\n"
    if args.write_report:
        path = ROOT / "docs/data/forsterite-katsura-2022-reproduction.json"
        path.write_text(serialized, encoding="utf-8")
    if args.update_ledgers:
        update_global_ledgers(result)
    print(serialized, end="")


if __name__ == "__main__":
    main()
