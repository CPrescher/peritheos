"""Full-inventory diagnostic with explicitly estimated errors and sensitivity.

Repository-only: python -m scripts.refit_grossular_angel_2022_full_sensitivity.
The predeclared plan fixes 29 cases; no weights are tuned to published values.
Original observations, blank source errors and earlier fits remain unchanged.
"""

from __future__ import annotations

import csv
import hashlib
import itertools
import json
import math

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, observations
from scripts.audit_pyrope_angel_2022_observations import summarize
from scripts.refit_grossular_angel_2022_expanded import native_comparisons
from scripts.refit_grossular_angel_2022_normalized import (
    prepared_rows as reference_rows,
)
from scripts.refit_grossular_angel_2022_subset import solve

PLAN = DATA / "full-sensitivity-plan.json"
BASELINE = "interpolated_i100_v100_t100"
FACTORS = (
    "isaak_multiplier",
    "thieblot_volume_multiplier",
    "thieblot_temperature_multiplier",
)


def scenarios():
    def scenario(i, v, t, profile="interpolated"):
        return dict(
            id=f"{profile}_i{round(100 * i):03d}_v{round(100 * v):03d}_t{round(100 * t):03d}",
            isaak_profile=profile,
            **dict(zip(FACTORS, [i, v, t])),
        )

    cases = [scenario(1, 1, 1)]
    cases.extend(
        scenario(*factors)
        for factors in itertools.product([0.5, 1, 2], repeat=3)
        if factors != (1, 1, 1)
    )
    cases.extend(
        [scenario(1, 1, 1, "constant_low"), scenario(1, 1, 1, "constant_high")]
    )
    return cases


def estimated_isaak_width(temperature, profile):
    if profile == "constant_low":
        return 0.71
    if profile == "constant_high":
        return 1.20
    if profile != "interpolated":
        raise ValueError(f"Unknown estimated-error profile: {profile}")
    if temperature <= 800:
        return 0.71 + (temperature - 300) * (1.02 - 0.71) / 500
    # Includes the explicitly declared extrapolation to 1350 K.
    return 1.02 + (temperature - 800) * (1.20 - 1.02) / 500


def prepared_rows(case):
    rows = reference_rows()
    for r in observations("isaak-1992-table-3b.csv"):
        if r["ks_gpa_reported_error"]:
            continue  # All three printed-error rows are already in the reference.
        temperature = float(r["temperature_k"])
        rows.append(
            dict(
                dataset="isaak_ks_estimated",
                row=int(r["source_row"]),
                P=0.0001,
                T=temperature,
                Y=float(r["ks_gpa"]),
                sY=case["isaak_multiplier"]
                * estimated_isaak_width(temperature, case["isaak_profile"]),
                sP=None,
                sT=None,
                KS=True,
                scale=None,
            )
        )
    for r in observations("figure4-thieblot-selected.csv"):
        temperature = float(r["temperature_k"])
        rows.append(
            dict(
                dataset="thieblot_vt_estimated",
                row=int(r["drawing_index"]),
                P=0.0001,
                T=temperature,
                Y=float(r["molar_volume_cm3_mol"]),
                sY=math.sqrt(
                    (case["thieblot_volume_multiplier"] * 0.03) ** 2 + 0.003**2 / 3
                ),
                sT=math.sqrt(
                    (case["thieblot_temperature_multiplier"] * 0.01 * temperature) ** 2
                    + 2**2 / 3
                ),
                sP=None,
                KS=False,
                scale=None,
            )
        )
    return rows


def write_csv(name, rows):
    with (DATA / name).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=rows[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(rows)


def run():
    plan = json.loads(PLAN.read_text())
    results, residual_rows, summary_rows = [], [], []
    baseline_native = None
    for case in scenarios():
        rows = prepared_rows(case)
        fits = [
            solve(
                rows,
                start,
                names=plan["free_parameters"],
                scaling=np.array(plan["parameter_scaling"]),
                lower=np.array(plan["bounds"]["lower"]),
                upper=np.array(plan["bounds"]["upper"]),
            )
            for start in plan["starts"]
        ]
        if case["id"] == BASELINE:
            baseline_native = native_comparisons(rows, fits[0])
            write_csv("full-baseline-native-comparisons.csv", baseline_native)
            write_csv(
                "full-baseline-inputs.csv",
                [
                    r
                    | dict(
                        added_estimated_errors=r["dataset"].endswith("_estimated"),
                        observation_kind="digitized_processed_observation"
                        if r["dataset"] == "thieblot_vt_estimated"
                        else "tabulated_or_table_derived",
                    )
                    for r in rows
                ],
            )
        for r in fits[0]["rows"]:
            residual_rows.append(dict(case=case["id"]) | r)
        for fit in fits:
            del fit[
                "rows"
            ]  # One row receipt per case; no duplicate second-start table.
        a, b = (np.array(list(f["parameters"].values())) for f in fits)
        agreement = float(np.max(np.abs((a - b) / np.array(plan["parameter_scaling"]))))
        results.append(
            case | dict(fits=fits, maximum_scaled_start_difference=agreement)
        )
        summary_rows.append(
            case
            | fits[0]["parameters"]
            | dict(
                conditional_reduced_objective=fits[0]["conditional_reduced_chi_square"],
                both_starts_converged=all(f["numerically_converged"] for f in fits),
                either_start_at_bound=any(bool(f["at_parameter_bound"]) for f in fits),
            )
        )
        print(
            f"{case['id']}: converged={summary_rows[-1]['both_starts_converged']}, start difference={agreement:.3g}",
            flush=True,
        )
    write_csv("full-sensitivity-residuals.csv", residual_rows)
    write_csv("full-sensitivity-summary.csv", summary_rows)
    names = plan["free_parameters"]
    baseline = results[0]["fits"][0]["parameters"]
    p, _ = source("grossular")
    published = dict(zip(names, [p[2], p[3], p[11], p[18], p[53], p[54]]))
    previous = json.loads((DATA / "normalized-refits.json").read_text())["fits"][0][
        "parameters"
    ]

    def envelope(selected):
        return {
            name: dict(
                minimum=min(r["fits"][0]["parameters"][name] for r in selected),
                maximum=max(r["fits"][0]["parameters"][name] for r in selected),
            )
            for name in names
        }

    one_factor = {}
    for factor in FACTORS:
        selected = [
            r
            for r in results
            if r["isaak_profile"] == "interpolated"
            and all(r[key] == 1 for key in FACTORS if key != factor)
        ]
        one_factor[factor] = dict(
            case_ids=[r["id"] for r in selected], envelope=envelope(selected)
        )
    groups = sorted(set(r["dataset"] for r in baseline_native))
    report = dict(
        audit_date="2026-11-12",
        scope="user_requested_full_inventory_estimated_error_sensitivity",
        fitted_rows=219,
        original_reference_rows=193,
        newly_estimated_error_rows=26,
        nominal_inventory_channels=220,
        independent_normalization_anchor_rows_added=0,
        author_workflow_reproduced=False,
        physical_validation=False,
        original_errors_recovered=False,
        source_coefficients_changed=False,
        plan=PLAN.name,
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        baseline_case=BASELINE,
        scientific_cases=len(results),
        starts_per_case=2,
        all_starts_converged=all(
            f["numerically_converged"] for r in results for f in r["fits"]
        ),
        any_start_at_bound=any(
            bool(f["at_parameter_bound"]) for r in results for f in r["fits"]
        ),
        maximum_scaled_start_difference=max(
            r["maximum_scaled_start_difference"] for r in results
        ),
        comparison={
            name: dict(
                published=value,
                reference_193=previous[name],
                baseline_219=baseline[name],
                baseline_minus_published=baseline[name] - value,
                percent_difference=100 * (baseline[name] - value) / value,
            )
            for name, value in published.items()
        },
        all_case_envelope=envelope(results),
        multiplier_grid_envelope=envelope(results[:27]),
        one_factor_sensitivity=one_factor,
        shape_control_case_ids=[
            r["id"] for r in results if r["isaak_profile"] != "interpolated"
        ],
        baseline_native_residuals={
            name: dict(
                unit=next(r["unit"] for r in baseline_native if r["dataset"] == name),
                source=summarize(
                    [
                        r["observed_minus_source"]
                        for r in baseline_native
                        if r["dataset"] == name
                    ]
                ),
                fit=summarize(
                    [
                        r["observed_minus_fit"]
                        for r in baseline_native
                        if r["dataset"] == name
                    ]
                ),
            )
            for name in groups
        },
        native_comparison_scope="Observed minus forward V/KS at measured P/T on identical baseline rows. Both curves share fitted modulus multipliers. All volume normalizations fixed as in 193-row reference; pressure-coordinate objectives are retained separately. Ranges are weight sensitivity, not confidence intervals.",
        cases=results,
        input_sha256={
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                "normalized-refits.json",
                "normalized-refit-plan.json",
                "normalization-recovery.json",
                "isaak-1992-table-3b.csv",
                "figure4-thieblot-selected.csv",
                "figure4-recovery.json",
                "milani-2017-table-2.csv",
                "milani-2017-table-3.csv",
                "du-2015-preprint-table-2.csv",
                "kono-2010-table-1.csv",
                "bass-1989-ambient.csv",
                "gwanmesia-2014-table-1.csv",
                "bosenick-geiger-1997-table-2.csv",
                "rodehorst-2002-table-4a.csv",
                "isaak-1992-table-2.csv",
                "conrad-1999-table-6.csv",
            ]
        },
    )
    (DATA / "full-sensitivity-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
