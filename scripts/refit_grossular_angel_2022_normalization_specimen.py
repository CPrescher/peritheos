"""Declared normalization/specimen factorial diagnostic, never a catalog refit.

Repository-only: python -m scripts.refit_grossular_angel_2022_normalization_specimen.
Twenty fixed cases isolate central-value changes from previous error sensitivity.
"""

from __future__ import annotations

import hashlib
import json
import math
from collections import Counter

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, observations
from scripts.audit_pyrope_angel_2022_observations import summarize
from scripts.refit_grossular_angel_2022_expanded import native_comparisons
from scripts.refit_grossular_angel_2022_full_sensitivity import (
    prepared_rows as full_rows,
)
from scripts.refit_grossular_angel_2022_full_sensitivity import scenarios, write_csv
from scripts.refit_grossular_angel_2022_subset import solve

PLAN = DATA / "normalization-specimen-plan.json"
BASELINE = "plotted__as_plotted"
DIFFRACTION = [
    "milani_vt",
    "du_vt",
    "bosenick_vt",
    "rodehorst_vt",
    "thieblot_vt_estimated",
]


def cases():
    plan = json.loads(PLAN.read_text())
    return [
        dict(id=f"{n['id']}__{s['id']}", normalization=n["id"], isaak=s["id"])
        for n in plan["normalizations"]
        for s in plan["isaak_expansion_treatments"]
    ]


def line_at(t, t1, v1, t2, v2):
    return v1 + (t - t1) * (v2 - v1) / (t2 - t1)


def denominators(mode):
    original = json.loads((DATA / "normalization-recovery.json").read_text())["groups"]
    values = {
        "milani_vt": original["milani"]["inferred_normalization_denominator"],
        "du_vt": original["du"]["inferred_normalization_denominator"],
        "bosenick_vt": 125.293,
    }
    if mode == "measured_room":
        values.update(milani_vt=1663.98, du_vt=1664.2)
    elif mode == "linear_298":
        values.update(
            milani_vt=line_at(298, 293, 1663.98, 301, 1664.01),
            du_vt=line_at(298, 296.1, 1664.2, 393.7, 1666.0),
            bosenick_vt=line_at(298, 270, 125.226, 295, 125.293),
        )
    elif mode not in {"plotted", "free_milani_du", "free_diffraction"}:
        raise ValueError(f"Unknown normalization: {mode}")
    return values


def isaak_curves():
    original = observations("isaak-1992-table-2.csv")
    first = original[0]
    ratio300 = 1 + 0.01 * float(first["observed_volume_strain_x1e2_token"])
    ratio298 = line_at(298, 293, 1, 300, ratio300)
    temperatures = [float(r["temperature_k"]) for r in original]
    alpha = [
        1e-6 * float(r["estimated_endmember_alpha_x1e6_per_k_token"]) for r in original
    ]
    previous_t = 298
    previous_alpha = line_at(298, temperatures[0], alpha[0], temperatures[1], alpha[1])
    integral = 0
    curves = []
    for r, temperature, a in zip(original, temperatures, alpha):
        integral += 0.5 * (a + previous_alpha) * (temperature - previous_t)
        plotted = 125.35 * (1 + 0.01 * float(r["observed_volume_strain_x1e2_token"]))
        curves.append(
            dict(
                source_row=int(r["source_row"]),
                temperature_k=temperature,
                observed_gr2b_as_plotted=plotted,
                observed_gr2b_reference_298=plotted / ratio298,
                calculated_endmember_from_alpha=125.35 * math.exp(integral),
                integrated_endmember_alpha=integral,
                gr2b_ratio_298_to_293=ratio298,
            )
        )
        previous_t, previous_alpha = temperature, a
    return curves


def prepared_case(case):
    plan = json.loads(PLAN.read_text())
    rows = full_rows(scenarios()[0])
    curves = {r["source_row"]: r for r in isaak_curves()}
    old, new = denominators("plotted"), denominators(case["normalization"])
    for r in rows:
        if r["dataset"] in new:
            factor = old[r["dataset"]] / new[r["dataset"]]
            r["Y"] *= factor
            r["sY"] *= factor
        if r["dataset"] == "isaak_vt":
            curve = curves[r["row"]]
            if case["isaak"] == "reference_298":
                r["Y"] = curve["observed_gr2b_reference_298"]
                r["sY"] /= curve["gr2b_ratio_298_to_293"]
            elif case["isaak"] == "calculated_endmember_alpha":
                r["Y"] = curve["calculated_endmember_from_alpha"]
            elif case["isaak"] not in {"as_plotted", "omit_expansion"}:
                raise ValueError(f"Unknown Isaak treatment: {case['isaak']}")
    if case["isaak"] == "omit_expansion":
        rows = [r for r in rows if r["dataset"] != "isaak_vt"]
    names = list(plan["base_free_parameters"])
    starts = [list(s) for s in plan["starts"]]
    scaling = list(plan["base_parameter_scaling"])
    lower, upper = (list(plan["base_bounds"][k]) for k in ["lower", "upper"])
    free = (
        DIFFRACTION
        if case["normalization"] == "free_diffraction"
        else (DIFFRACTION[:2] if case["normalization"] == "free_milani_du" else [])
    )
    for dataset in free:
        index = len(names)
        names.append(f"scale_{dataset}")
        scaling.append(1)
        lower.append(0.98)
        upper.append(1.02)
        for start in starts:
            start.append(1)
        for r in rows:
            if r["dataset"] == dataset:
                r["scale"] = index
    config = dict(
        names=names,
        scaling=np.array(scaling),
        lower=np.array(lower),
        upper=np.array(upper),
    )
    return rows, starts, config


def normalization_offsets(results):
    """Equivalent raw denominators for fitted calculated-volume multipliers."""
    original = denominators("plotted") | {"rodehorst_vt": 11.850**3}
    resolution = json.loads((DATA / "normalization-recovery.json").read_text())[
        "groups"
    ]
    output = []
    for result in results:
        for dataset in DIFFRACTION:
            name = f"scale_{dataset}"
            parameters = result["fits"][0]["parameters"]
            if name not in parameters:
                continue
            multiplier = parameters[name]
            denominator = original.get(dataset)
            group = dataset.removesuffix("_vt")
            output.append(
                dict(
                    case=result["id"],
                    dataset=dataset,
                    calculated_volume_multiplier=multiplier,
                    percentage_offset=100 * (multiplier - 1),
                    equivalent_observed_offset_at_vm0_cm3_mol=125.35 * (multiplier - 1),
                    original_denominator=denominator,
                    equivalent_denominator=None
                    if denominator is None
                    else denominator * multiplier,
                    denominator_change=None
                    if denominator is None
                    else denominator * (multiplier - 1),
                    denominator_unit="cm3/mol"
                    if dataset == "bosenick_vt"
                    else ("A3" if denominator is not None else "unrecovered"),
                    plot_denominator_resolution=resolution[group][
                        "denominator_resolution_estimate"
                    ]
                    if group in resolution
                    else None,
                )
            )
    return output


def run():
    plan = json.loads(PLAN.read_text())
    all_results, summaries, residuals, native_rows = [], [], [], []
    write_csv("isaak-specimen-treatment-curves.csv", isaak_curves())
    for case in cases():
        rows, starts, config = prepared_case(case)
        fits = [solve(rows, start, **config) for start in starts]
        native = native_comparisons(rows, fits[0])
        for r in native:
            kind = (
                "calculated_from_estimated_endmember_alpha"
                if (
                    r["dataset"] == "isaak_vt"
                    and case["isaak"] == "calculated_endmember_alpha"
                )
                else "processed_input_observation"
            )
            native_rows.append(dict(case=case["id"], input_kind=kind) | r)
        residuals.extend(dict(case=case["id"]) | r for r in fits[0]["rows"])
        groups = {
            name: dict(
                unit=next(r["unit"] for r in native if r["dataset"] == name),
                source=summarize(
                    [r["observed_minus_source"] for r in native if r["dataset"] == name]
                ),
                fit=summarize(
                    [r["observed_minus_fit"] for r in native if r["dataset"] == name]
                ),
            )
            for name in Counter(r["dataset"] for r in native)
        }
        for fit in fits:
            del fit["rows"]
        a, b = (np.array(list(f["parameters"].values())) for f in fits)
        delta = float(np.max(np.abs((a - b) / config["scaling"])))
        result = case | dict(
            rows=len(rows),
            free_parameters=config["names"],
            normalization_denominators=denominators(case["normalization"]),
            fits=fits,
            maximum_scaled_start_difference=delta,
            native_residuals=groups,
        )
        all_results.append(result)
        summaries.append(
            case
            | dict(rows=len(rows), free_parameters=len(config["names"]))
            | {
                name: fits[0]["parameters"][name]
                for name in plan["base_free_parameters"]
            }
            | dict(
                conditional_reduced_objective=fits[0]["conditional_reduced_chi_square"],
                both_starts_converged=all(f["numerically_converged"] for f in fits),
                either_start_at_bound=any(bool(f["at_parameter_bound"]) for f in fits),
            )
        )
        print(
            f"{case['id']}: converged={summaries[-1]['both_starts_converged']}, start difference={delta:.3g}",
            flush=True,
        )
    baseline = all_results[0]["fits"][0]["parameters"]
    p, _ = source("grossular")
    published = dict(
        zip(plan["base_free_parameters"], [p[2], p[3], p[11], p[18], p[53], p[54]])
    )
    by_id = {r["id"]: r["fits"][0]["parameters"] for r in all_results}
    for result in all_results:
        parameters = result["fits"][0]["parameters"]
        norm_only = by_id[f"{result['normalization']}__as_plotted"]
        specimen_only = by_id[f"plotted__{result['isaak']}"]
        result["effects"] = {
            name: dict(
                published=value,
                minus_source=parameters[name] - value,
                minus_baseline=parameters[name] - baseline[name],
                normalization_effect_at_this_specimen=parameters[name]
                - specimen_only[name],
                specimen_effect_at_this_normalization=parameters[name]
                - norm_only[name],
                two_factor_interaction=parameters[name]
                - norm_only[name]
                - specimen_only[name]
                + baseline[name],
            )
            for name, value in published.items()
        }
    report = dict(
        audit_date="2026-11-12",
        scope="user_requested_normalization_specimen_factorial_diagnostic",
        scientific_cases=len(all_results),
        starts_per_case=2,
        baseline_case=BASELINE,
        author_workflow_reproduced=False,
        physical_validation=False,
        source_coefficients_changed=False,
        calculated_endmember_curves_are_measurements=False,
        plan=PLAN.name,
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        all_starts_converged=all(
            f["numerically_converged"] for r in all_results for f in r["fits"]
        ),
        any_start_at_bound=any(
            bool(f["at_parameter_bound"]) for r in all_results for f in r["fits"]
        ),
        maximum_scaled_start_difference=max(
            r["maximum_scaled_start_difference"] for r in all_results
        ),
        cases=all_results,
        comparison_scope="Paired effects at fixed estimated-error baseline. Cases with omitted GR-2B expansion have204rows; other cases have219. Free normalization cases add2or5nuisance parameters. Calculated endmember-alpha curves keep assigned comparison weights solely for the declared central-reduction test; unknown correction covariance is not propagated. Source and fit native predictions use identical per-case fitted nuisance multipliers.",
        input_sha256={
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                "full-sensitivity-refits.json",
                "full-sensitivity-plan.json",
                "normalization-recovery.json",
                "isaak-1992-table-2.csv",
                "isaak-1992-table-3b.csv",
                "figure4-thieblot-selected.csv",
                "milani-2017-table-2.csv",
                "milani-2017-table-3.csv",
                "du-2015-preprint-table-2.csv",
                "bosenick-geiger-1997-table-2.csv",
                "rodehorst-2002-table-4a.csv",
                "kono-2010-table-1.csv",
                "gwanmesia-2014-table-1.csv",
                "bass-1989-ambient.csv",
                "conrad-1999-table-6.csv",
            ]
        },
    )
    write_csv("normalization-specimen-summary.csv", summaries)
    write_csv("normalization-specimen-residuals.csv", residuals)
    write_csv("normalization-specimen-native-comparisons.csv", native_rows)
    write_csv("normalization-specimen-offsets.csv", normalization_offsets(all_results))
    (DATA / "normalization-specimen-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
