"""Reproduce Sakai (2011) NaCl-B2 scales from the primary printed tables."""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from functools import lru_cache
from itertools import product
from pathlib import Path

import numpy as np
from scipy.optimize import brentq, least_squares

from peritheos import get_eos_record, get_material_document

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUTPUT = ROOT / "docs/data/sakai-2011-nacl-b2-reproduction.json"
PREFIX = "nacl_b2_sakai_2011_"
DATASETS = ["nacl_b2_sakai_2011_table1", "nacl_b2_sata_2002_table1_sakai_input"]


def load_rows(filename):
    source = json.loads((DATA / "nacl-b2-sakai-2011-source.json").read_text())
    payload = (DATA / filename).read_bytes()
    assert hashlib.sha256(payload).hexdigest() == source["datasets"][filename]["sha256"]
    with (DATA / filename).open() as stream:
        return list(csv.DictReader(stream))


def reference_pressure(volume, v0, k0, kp, model):
    """Source equations III.C, independent of the catalog implementation."""
    volume = np.asarray(volume)
    if model == "bm3":
        f = ((v0 / volume) ** (2 / 3) - 1) / 2
        return 3 * k0 * f * (1 + 2 * f) ** 2.5 * (1 + 1.5 * (kp - 4) * f)
    x = (volume / v0) ** (1 / 3)
    return 3 * k0 * (1 - x) / x**2 * np.exp(1.5 * (kp - 1) * (1 - x))


def marker_pressure(volume, scale):
    if scale in ("matsui", "fei"):
        # Matsui (2009) Eq.6/Table II, including the cited Fei Vinet reference.
        k, kp = (273, 5.20) if scale == "matsui" else (277, 5.08)
        return reference_pressure(volume, 60.38, k, kp, "vinet")
    rid = (
        "platinum_dorogokupets_oganov_2007_vinet_4"
        if scale == "do"
        else "platinum_holmes_1989_vinet_1"
    )
    return np.asarray(get_eos_record(rid).pressure(volume, 300))


def observations(scale, selection="all_except_8_11"):
    current = load_rows("nacl-b2-sakai-2011-table1.csv")
    older = load_rows("nacl-b2-sata-2002-table1-sakai-input.csv")
    p = [float(r[f"pressure_{scale}_gpa"]) for r in current]
    v = [float(r["volume_a3"]) for r in current]
    e = [float(r["volume_error_a3"]) for r in current]
    chosen = []
    for row in older:
        if row["excluded_by_sakai"] == "1" or selection == "table1_only":
            continue
        if selection == "annealed_only" and row["annealed"] != "1":
            continue
        # Equal-weight mean of the two cubic Pt lattice parameters. The source
        # does not publish its exact peak-combination weighting: test this choice.
        a = (
            np.sqrt(3) * float(row["pt_d111_angstrom"])
            + 2 * float(row["pt_d200_angstrom"])
        ) / 2
        p.append(float(marker_pressure(a**3, scale)))
        v.append(float(row["volume_a3"]))
        e.append(float(row["volume_error_a3"]))
        chosen.append(int(row["source_row"]))
    return np.array(p), np.array(v), np.array(e), chosen


def fit(p, v, v0, model, errors=None):
    def residual(parameters):
        pressure = reference_pressure(v, v0, *parameters, model)
        if errors is None:
            return pressure - p
        # First-order volume-error propagation; no missing pressure sigma is
        # invented. This is only an alternative diagnostic objective.
        step = 1e-5
        slope = (
            reference_pressure(v + step, v0, *parameters, model)
            - reference_pressure(v - step, v0, *parameters, model)
        ) / (2 * step)
        return (pressure - p) / (np.abs(slope) * errors)

    result = least_squares(
        residual,
        [45, 4.5],
        bounds=([1, 1], [200, 12]),
        xtol=1e-12,
        ftol=1e-12,
        gtol=1e-12,
    )
    assert result.success
    pred = reference_pressure(v, v0, *result.x, model)
    return {
        "K0": float(result.x[0]),
        "K0_prime": float(result.x[1]),
        "rmse_gpa": float(np.sqrt(np.mean((pred - p) ** 2))),
        "observations": len(p),
    }


def staged_volume(p, v):
    # Figure 4 uses the conventional B1 cell 179.42 A^3 as numerical V01.
    # G=P/[3(1+2g)^(+5/2)] is consistent with Figure 4 and BM3; the
    # negative exponent printed inside the denominator is not.
    g = ((v / 179.42) ** (-2 / 3) - 1) / 2
    normalized = p / (3 * (1 + 2 * g) ** 2.5)
    coefficients = np.polyfit(g, normalized, 2)
    roots = np.roots(coefficients)
    candidates = [
        float(179.42 / (1 + 2 * r.real) ** 1.5)
        for r in roots
        if abs(r.imag) < 1e-10 and r.real > -0.5
    ]
    physical = [value for value in candidates if max(v) < value < 179.42]
    if len(physical) != 1:
        raise ValueError(f"Ambiguous g-G extrapolation: {candidates}")
    return {
        "V0": physical[0],
        "quadratic_coefficients": coefficients.tolist(),
        "all_real_candidate_volumes_a3": candidates,
    }


@lru_cache(maxsize=1)
def reproduce():
    result = {"doi": "10.1063/1.3573393", "records": {}, "calibrant_checks": {}}
    rows = load_rows("nacl-b2-sakai-2011-table1.csv")
    for scale in ("matsui", "fei", "do", "holmes"):
        observed = np.array([float(r[f"pressure_{scale}_gpa"]) for r in rows])
        recalculated = marker_pressure(
            np.array([float(r["pt_volume_a3"]) for r in rows]), scale
        )
        result["calibrant_checks"][scale] = {
            "max_abs_pressure_difference_gpa": float(max(abs(recalculated - observed))),
            "rmse_gpa": float(np.sqrt(np.mean((recalculated - observed) ** 2))),
        }
    for record in get_material_document("nacl_b2")["eos_records"]:
        rid = record["identifier"]
        if not rid.startswith(PREFIX):
            continue
        scale, model = rid.removeprefix(PREFIX).split("_pt_")
        params = record["eos"]["parameters"]
        v0, k, kp = (params[name] for name in ("V0", "K0", "K0_prime"))
        p, v, errors, chosen = observations(scale)
        source_curve = reference_pressure(v, v0, k, kp, model)
        fixed = fit(p, v, v0, model)
        stage = staged_volume(p, v)
        stage["final_fit"] = fit(p, v, stage["V0"], model)
        a_p, a_v, _, a_chosen = observations(scale, "annealed_only")
        t_p, t_v, _, _ = observations(scale, "table1_only")
        benchmark = reference_pressure(t_v, v0, k, kp, model)
        # Independent primary measurements: compatibility bounded by the printed
        # volume errors, using three error widths without asserting confidence.
        t_errors = np.array([float(r["volume_error_a3"]) for r in rows])
        bounds = (
            np.abs(
                reference_pressure(t_v - 3 * t_errors, v0, k, kp, model)
                - reference_pressure(t_v + 3 * t_errors, v0, k, kp, model)
            )
            / 2
        )
        rounded_tolerance = np.maximum(bounds, 0.1)
        result["records"][rid] = {
            "fixed_V0_refit": fixed,
            "staged_gG_refit": stage,
            "annealed_only_refit": {**fit(a_p, a_v, v0, model), "sata_rows": a_chosen},
            "table1_only_refit": fit(t_p, t_v, v0, model),
            "volume_error_weighted_refit": fit(p, v, v0, model, errors),
            "selected_sata_rows": chosen,
            "published_rmse_gpa": float(np.sqrt(np.mean((source_curve - p) ** 2))),
            "table1_published_rmse_gpa": float(
                np.sqrt(np.mean((benchmark - t_p) ** 2))
            ),
            "table1_max_residual_gpa": float(max(abs(benchmark - t_p))),
            "table1_points_within_three_printed_volume_error_widths": int(
                np.sum(abs(benchmark - t_p) <= rounded_tolerance)
            ),
            "catalog_max_difference_from_independent_equation_gpa": float(
                max(abs(np.asarray(get_eos_record(rid).pressure(v)) - source_curve))
            ),
            "measured_table1_pressure_range_gpa": [float(min(t_p)), float(max(t_p))],
            "reconstructed_joint_pressure_range_gpa": [float(min(p)), float(max(p))],
            "coefficient_differences": {
                "K0": fixed["K0"] - k,
                "K0_prime": fixed["K0_prime"] - kp,
            },
        }
    base = "nacl_b2_sakai_2011_matsui_pt_bm3"
    reference_volume = float(get_eos_record(base).volume(364))
    differences = {}
    for scale in ("matsui", "fei", "do", "holmes"):
        record = get_eos_record(f"{PREFIX}{scale}_pt_bm3")
        differences[scale] = float(record.pressure(reference_volume) - 364)
    rounding_bounds = {}
    # Each Table II BM3 coefficient is printed to two decimal places.
    # Vary both reference and target curves over half-last-digit intervals.
    base_parameters = (37.73, 47.00, 4.10)
    for scale in ("fei", "do", "holmes"):
        parameters = get_eos_record(f"{PREFIX}{scale}_pt_bm3")
        raw = next(
            r
            for r in get_material_document("nacl_b2")["eos_records"]
            if r["identifier"] == parameters.identifier
        )["eos"]["parameters"]
        target = tuple(raw[key] for key in ("V0", "K0", "K0_prime"))
        outcomes = []
        for offsets in product((-0.005, 0.005), repeat=6):
            b = tuple(
                value + offset for value, offset in zip(base_parameters, offsets[:3])
            )
            t = tuple(value + offset for value, offset in zip(target, offsets[3:]))
            volume = brentq(lambda v: reference_pressure(v, *b, "bm3") - 364, 10, 20)
            outcomes.append(float(reference_pressure(volume, *t, "bm3") - 364))
        rounding_bounds[scale] = [min(outcomes), max(outcomes)]
    result["independent_published_comparisons"] = {
        "source_location": "p.084912-5 and Figure 5; extrapolation beyond measurements",
        "matsui_bm3_pressure_gpa": 364,
        "volume_a3": reference_volume,
        "bm3_scale_difference_gpa": differences,
        "source_reported_difference_gpa": {"fei": -4, "do": 12, "holmes": 37},
        "half_last_digit_coefficient_rounding_bounds_gpa": rounding_bounds,
        "matsui_vinet_minus_bm3_gpa": float(
            get_eos_record(f"{PREFIX}matsui_pt_vinet").pressure(reference_volume) - 364
        ),
        "reported_matsui_model_difference_gpa": 5,
    }
    return result


def ledger_outcome(record):
    result = reproduce()["records"][record["identifier"]]
    parameters = []
    for name in ("K0", "K0_prime"):
        published = record["eos"]["parameters"][name]
        value = result["fixed_V0_refit"][name]
        error = record["parameter_errors"][name]
        parameters.append(
            {
                "parameter": name,
                "published": published,
                "published_error": error,
                "refit": value,
                "refit_error": None,
                "difference": value - published,
                "relative_difference": abs(value - published) / published,
                "within_combined_2sigma": None,
                "within_reported_error": abs(value - published) <= error,
                "similar": abs(value - published) / published < 0.15,
            }
        )
    return {
        "status": "parity"
        if all(p["within_reported_error"] for p in parameters)
        else "similar"
        if all(p["similar"] for p in parameters)
        else "parity_not_achieved",
        "parity_basis": "conditional_fixed_V0_diagnostic_only",
        "dataset_identifiers": DATASETS,
        "observations": result["fixed_V0_refit"]["observations"],
        "fit_kind": "primary_table_fixed_V0_bm3_or_vinet",
        "objective": "unweighted_pressure_residuals",
        "selection": "27 Sakai rows plus 27 Sata rows after explicit exclusions 8 and 11; annealed-only selection retained as sensitivity",
        "free_parameters": ["K0", "K0_prime"],
        "fixed_parameters": ["V0"],
        "parameters": parameters,
        "rmse_gpa": result["fixed_V0_refit"]["rmse_gpa"],
        "published_rmse_gpa": result["published_rmse_gpa"],
        "qualification": "Conditional diagnostic, not exact source-regression recovery: source weights and Pt peak-combination rule are unspecified. Sata pressures reconstructed from equal-weight cubic Pt lattice parameters. Staged V0 and selection/objective sensitivity documented in sakai-2011-nacl-b2 reproduction; published parameters unchanged.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    result = reproduce()
    for entry in result["records"].values():
        assert entry["catalog_max_difference_from_independent_equation_gpa"] < 1e-10
    comparison = result["independent_published_comparisons"]
    for scale, expected in comparison["source_reported_difference_gpa"].items():
        low, high = comparison["half_last_digit_coefficient_rounding_bounds_gpa"][scale]
        assert low - 0.5 <= expected <= high + 0.5

    # Round diagnostics for readable output in the contributor environment.
    def rounded(value):
        if isinstance(value, float):
            return round(value, 8)
        if isinstance(value, dict):
            return {k: rounded(v) for k, v in value.items()}
        if isinstance(value, list):
            return [rounded(v) for v in value]
        return value

    text = json.dumps(rounded(result), indent=2, sort_keys=True) + "\n"
    if args.check:
        return int(not OUTPUT.exists() or OUTPUT.read_text() != text)
    OUTPUT.write_text(text)
    print(text, end="")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
