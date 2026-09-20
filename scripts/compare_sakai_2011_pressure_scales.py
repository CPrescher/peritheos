"""Compare Pt-scale conversion with the separately fitted Sakai (2011) curves."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from peritheos import (
    get_eos_record,
    get_eos_record_document,
    recalculate_eos_pressure_scale,
)
from scripts.reproduce_sakai_2011_nacl_b2 import PREFIX, load_rows, marker_pressure

ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "docs/data/sakai-2011-pressure-scale-comparison.json"
STANDARDS = {
    "matsui": "platinum_matsui_2009_vinet_300k",
    "fei": "platinum_fei_2007_vinet_300k",
    "do": "platinum_dorogokupets_oganov_2007_vinet_4",
    "holmes": "platinum_holmes_1989_vinet_1",
}
RESOLVED = ("matsui", "fei", "do")


def metrics(volume, converted, published):
    difference = converted - published
    largest = int(np.argmax(np.abs(difference)))
    return {
        "sampled_states": int(volume.size),
        "sample_volume_range_a3": [float(volume.min()), float(volume.max())],
        "published_target_pressure_range_gpa": [
            float(published.min()),
            float(published.max()),
        ],
        "converted_pressure_range_gpa": [
            float(converted.min()),
            float(converted.max()),
        ],
        "max_abs_difference_gpa": float(np.max(np.abs(difference))),
        "max_abs_relative_difference_percent": float(
            np.max(np.abs(difference / published)) * 100
        ),
        "signed_difference_range_gpa": [
            float(difference.min()),
            float(difference.max()),
        ],
        "rms_difference_gpa": float(np.sqrt(np.mean(difference**2))),
        "largest_difference_at": {
            "sample_volume_a3": float(volume[largest]),
            "published_target_pressure_gpa": float(published[largest]),
            "converted_pressure_gpa": float(converted[largest]),
        },
    }


def compare_pair(source_scale, target_scale, model, points=1001):
    source_id = PREFIX + source_scale + "_pt_" + model
    target_id = PREFIX + target_scale + "_pt_" + model
    source = get_eos_record(source_id)
    target = get_eos_record(target_id)
    source_standard = get_eos_record(STANDARDS[source_scale])
    target_standard = get_eos_record(STANDARDS[target_scale])
    # Intersect the two sample-fit envelopes in NaCl volume, not in pressure:
    # each pressure scale assigns a different pressure to the same state.
    source_range = get_eos_record_document(source_id)["validity"]["pressure_gpa"]
    target_range = get_eos_record_document(target_id)["validity"]["pressure_gpa"]
    low = max(source.volume(source_range[1]), target.volume(target_range[1]))
    high = min(source.volume(source_range[0]), target.volume(target_range[0]))
    volume = np.linspace(low, high, points)
    converted_result = recalculate_eos_pressure_scale(
        source_id, volume, STANDARDS[target_scale], 300.0
    )
    converted = np.asarray(converted_result.target_pressure_gpa)
    published = np.asarray(target.pressure(volume, 300.0))
    pt_volume = source_standard.volume(source.pressure(volume, 300.0), 300.0)
    strict = (
        source.within_validity(volume, 300.0)
        & target.within_validity(volume, 300.0)
        & source_standard.within_validity(pt_volume, 300.0)
        & target_standard.within_validity(pt_volume, 300.0)
    )
    result = {
        "source_eos_record": source_id,
        "published_target_eos_record": target_id,
        "target_standard_eos_record": STANDARDS[target_scale],
        "calibration_path": list(converted_result.calibration_path),
        "comparison_status": (
            "exact_calibration_links"
            if target_scale in RESOLVED
            else "diagnostic_only_unresolved_holmes_variant"
        ),
        "sample_envelope": metrics(volume, converted, published),
        "all_linked_envelopes": None,
        "states_outside_linked_envelopes": int(np.count_nonzero(~strict)),
    }
    if np.any(strict):
        checked = recalculate_eos_pressure_scale(
            source_id,
            volume[strict],
            STANDARDS[target_scale],
            300.0,
            check_validity=True,
        )
        np.testing.assert_allclose(
            checked.target_pressure_gpa, converted[strict], atol=1e-10, rtol=0
        )
        result["all_linked_envelopes"] = metrics(
            volume[strict], converted[strict], published[strict]
        )
    return result


def reproduce(points=1001):
    comparisons = {}
    for model in ("bm3", "vinet"):
        # Full directed matrix for the three resolved calibrations, including
        # identity controls. Matsui-to-Holmes is an explicitly qualified diagnostic.
        for source in RESOLVED:
            for target in RESOLVED:
                comparisons[f"{source}_to_{target}_{model}"] = compare_pair(
                    source, target, model, points
                )
        comparisons[f"matsui_to_holmes_{model}"] = compare_pair(
            "matsui", "holmes", model, points
        )
    rows = load_rows("nacl-b2-sakai-2011-table1.csv")
    pt_volume = np.array([float(r["pt_volume_a3"]) for r in rows])
    checks = {}
    for scale in ("matsui", "fei"):
        standard = get_eos_record(STANDARDS[scale])
        prediction = np.asarray(standard.pressure(pt_volume, 300.0))
        primary = np.array([float(r[f"pressure_{scale}_gpa"]) for r in rows])
        checks[scale] = {
            "max_difference_from_independent_vinet_gpa": float(
                np.max(np.abs(prediction - marker_pressure(pt_volume, scale)))
            ),
            "max_difference_from_sakai_marker_pressures_gpa": float(
                np.max(np.abs(prediction - primary))
            ),
        }
    checks["matsui"]["table3_300k_0p7V0_pressure_gpa"] = float(
        get_eos_record(STANDARDS["matsui"]).pressure(0.7 * 60.38, 300.0)
    )
    return {
        "method": "Convert each published NaCl curve through a virtual Pt volume at 300 K; compare with the same-form published target curve at identical NaCl volume. No refit or coefficient replacement.",
        "sampling": f"{points} equally spaced NaCl volumes per pair across the intersection of both published sample envelopes; maxima are sampled, not analytic bounds.",
        "qualifications": [
            "Matsui is the display baseline, not a new recommended material default.",
            "The all_linked_envelopes subset passes the public API with check_validity=True; full sample-envelope comparisons explicitly include calibrant extrapolation.",
            "Holmes uses the exact bundled Equation (11) as a diagnostic. Sakai's precise Holmes variant is unresolved, so its published NaCl records remain unlinked and cannot be used as an automatically resolved source calibration.",
            "These are deterministic differences between rounded published curves, not experimental errors or confidence bounds. The variants share observations and are not independent experiments.",
        ],
        "calibrant_checks": checks,
        "comparisons": comparisons,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()

    def rounded(value):
        if isinstance(value, float):
            return round(value, 8)
        if isinstance(value, dict):
            return {k: rounded(v) for k, v in value.items()}
        if isinstance(value, list):
            return [rounded(v) for v in value]
        return value

    text = json.dumps(rounded(reproduce()), indent=2, sort_keys=True) + "\n"
    if args.check:
        return int(not OUTPUT.exists() or OUTPUT.read_text() != text)
    OUTPUT.write_text(text)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
