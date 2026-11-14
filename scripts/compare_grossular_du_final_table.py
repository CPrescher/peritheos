"""Compare the separately recovered final grossular Table 2 with the preprint.

Repository-only: python -m scripts.compare_grossular_du_final_table.
Reads factual transcriptions, never fits or changes source observations.
PDF versions and visual source checks are recorded in the evidence JSON.
"""

from __future__ import annotations

import hashlib
import json

from scripts.audit_grossular_angel_2022_observations import DATA, observations

PREPRINT = "du-2015-preprint-table-2.csv"
FINAL = "du-2015-final-table-2.csv"
EVIDENCE = "du-final-table-source-evidence.json"
OUTPUT = "du-final-table-comparison.json"
NUMERIC = [
    "temperature_k",
    "temperature_error_k",
    "lattice_a_angstrom",
    "lattice_a_sd_angstrom",
    "cell_volume_a3",
    "cell_volume_sd_a3",
]


def compare(before, after, alpha):
    """Keep row identity and all original CSV fields, including printed tokens."""
    if len(before) != len(after) or len(before) != len(alpha):
        raise ValueError("Source row counts differ")
    differences = []
    for old, new, calculated in zip(before, after, alpha):
        fields = [key for key in old if old[key] != new[key]]
        if new["calculated_alpha_x1e5_per_k_token"] != calculated:
            fields.append("calculated_alpha_x1e5_per_k_token")
        if fields:
            differences.append(dict(source_row=old["source_row"], fields=fields))
    return dict(
        rows=len(before),
        differing_rows=differences,
        identical=not differences,
        maximum_absolute_differences={
            key: max(abs(float(a[key]) - float(b[key])) for a, b in zip(before, after))
            for key in NUMERIC
        },
        shared_fit_fields=list(before[0]),
        calculated_alpha_compared=True,
    )


def report():
    evidence = json.loads((DATA / EVIDENCE).read_text())
    return dict(
        checked_on="2026-11-13",
        scope="Grossular Table2 and associated definitions only, not whole-paper identity",
        files={
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [PREPRINT, FINAL, EVIDENCE]
        },
        comparison=compare(
            observations(PREPRINT),
            observations(FINAL),
            evidence["preprint_alpha_tokens_visually_transcribed"],
        ),
        conclusion="Final publication adds no changed grossular Table2 fit input; version changes in this block cannot explain the remaining theta discrepancy",
        fits_rerun=False,
        author_workflow_reproduced=False,
        author_contacted=False,
        conrad_2percent="Retained diagnostic, no source-error change",
    )


if __name__ == "__main__":
    result = report()
    (DATA / OUTPUT).write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result["comparison"], indent=2))
