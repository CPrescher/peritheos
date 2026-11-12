"""Declared 193-row source-processing diagnostic, not a complete source replay.

Repository-only: python -m scripts.refit_grossular_angel_2022_normalized.
Weights, missing inputs and mixed-composition qualifications are predeclared in
normalized-refit-plan.json. Never writes a catalog record or changes source data.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import numpy as np

from scripts.audit_garnet_angel_2022 import source
from scripts.audit_grossular_angel_2022_observations import DATA, observations
from scripts.audit_pyrope_angel_2022_observations import summarize
from scripts.refit_grossular_angel_2022_expanded import native_comparisons
from scripts.refit_grossular_angel_2022_expanded import prepared_rows as earlier_rows
from scripts.refit_grossular_angel_2022_subset import solve

PLAN = DATA / "normalized-refit-plan.json"


def conrad_width_interval(s11, s12):
    """Admissible KS width bounds for any covariance of the two elastic moduli."""
    return abs(s11 - 2 * s12) / 3, (s11 + 2 * s12) / 3


def prepared_rows():
    rows = earlier_rows()
    normalization = json.loads((DATA / "normalization-recovery.json").read_text())[
        "groups"
    ]
    for dataset, filename, group in [
        ("milani_vt", "milani-2017-table-2.csv", "milani"),
        ("du_vt", "du-2015-preprint-table-2.csv", "du"),
    ]:
        raw = {int(r["source_row"]): r for r in observations(filename)}
        denominator = normalization[group]["inferred_normalization_denominator"]
        for r in rows:
            if r["dataset"] == dataset:
                r["Y"] = 125.35 * float(raw[r["row"]]["cell_volume_a3"]) / denominator
                r["sY"] = (
                    125.35 * float(raw[r["row"]]["cell_volume_sd_a3"]) / denominator
                )
                r["scale"] = None
    for r in rows:
        if r["KS"] and r["scale"] is not None:
            r["scale"] -= 2  # Only the two modulus nuisance factors remain free.

    def add(dataset, r, y, sy, *, t=None, st=None, ks=False, sp=None):
        rows.append(
            dict(
                dataset=dataset,
                row=int(r["source_row"]),
                P=float(r.get("pressure_gpa") or 0.0001),
                T=float(r["temperature_k"]) if t is None else t,
                Y=y,
                sY=sy,
                sT=st,
                sP=sp,
                KS=ks,
                scale=None,
            )
        )

    for r in observations("bosenick-geiger-1997-table-2.csv")[1:]:
        t = float(r["temperature_k"])
        add(
            "bosenick_vt",
            r,
            125.35 * float(r["molar_volume_cm3_mol"]) / 125.293,
            125.35 * float(r["molar_volume_sd_cm3_mol"]) / 125.293,
            st=3 if t >= 100 else 5,
        )
    for r in observations("rodehorst-2002-table-4a.csv"):
        if float(r["temperature_k"]) == 298:
            continue
        a, sa = float(r["lattice_a_angstrom"]), float(r["lattice_a_sd_angstrom"])
        add(
            "rodehorst_vt",
            r,
            125.35 * (a / 11.850) ** 3,
            3 * 125.35 * a**2 * sa / 11.850**3,
        )
    for r in observations("isaak-1992-table-2.csv"):
        add(
            "isaak_vt",
            r,
            125.35 * (1 + 0.01 * float(r["observed_volume_strain_x1e2_token"])),
            0.01,
            st=2,
        )
    for r in observations("conrad-1999-table-6.csv"):
        if r["selection_angel_2022"] != "selected":
            continue
        k = (float(r["c11_gpa"]) + 2 * float(r["c12_gpa"])) / 3
        _, bound = conrad_width_interval(float(r["c11_sd_gpa"]), float(r["c12_sd_gpa"]))
        add("conrad_ks", r, k, bound, t=298, sp=0.1, ks=True)
    return rows


def run():
    plan = json.loads(PLAN.read_text())
    rows = prepared_rows()
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
    native = native_comparisons(rows, fits[0])
    p, _ = source("grossular")
    published = dict(
        zip(plan["free_parameters"], [p[2], p[3], p[11], p[18], p[53], p[54]])
    )
    previous = json.loads((DATA / "expanded-refits.json").read_text())["fits"][0][
        "parameters"
    ]
    report = dict(
        audit_date="2026-11-12",
        scope="conditional_193_row_source_processing_diagnostic",
        author_workflow_reproduced=False,
        full_220_channel_fit=False,
        physical_validation=False,
        source_coefficients_changed=False,
        mixed_composition_isaak_expansion_included=True,
        plan=PLAN.name,
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        channels=len(rows),
        dataset_counts=dict(Counter(r["dataset"] for r in rows)),
        fits=fits,
        parameter_comparison={
            name: dict(
                published=value,
                earlier_147_fit=previous[name],
                normalized_193_fit=fits[0]["parameters"][name],
                normalized_minus_published=fits[0]["parameters"][name] - value,
                percent_difference=100 * (fits[0]["parameters"][name] - value) / value,
            )
            for name, value in published.items()
        },
        native_residuals={
            name: dict(
                unit=next(r["unit"] for r in native if r["dataset"] == name),
                source=summarize(
                    [r["observed_minus_source"] for r in native if r["dataset"] == name]
                ),
                fit=summarize(
                    [r["observed_minus_fit"] for r in native if r["dataset"] == name]
                ),
            )
            for name in Counter(r["dataset"] for r in rows)
        },
        native_residual_scope="Observed minus calculated forward V/KS at measured P/T on identical rows. Same fitted KS multipliers for both curves; all V multipliers fixed by declared normalization. Source coefficient effects, not original author workflow. Pressure objective retained separately in fit rows.",
        native_comparison_csv="normalized-native-comparisons.csv",
        input_sha256={
            name: hashlib.sha256((DATA / name).read_bytes()).hexdigest()
            for name in [
                "milani-2017-table-2.csv",
                "milani-2017-table-3.csv",
                "du-2015-preprint-table-2.csv",
                "kono-2010-table-1.csv",
                "isaak-1992-table-3b.csv",
                "bass-1989-ambient.csv",
                "gwanmesia-2014-table-1.csv",
                "bosenick-geiger-1997-table-2.csv",
                "rodehorst-2002-table-4a.csv",
                "isaak-1992-table-2.csv",
                "conrad-1999-table-6.csv",
                "normalization-recovery.json",
                "remaining-channel-audit.json",
                "expanded-refits.json",
            ]
        },
    )
    with (DATA / report["native_comparison_csv"]).open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=native[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(native)
    (DATA / "normalized-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
