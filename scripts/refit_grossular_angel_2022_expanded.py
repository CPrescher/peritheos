"""Predeclared 147-channel diagnostic with the recovered Gwanmesia moduli.

Repository-only: python -m scripts.refit_grossular_angel_2022_expanded.
Never writes catalog records. See grossular/expanded-refit-plan.json.
"""

from __future__ import annotations

import csv
import hashlib
import json
from collections import Counter

import numpy as np

from scripts.audit_garnet_angel_2022 import FOLDER, source
from scripts.audit_grossular_angel_2022_observations import (
    DATA,
    GAS_CONSTANT,
    observations,
)
from scripts.audit_pyrope_angel_2022_observations import properties, summarize
from scripts.refit_grossular_angel_2022_subset import (
    NAMES,
    solve,
)
from scripts.refit_grossular_angel_2022_subset import (
    prepared_rows as previous_rows,
)

PLAN = DATA / "expanded-refit-plan.json"


def prepared_rows():
    rows = previous_rows()
    for i, r in enumerate(observations("gwanmesia-2014-table-1.csv"), 1):
        pressure, ks = float(r["pressure_gpa"]), float(r["ks_gpa"])
        rows.append(
            dict(
                dataset="gwanmesia_ks",
                row=i,
                P=pressure,
                T=float(r["temperature_c"]) + 273.15,
                Y=ks,
                sY=0.008 * ks,
                sP=0.007 * pressure,
                sT=None,
                KS=True,
                scale=7,
            )
        )
    return rows


def native_comparisons(rows, fit):
    """Compare independent forward properties at measured P/T on identical rows.

    Source and fit share fitted nuisance scales to isolate coefficient effects.
    These native-coordinate residuals are not the pressure-coordinate objective.
    """
    p, _ = source("grossular")
    c = list(fit["parameters"].values())
    fitted = p | dict(zip([2, 3, 11, 18], c[:4]))
    result = []
    for r in rows:
        field = "ks_gpa" if r["KS"] else "volume_cm3_mol"
        multiplier = c[r["scale"]] if r["scale"] is not None else 1
        source_y = multiplier * properties(p, r["P"], r["T"], GAS_CONSTANT)[field]
        fit_y = multiplier * properties(fitted, r["P"], r["T"], GAS_CONSTANT)[field]
        result.append(
            dict(
                dataset=r["dataset"],
                source_row=r["row"],
                pressure_gpa=r["P"],
                temperature_k=r["T"],
                coordinate="KS" if r["KS"] else "V",
                unit="GPa" if r["KS"] else "cm3/mol",
                observed=r["Y"],
                source=source_y,
                fit=fit_y,
                observed_minus_source=r["Y"] - source_y,
                observed_minus_fit=r["Y"] - fit_y,
                common_nuisance_multiplier=multiplier,
            )
        )
    return result


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
    p, _ = source("grossular")
    previous = json.loads((DATA / "subset-refits.json").read_text())["fits"][0]
    published = dict(zip(NAMES[:4], [p[2], p[3], p[11], p[18]])) | {
        "scale_kono_ks": p[53],
        "scale_gwanmesia_ks": p[54],
    }
    comparison = {}
    for name, value in published.items():
        fitted = fits[0]["parameters"][name]
        old = previous["parameters"].get(name)
        comparison[name] = dict(
            published=value,
            earlier_95_channel_fit=old,
            expanded_147_channel_fit=fitted,
            expanded_minus_published=fitted - value,
            percent_difference=100 * (fitted - value) / value,
            absolute_distance_decreased=None
            if old is None
            else abs(fitted - value) < abs(old - value),
        )
    native = native_comparisons(rows, fits[0])
    summaries = {}
    for name in Counter(r["dataset"] for r in rows):
        group = [r for r in native if r["dataset"] == name]
        summaries[name] = {
            "unit": group[0]["unit"],
            "source": summarize([r["observed_minus_source"] for r in group]),
            "fit": summarize([r["observed_minus_fit"] for r in group]),
        }
    report = dict(
        audit_date="2026-11-12",
        scope="conditional_147_channel_expanded_refit",
        author_workflow_reproduced=False,
        source_coefficients_changed=False,
        catalog_record_written=False,
        plan=PLAN.name,
        plan_sha256=hashlib.sha256(PLAN.read_bytes()).hexdigest(),
        channels=len(rows),
        dataset_counts=dict(Counter(r["dataset"] for r in rows)),
        fits=fits,
        parameter_comparison=comparison,
        native_residual_comparison={
            "scope": "Observed minus forward-calculated V or KS at the measured P/T; identical fitted nuisance multipliers for source and fitted coefficients, isolating coefficient effects. These are not the minimized pressure residuals, and source nuisance multipliers here are not the original published ones.",
            "groups": summaries,
            "csv": "expanded-native-comparisons.csv",
        },
        gwanmesia_source_original_scale_comparison=json.loads(
            (DATA / "followup-recovery.json").read_text()
        )["gwanmesia"]["scaled_residual_gpa"],
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
                "gwanmesia-2014-source.json",
                "subset-refits.json",
                "followup-recovery.json",
            ]
        },
        official_source_sha256=hashlib.sha256(
            (FOLDER / "sources/grossular-official.eos").read_bytes()
        ).hexdigest(),
    )
    with (DATA / "expanded-native-comparisons.csv").open(
        "w", newline="", encoding="utf-8"
    ) as stream:
        writer = csv.DictWriter(stream, fieldnames=native[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(native)
    (DATA / "expanded-refits.json").write_text(
        json.dumps(report, indent=2, allow_nan=False) + "\n"
    )
    return report


if __name__ == "__main__":
    run()
