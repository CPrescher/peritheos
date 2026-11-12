"""Extract numerical normalization evidence from final Figure 4, not fit inputs.

Repository-only; requires PyMuPDF and a privately obtained final PDF. Original
tables and earlier digitization receipts are never rewritten.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

from scripts.audit_grossular_angel_2022_observations import DATA, observations
from scripts.audit_pyrope_angel_2022_observations import summarize
from scripts.recover_grossular_angel_2022_figure4 import PDF_SHA256


def recover(pdf):
    import fitz

    if hashlib.sha256(pdf.read_bytes()).hexdigest() != PDF_SHA256:
        raise ValueError("Extraction qualified only for the recorded final PDF")
    with fitz.open(pdf) as document:
        drawings = document[7].get_drawings()
    axes = json.loads((DATA / "figure4-recovery.json").read_text())
    tc, vc = axes["temperature_axis"], axes["volume_axis"]
    groups = [
        (
            "bosenick",
            range(434, 453),
            "bosenick-geiger-1997-table-2.csv",
            "molar_volume_cm3_mol",
            (0, 0, 1),
            [125.293],
        ),
        (
            "du",
            range(460, 473),
            "du-2015-preprint-table-2.csv",
            "cell_volume_a3",
            (0.502, 0, 1),
            [1664.2],
        ),
        (
            "milani",
            range(473, 508),
            "milani-2017-table-2.csv",
            "cell_volume_a3",
            (1, 0.502, 0),
            [1663.98, 1664.01, 1664.15],
        ),
        (
            "rodehorst",
            range(508, 518),
            "rodehorst-2002-table-4a.csv",
            "lattice_a_angstrom",
            (0, 1, 1),
            [11.850**3],
        ),
    ]
    extracted, result = [], {}
    for name, indices, filename, key, color, candidates in groups:
        original = observations(filename)
        rows = []
        for index in indices:
            drawing = drawings[index]
            rect = drawing["rect"]
            if not (
                len(drawing["items"]) == 4
                and 3 < rect.width < 3.4
                and 3 < rect.height < 3.4
                and np.allclose(drawing["color"], color, atol=1e-6)
                and np.allclose(drawing["fill"], color, atol=1e-6)
            ):
                raise ValueError(f"Unexpected marker {index}")
            x, y = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
            t, v = tc["slope"] * x + tc["intercept"], vc["slope"] * y + vc["intercept"]
            matches = [r for r in original if abs(float(r["temperature_k"]) - t) < 2]
            if len(matches) > 1:
                same_x = [
                    i
                    for i in indices
                    if abs((drawings[i]["rect"].x0 + drawings[i]["rect"].x1) / 2 - x)
                    < 0.01
                ]
                # Repeats have identical T; match ascending measured/plotted V.
                rank = sorted(same_x, key=lambda i: -drawings[i]["rect"].y0).index(
                    index
                )
                r = sorted(matches, key=lambda r: float(r[key]))[rank]
            elif len(matches) == 1:
                r = matches[0]
            else:
                raise ValueError(f"No original temperature matches path {index}")
            raw = float(r[key]) ** (3 if name == "rodehorst" else 1)
            rows.append(
                dict(
                    dataset=name,
                    drawing_index=index,
                    source_row=int(r["source_row"]),
                    original_temperature_k=float(r["temperature_k"]),
                    digitized_temperature_k=t,
                    digitized_molar_volume_cm3_mol=v,
                    original_volume=raw,
                    original_volume_unit="cm3/mol" if name == "bosenick" else "A3",
                )
            )
        raw = np.array([r["original_volume"] for r in rows])
        plotted = np.array([r["digitized_molar_volume_cm3_mol"] for r in rows])
        multiplier = float(raw @ plotted / (raw @ raw))
        result[name] = dict(
            markers=len(rows),
            source_table=filename,
            source_table_sha256=hashlib.sha256(
                (DATA / filename).read_bytes()
            ).hexdigest(),
            inferred_normalization_denominator=125.35 / multiplier,
            denominator_unit=rows[0]["original_volume_unit"],
            denominator_resolution_estimate=0.003 * (125.35 / multiplier) / 125.35,
            profiled_plot_residual_cm3_mol=summarize(plotted - multiplier * raw),
            fixed_anchor_comparisons=[
                dict(
                    denominator=a,
                    residual_cm3_mol=summarize(plotted - 125.35 * raw / a),
                )
                for a in candidates
            ],
            source_rows_plotted=[r["source_row"] for r in rows],
        )
        extracted.extend(rows)
    receipt = dict(
        recovered_on="2026-11-12",
        source="Angel et al. (2022), final Figure 4, PDF page 8",
        source_doi="10.1007/s00410-022-01918-5",
        pdf_sha256=PDF_SHA256,
        method="Same affine axes and 2 K / 0.003 cm3/mol resolution bounds as figure4-recovery.json. Match primary row identities by T and ordered V for repeats. Profile one scalar Vplot/rawV per dataset, without fitting EOS coefficients.",
        qualification="Inferred denominators describe plotted processing, not recovered measured anchors or exact author fit inputs. Resolution errors are shared and not divided by sqrt(N). Experimental and preliminary-EOS covariance remain unknown.",
        groups=result,
        csv="figure4-normalization-points.csv",
        findings=[
            "Nineteen blue Bosenick markers match original rows 2-20 with fixed Vm295=125.293 cm3/mol denominator; original 295 K row is not present in these marker paths. This supports an anchor omission but does not prove the fit mask; S1 says 20 used.",
            "All 35 Milani and all 13 Du markers map to original rows using a constant denominator inferred from the plotted series. This narrows normalization without uniquely identifying its original preliminary calculation.",
            "Rodehorst's 156 K marker disagrees with the literal Table 4A a=11.842 A by about 0.03 cm3/mol under a common normalization. Preserve the original table; do not repair it to match the plot.",
        ],
    )
    with (DATA / receipt["csv"]).open("w", newline="", encoding="utf-8") as stream:
        writer = csv.DictWriter(stream, fieldnames=extracted[0], lineterminator="\n")
        writer.writeheader()
        writer.writerows(extracted)
    (DATA / "normalization-recovery.json").write_text(
        json.dumps(receipt, indent=2, allow_nan=False) + "\n"
    )
    return receipt


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path, required=True)
    recover(parser.parse_args().pdf)
