#!/usr/bin/env python3
"""Export visually checked Co/Ni/Ag cells from AIP Handbook (1972).

Repository maintenance tool: source transcription only, no EOS refitting or
catalogue writes. The original full book remains at its user-supplied location.
"""

import argparse
import csv
import hashlib
import json
from collections import Counter
from decimal import Decimal
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/dewaele-2008-ruby-construction"
PDF_SHA256 = "dda86a046830f63d95933b01e5b56bce9a46572a038ba37e68ceee9f2e059dad"

# Table 4d-12, printed pp.4-100/4-101, PDF pages 621/622 (one-based).
# P in kbar; Co, Ni and Ag values are the printed V/V0, not atomic volumes.
CELLS = """15 .993 .992 .987
20 .990 .990 .982
25 .988 .987 .978
30 .985 .985 .974
35 .983 .982 .970
40 .981 .980 .966
45 .978 .978 .963
50 .976 .975 .959
60 .972 .971 .952
70 .967 .966 .945
80 .963 .962 .939
90 .959 .958 .932
100 .955 .954 .926
120 .947 .946 .915
140 .940 .938 .905
160 .932 .931 .895
180 .925 .924 .886
200 .919 .917 .877
220 .912 .911 .869
240 .906 .904 .861
260 .900 .898 .854
280 .894 .893 .847
300 .888 .887 .840
320 .883 .881 .834
340 .877 .876 .828
360 .872 .871 .822
380 .867 .866 .816
400 .862 .861 .811
420 .857 .857 .806
440 .853 .852 .800
460 .848 .848 .796
480 .844 .843 .791
500 .839 .839 .786
550 .829 .829 .775
600 .819 .820 .765
650 .810 .811 .756
700 .801 .802 .747
750 .792 .794 .739
800 .784 .786 .731
850 .777 .779 .724
900 .769 .772 .717
950 .762 .765 .710
1000 .755 .758 .704
1200 .730 .735 .681
1400 - - .662
1600 - - .645
1800 - - .630
2000 - - .617"""


def write_json(name, value):
    (OUT / name).write_text(json.dumps(value, indent=2) + "\n")


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pdf", type=Path)
    args = parser.parse_args()
    if args.pdf:
        assert hashlib.sha256(args.pdf.read_bytes()).hexdigest() == PDF_SHA256

    source = {}
    for line in CELLS.splitlines():
        pressure, *values = line.split()
        source[int(pressure)] = values
    # Retain every printed pressure row, including explicit empty cells. The
    # 5-kbar row is printed on the Co/Ni page but absent from the Ag page.
    co_ni_grid = [5, 10, *source, 2500, 3000, 3500, 4000, 4500]
    ag_grid = [10, *source]
    rows = []
    for column, material in enumerate(["Co", "Ni", "Ag"]):
        grid = ag_grid if material == "Ag" else co_ni_grid
        for pressure in grid:
            value = source.get(pressure, ["-", "-", "-"])[column]
            if value != "-":
                value = "0" + value
            page = "4-101" if material == "Ag" else "4-100"
            rows.append(
                {
                    "material": material,
                    "table": "4d-12",
                    "printed_page": page,
                    "pdf_page_one_based": 622 if material == "Ag" else 621,
                    "P_kbar_source": f"{pressure:,}",
                    "P_GPa_derived": str(Decimal(pressure) / 10),
                    "V_over_V0_source": "" if value == "-" else value,
                    "source_cell_status": "blank_dot_leaders"
                    if value == "-"
                    else "populated",
                    "T_C_source": "25",
                    "T_K_derived": "298.15",
                    "kind": "published_shock_reduced_calculated_isotherm",
                }
            )
    points = [row for row in rows if row["source_cell_status"] == "populated"]
    counts = Counter(row["material"] for row in points)
    assert counts == {"Co": 44, "Ni": 44, "Ag": 48}
    assert len(rows) == 159 and len(points) == 136
    for material in counts:
        values = [
            Decimal(row["V_over_V0_source"])
            for row in points
            if row["material"] == material
        ]
        assert all(0 < value < 1 for value in values)
        assert all(b < a for a, b in zip(values, values[1:]))
    write_csv("aip-1972-table4d12-co-ni-ag-cells.csv", rows)
    write_csv("aip-1972-table4d12-co-ni-ag-isotherms.csv", points)

    write_json(
        "aip-1972-co-ni-ag-source.json",
        {
            "citation": "G. C. Kennedy and R. N. Keeler, Compressibility, in American Institute of Physics Handbook, 3rd ed., D. E. Gray (ed.), McGraw-Hill, 1972",
            "isbn": "007001485X",
            "pdf_sha256": PDF_SHA256,
            "user_supplied_pdf": str(args.pdf.resolve()) if args.pdf else None,
            "pdf_total_pages": 2359,
            "locations": {
                "title_and_edition_pdf_page": 4,
                "copyright_1972_pdf_page": 5,
                "methods_printed_pages": ["4-96", "4-97", "4-98"],
                "methods_pdf_pages": [617, 618, 619],
                "Co_Ni_table_printed_page": "4-100",
                "Co_Ni_table_pdf_page": 621,
                "Ag_table_printed_page": "4-101",
                "Ag_table_pdf_page": 622,
            },
            "table": "4d-12, Relative volumes of solids at 25 C",
            "units_and_normalization": "P in kbar; V/V0 dimensionless; GPa=P_kbar/10; T=25 C=298.15 K",
            "phase_annotation": "Table labels elements only; phase identity not explicitly supplied for these columns. Do not infer a crystallographic model from the table alone.",
            "kind": "Published calculated isotherms reduced from experimental shock data, not original shock observations or independent static volume measurements",
            "method": {
                "shock_data_ancestor": "M. van Thiel, Compendium of Shock Wave Data, UCRL-50108, 1966 (Sec.4d-2 ref.1)",
                "experimental_laboratories": "Livermore, Los Alamos, Ballistics Research Laboratory and Soviet Union; selected for completeness, accuracy and absence of large error sources (p.4-96)",
                "equation_4d3": "PH(V)-P=gamma(V)*(EH(V)-E)/V",
                "equation_4d4": "EH(V)=E0+[PH(V)+P0]*(V0-V)/2",
                "gamma_law": "Volume-only gamma linked to curvature of the 0 K pressure curve by Eq.4d-6, t=1 (Dugdale-MacDonald); not constant gamma or constant rho*gamma",
                "thermal_step": "Solve 0 K curve and gamma(V), then add the 25 C thermal energy correction via Eq.4d-7",
                "strength": "Finite yield-strength corrections neglected for these metals; source applies a strength correction to Al2O3 separately",
                "phase_limit": "Source warns about phase changes and melting along the Hugoniot; the approximation is not experimental proof of a single phase",
            },
            "uncertainty": {
                "source_accuracy_statement": "p.4-98 estimates calculated 25 C isotherms probably accurate to +/-5 percent in pressure, certainly better than +/-10 percent",
                "semantics": "Qualitative method-level estimate, not per-row standard deviations, confidence intervals or covariance; no per-cell uncertainty supplied",
                "OCR_caution": "Printed 5 percent is repeatedly misread as .5 percent in the OCR. Checked visually; do not report 0.5 percent. Likewise printed 25 C can appear as 2.5 C.",
                "printed_volume_precision": "Three decimals in V/V0; numerical display precision is not a supplied measurement error",
            },
            "counts": {
                "all_printed_cells": len(rows),
                "populated": len(points),
                "blank": len(rows) - len(points),
                "populated_by_material": dict(counts),
            },
            "P_coverage_GPa": {"Co": [1.5, 120], "Ni": [1.5, 120], "Ag": [1.5, 200]},
            "visual_checks": [
                "Complete relevant pages inspected, including table titles and column alignment",
                "Ag at 400 kbar is 0.811; OCR 0.1811 rejected using printed pixels",
                "Co/Ni lower 5/10-kbar cells are empty; no zero-pressure anchor synthesized",
            ],
            "unrecovered": [
                "Original selected shock experiments and pointwise reduction arrays from UCRL-50108",
                "Per-row covariance or statistically defined uncertainties",
                "Dewaele sample mask, interpolation grid, fitting weights, averaging and original worksheet",
            ],
            "scope": "Source recovery only; no new refit, graph correction, EOS recommendation or catalogue change",
        },
    )
    print(
        json.dumps(
            {
                "populated": dict(counts),
                "total": len(points),
                "empty_cells_preserved": len(rows) - len(points),
            }
        )
    )


if __name__ == "__main__":
    main()
