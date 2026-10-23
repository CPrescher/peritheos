"""Recover factual tables and Figure 4 vector checkpoints from the supplied PDF.

Repository maintenance tool; requires pdfplumber (not a package dependency).
The copyrighted PDF is never copied into the repository.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import re
from pathlib import Path

import pdfplumber

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"

# Table VI, p. 104102-10: unparenthesized (298 K), then parenthesized
# (static lattice), in atomic A^3, GPa, dimensionless K0-prime.
TABLE6 = """fe lda 10.631 237 1.88 10.586 239 1.88
fe pbe 11.546 181 5.31 11.483 186 5.27
fe pbesol 10.931 206 6.56 10.881 213 6.50
iron lda 9.726 329 4.75 9.697 334 4.74
iron pbe 10.532 224 5.46 10.486 229 5.44
iron pbesol 9.913 297 4.94 9.880 301 4.93
cobalt_hcp lda 10.124 281 2.85 10.071 286 2.84
cobalt_hcp pbe 10.992 208 4.85 10.914 215 4.81
cobalt_hcp pbesol 10.542 243 4.21 10.478 249 4.18
nickel lda 10.099 244 5.04 10.018 254 4.99
nickel pbe 11.008 184 5.24 10.893 195 5.17
nickel pbesol 10.477 217 5.13 10.383 227 5.08
zinc_hcp lda 13.805 82 5.90 13.585 89 5.80
zinc_hcp pbe 15.760 50 6.19 15.365 57 6.04
zinc_hcp pbesol 14.468 71 5.91 14.205 79 5.80
molybdenum lda 15.145 295 4.01 15.063 301 3.99
molybdenum pbe 15.947 255 4.17 15.849 262 4.14
molybdenum pbesol 15.410 280 4.05 15.323 287 4.02
silver lda 16.178 131 5.93 16.046 137 5.87
silver pbe 18.137 82 6.23 17.907 88 6.13
silver pbesol 16.856 110 6.05 16.694 116 5.98
platinum lda 14.898 293 5.65 14.823 302 5.61
platinum pbe 15.899 219 6.20 15.793 228 6.14
platinum pbesol 15.204 274 5.64 15.123 282 5.60
gold lda 16.724 187 5.92 16.603 195 5.85
gold pbe 18.173 132 6.15 17.988 140 6.07
gold pbesol 17.261 169 5.69 17.122 176 5.62"""

# Table IV: V0 fixed; free K0/K0' and source 95% widths; adopted-US-K0
# alternatives (K0 fixed, K0' fitted and its 95% width).
TABLE4 = {
    "fe": [11.760, 171, 4.34, 11, 2.10, [[165, 5.47, 0.51]]],
    "cobalt_hcp": [
        11.096,
        190.5,
        4.38,
        2.9,
        0.18,
        [[184, 4.73, 0.06], [190, 4.37, 0.04]],
    ],
    "nickel": [10.942, 176.7, 5.23, 2.5, 0.09, [[183, 4.99, 0.02]]],
    "zinc_hcp": [15.213, 60.5, 5.69, 0.9, 0.07, [[60, 5.73, 0.01]]],
    "molybdenum": [15.558, 274.1, 3.58, 3.2, 0.12, [[261, 4.06, 0.05]]],
    "silver": [17.021, 103.9, 5.78, 1.4, 0.09, [[101, 5.97, 0.02]]],
}

RSW_COMPARISONS = {
    "cobalt_hcp": [[184, 4.93, 0.04, 42], [190, 4.67, 0.03, 42]],
    "nickel": [[183, 5.13, 0.01, 42]],
    "molybdenum": [[261, 4.18, 0.02, 43]],
    "silver": [[101, 6.19, 0.01, 42]],
    "column_order": [
        "K0_fixed_gpa",
        "K0_prime",
        "K0_prime_error_95pct",
        "source_reference_number",
    ],
    "zinc_hcp": {
        "K0_fixed_gpa": 60,
        "K0_prime_range": [5.60, 5.95],
        "source_references": [45, 42],
    },
    "status": "source comparison values preserved; not promoted to independently validated executable records; exact reduced-shock inputs and graph benchmark not recovered",
}

CALIBRATION_LINEAGE = {
    "table4_input": "ruby_dewaele_2004",
    "new_calibration_output": "ruby_dewaele_2008",
    "ruby_output_parameters": {"A_gpa": 1920, "B": 9.61},
    "ten_metal_comparison": ["Co", "Ni", "Mo", "Ag", "Al", "Cu", "Pt", "Ta", "W", "Au"],
    "excluded_metal": {
        "symbol": "Zn",
        "reason": "uncertain reduced shock EOS, melting and thermal reduction",
    },
    "network_reproduction": "not reproduced: exact complete comparison grid, weights, covariances and reduced-shock inputs unavailable",
    "warning": "The calibration produced by this paper is not the input calibration for its Table IV fits.",
}


def token_number(token):
    return float(token.strip("()共兲")) if token else None


def word_rows(page):
    groups = {}
    for word in sorted(page.extract_words(), key=lambda w: (w["top"], w["x0"])):
        word["text"] = word["text"].replace("(cid:1)", "(").replace("(cid:2)", ")")
        top = next(
            (y for y in groups if abs(y - word["top"]) < 2.2), round(word["top"], 1)
        )
        groups.setdefault(top, []).append(word)
    return groups


def write_csv(path, rows):
    with path.open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def table2(page):
    groups = word_rows(page)
    runs = [None, None]
    bounds = [45, 78, 110, 158, 204, 240, 277, 327, 370, 408, 446, 488, 526, 564]
    rows = []
    for top, words in sorted(groups.items()):
        if top < 175 or top > 752:
            continue
        cells = [""] * 13
        for word in words:
            match = re.fullmatch(r"Run\s*(\d)", word["text"])
            if match:
                runs[0 if word["x0"] < 204 else 1] = int(match[1])
                continue
            x = (word["x0"] + word["x1"]) / 2
            for i, (left, right) in enumerate(zip(bounds, bounds[1:])):
                if left <= x < right:
                    cells[i] += word["text"]
        for side, specs in (
            (0, [("silver", 2, None, None), ("molybdenum", 3, None, None)]),
            (
                1,
                [
                    ("nickel", 6, None, None),
                    ("zinc_hcp", 7, 8, 9),
                    ("cobalt_hcp", 10, 11, 12),
                ],
            ),
        ):
            pc = 0 if side == 0 else 4
            if not re.fullmatch(r"[0-9.共兲()]+", cells[pc]):
                continue
            for material, vc, ac, cc in specs:
                if not cells[vc]:
                    continue
                p = token_number(cells[pc + 1])
                rows.append(
                    {
                        "material": material,
                        "run": runs[side],
                        "pdf_top_pt": top,
                        "pressure_mao_token": cells[pc],
                        "pressure_revised_token": cells[pc + 1],
                        "volume_token": cells[vc],
                        "a_token": cells[ac] if ac is not None else "",
                        "c_token": cells[cc] if cc is not None else "",
                        "ruby_pressure_mao_gpa": token_number(cells[pc]),
                        "ruby_pressure_dewaele_gpa": p,
                        "atomic_volume_a3": token_number(cells[vc]),
                        "a_angstrom": token_number(cells[ac]) if ac is not None else "",
                        "c_angstrom": token_number(cells[cc]) if cc is not None else "",
                        "pressure_from_helium_eos": "共" in cells[pc]
                        or "(" in cells[pc],
                        "fit_included": not (side == 0 and runs[side] == 2 and p > 100),
                    }
                )
    return rows


def table3(page):
    groups = word_rows(page)
    bounds = [112, 150, 195, 255, 300, 350, 406, 455, 500]
    rows = []
    run = 1
    for top, words in sorted(groups.items()):
        text = "".join(w["text"] for w in words)
        if text == "Run3":
            run = 3
        if text == "Run2":
            run = 2
        cells = [""] * 8
        if top < 150 or top > 650:
            continue
        for word in words:
            x = (word["x0"] + word["x1"]) / 2
            for i, (a, b) in enumerate(zip(bounds, bounds[1:])):
                if a <= x < b:
                    cells[i] += word["text"]
        if re.fullmatch(r"[0-9.]+", cells[4]) and re.fullmatch(r"[0-9.]+", cells[5]):
            rows.append(
                {
                    "run": "ref15_run5",
                    "pdf_top_pt": top,
                    "pressure_mao_gpa": "",
                    "pressure_revised_gpa": "",
                    "atomic_volume_a3": float(cells[5]),
                    "c_over_a": cells[6],
                    "phase": "hcp",
                    "pressure_w_gpa": float(cells[4]),
                    "a_w_angstrom": float(cells[7]),
                }
            )
        if not re.fullmatch(r"[0-9.]+", cells[0]) or not re.fullmatch(
            r"[0-9.]+", cells[2]
        ):
            continue
        rows.append(
            {
                "run": run,
                "pdf_top_pt": top,
                "pressure_mao_gpa": float(cells[0]),
                "pressure_revised_gpa": float(cells[1]),
                "atomic_volume_a3": float(cells[2]),
                "c_over_a": cells[3],
                "phase": "hcp" if cells[3] else "bcc",
                "pressure_w_gpa": "",
                "a_w_angstrom": "",
            }
        )
    return rows


def figure4(page):
    # Curve IDs inspected against the complete rendered Figure 4 and legend.
    # Thickness 1.092 pt = with semicore; line/color identifies functional.
    panels = {
        "iron": (
            129.164,
            240.335,
            72.705,
            171.203,
            250,
            7,
            12,
            {"lda": 25, "pbe": 23, "pbesol": 93},
        ),
        "fe": (129.164, 240.335, 72.705, 171.203, 250, 7, 12, {"pbe": 1, "pbesol": 90}),
        "cobalt_hcp": (
            259.779,
            370.951,
            72.705,
            171.203,
            140,
            7.5,
            11.5,
            {"lda": 119, "pbe": 116, "pbesol": 117},
        ),
        "nickel": (
            390.394,
            501.565,
            72.705,
            171.203,
            200,
            7,
            11,
            {"lda": 153, "pbe": 149, "pbesol": 151},
        ),
        "zinc_hcp": (
            129.164,
            240.335,
            188.427,
            286.926,
            175,
            8.5,
            15.5,
            {"lda": 219, "pbe": 217, "pbesol": 263},
        ),
        "molybdenum": (
            259.779,
            370.951,
            188.427,
            286.926,
            150,
            11,
            16,
            {"lda": 274, "pbe": 272, "pbesol": 304},
        ),
        "silver": (
            390.394,
            501.565,
            188.427,
            286.926,
            175,
            10.5,
            18.5,
            {"lda": 312, "pbe": 311, "pbesol": 347},
        ),
        "platinum": (
            128.218,
            240.335,
            304.149,
            402.866,
            150,
            11.5,
            16,
            {"lda": 361, "pbe": 359, "pbesol": 397},
        ),
        "gold": (
            259.779,
            370.951,
            304.149,
            402.648,
            150,
            12,
            18.5,
            {"lda": 412, "pbe": 408, "pbesol": 410},
        ),
    }
    rows = []
    for material, spec in panels.items():
        x0, x1, top, bottom, pmax, vmin, vmax, curves = spec
        for functional, curve_id in curves.items():
            curve = page.curves[curve_id]
            if abs(curve["linewidth"] - 1.092) > 1e-6:
                raise ValueError("Unexpected curve style")
            for vertex, (x, y) in enumerate(curve["pts"]):
                pressure = (x - x0) / (x1 - x0) * pmax
                if pressure < 5:
                    continue  # finite-compression checkpoints, not P(V0).
                rows.append(
                    {
                        "material": material,
                        "functional": functional,
                        "curve_index": curve_id,
                        "vertex": vertex,
                        "pdf_x_pt": x,
                        "pdf_top_pt": y,
                        "pressure_gpa": pressure,
                        "atomic_volume_a3": vmax
                        - (y - top) / (bottom - top) * (vmax - vmin),
                        "volume_reading_halfwidth_a3": 0.04,
                        "pressure_reading_halfwidth_gpa": 1.0,
                    }
                )
    return rows, panels


def recover(pdf_path):
    with pdfplumber.open(pdf_path) as pdf:
        observations = table2(pdf.pages[2])
        iron = table3(pdf.pages[3])
        checkpoints, axes = figure4(pdf.pages[8])
    theory = []
    for line in TABLE6.splitlines():
        material, functional, *tokens = line.split()
        theory.append(
            {
                "material": material,
                "functional": functional,
                "at_298k_tokens": tokens[:3],
                "static_lattice_tokens": tokens[3:],
            }
        )
    write_csv(DATA / "dewaele-2008-table2-complete.csv", observations)
    write_csv(DATA / "dewaele-2008-table3-fe.csv", iron)
    write_csv(DATA / "dewaele-2008-figure4-vector.csv", checkpoints)
    for material in {row["material"] for row in checkpoints}:
        write_csv(
            DATA / f"{material}-dewaele-2008-figure4-vector.csv",
            [row for row in checkpoints if row["material"] == material],
        )
    manifest = {
        "doi": "10.1103/PhysRevB.78.104102",
        "source_pdf_sha256": hashlib.sha256(pdf_path.read_bytes()).hexdigest(),
        "source_version": "Publisher version of record, supplied by user; 13 pages",
        "tables": {"IV": TABLE4, "VI": theory},
        "table4_reduced_shock_wave_comparisons": RSW_COMPARISONS,
        "calibration_lineage": CALIBRATION_LINEAGE,
        "figure4_axes": axes,
        "figure4_note": "Vector path vertices of with-semicore static-lattice calculations; not measured or raw DFT observations. Axis mappings and styles visually checked. Reading halfwidths are diagnostic bounds, not standard errors.",
        "files": {
            p.name: hashlib.sha256(p.read_bytes()).hexdigest()
            for p in DATA.glob("dewaele-2008-*.csv")
        },
        "thermal_pressure_298k_gpa": {
            "fe": 1.0,
            "iron": 1.0,
            "cobalt_hcp": 1.5,
            "nickel": 2.0,
            "zinc_hcp": 1.4,
            "molybdenum": 1.6,
            "silver": 1.1,
            "platinum": 1.5,
            "gold": 1.4,
        },
        "qualifications": [
            "Table IV V0 is fixed after a separate 0-5 GPa determination; exact low-pressure weighting unavailable.",
            "Table II preserves run-2 Ag/Mo rows above 100 GPa but excludes them from Table IV replay.",
            "Tables IV and VI have inconsistent alpha-Fe V0 (11.760 versus 11.770); Table VI also rounds Co/Ni comparison volumes differently. Comparison rows do not overwrite Table IV.",
            "Table VI parentheses are static-lattice coefficients, not uncertainties; 298 K correction includes zero-point motion.",
            "Alpha-Fe LDA has no separately identifiable with-semicore finite-compression curve in Figure 4; preserve its two theory coefficient sets as deferred candidates pending an independent benchmark.",
        ],
    }
    (DATA / "dewaele-2008-source.json").write_text(
        json.dumps(manifest, indent=2) + "\n"
    )
    print(
        f"Recovered {len(observations)} material observations, {len(iron)} Fe rows, {len(checkpoints)} vector vertices, 54 theory triplets"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    recover(parser.parse_args().pdf)
