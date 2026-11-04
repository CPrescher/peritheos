#!/usr/bin/env python3
"""Recover newly located upstream tables; repository-only, no catalogue writes.

Requires pdfplumber. PDFs are local evidence and are not redistributed.
Run from the repository root with the two PDFs in the dedicated sources folder.
"""

import csv
import hashlib
import json
import re
from pathlib import Path

import pdfplumber

ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "output/dewaele-2008-ruby-construction"
SOURCES = OUT / "sources"
EXPECTED = {
    "hixson-1992.pdf": "ea920df3a1d80ccb992e6077d760d5b2a93e564079d61bb0c5f0c14a0ba14eb6",
    "holmes-1989.pdf": "4b51ba86f70dfee95777fa0ac874a21d37d6b45a9bbfbc5514a9ba7943399248",
}


def write_csv(name, rows):
    with (OUT / name).open("w", newline="") as stream:
        writer = csv.DictWriter(stream, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_json(name, data):
    (OUT / name).write_text(json.dumps(data, indent=2) + "\n")


def main():
    for name, expected in EXPECTED.items():
        assert hashlib.sha256((SOURCES / name).read_bytes()).hexdigest() == expected
    with pdfplumber.open(SOURCES / "hixson-1992.pdf") as pdf:
        pages = [page.extract_text(layout=True) for page in pdf.pages]

    pvt = []
    corrections = []
    for metal, page, table, last_eta in [("Mo", 7, "VI", 0.38), ("W", 8, "VII", 0.36)]:
        text = pages[page - 1].split("TABLE " + table + ".", 1)[1]
        lines = [
            line
            for line in text.splitlines()
            if re.match(r"\s*[-.]?\s*[0O]\.\w+\s+[-\d]", line)
        ]
        expected_count = round((last_eta + 0.02) / 0.02) + 1
        assert len(lines) == expected_count, (metal, lines)
        for index, line in enumerate(lines):
            # Row labels contain OCR letter substitutions. The visually verified
            # printed eta grid is -0.0200, 0.0000, ..., the stated last_eta.
            eta = round(-0.02 + 0.02 * index, 4)
            tokens = re.sub(r"-\s+", "-", line.strip()).split()[1:]
            assert len(tokens) == 10, (metal, eta, tokens)
            values = [float(token) for token in tokens]
            if metal == "Mo" and eta == 0.02:
                assert values[7] == 7.14
                values[7] = 7.74
                corrections.append(
                    {
                        "table": table,
                        "eta": eta,
                        "T_K": 800,
                        "OCR": "7.14",
                        "printed": "7.74",
                    }
                )
            if metal == "W" and eta == 0.18:
                assert values[4] == 91.88
                values[4] = 91.84
                corrections.append(
                    {
                        "table": table,
                        "eta": eta,
                        "T_K": 500,
                        "OCR": "91.88",
                        "printed": "91.84",
                    }
                )
            for temperature, pressure in zip(range(100, 1001, 100), values):
                pvt.append(
                    {
                        "material": metal,
                        "table": table,
                        "printed_page": 1719 + page,
                        "eta_1_minus_V_over_V0": f"{eta:.4f}",
                        "V_over_V0_derived": f"{1 - eta:.4f}",
                        "T_K": temperature,
                        "P_GPa": f"{pressure:.2f}",
                        "kind": "published_model_calculation",
                    }
                )
    write_csv("hixson-1992-tables6-7-pvt.csv", pvt)
    write_csv(
        "hixson-1992-300K-isotherms.csv", [row for row in pvt if row["T_K"] == 300]
    )

    density = []
    for metal, page, table, rho0 in [("Mo", 6, "IV", 10.215), ("W", 7, "V", 19.256)]:
        text = pages[page - 1].split("TABLE " + table + ".", 1)[1]
        for line in text.splitlines():
            if not re.match(r"\s*\d+\.0\s+\d+\.\d{3}\s+\d+\s", line):
                continue
            tokens = line.split()
            pressure = float(tokens[0])
            rho_token = tokens[-2]
            if metal == "W" and pressure == 180:
                assert rho_token == "26.3!7"
                rho_token = "26.357"
                corrections.append(
                    {
                        "table": table,
                        "P_GPa": pressure,
                        "column": "rho_293K",
                        "OCR": "26.3!7",
                        "printed": rho_token,
                    }
                )
            rho = float(rho_token)
            density.append(
                {
                    "material": metal,
                    "table": table,
                    "printed_page": 1719 + page,
                    "T_K": 293,
                    "P_GPa": f"{pressure:.1f}",
                    "rho_g_cm3": rho_token,
                    "rho0_g_cm3": f"{rho0:.3f}",
                    "V_over_V0_derived": f"{rho0 / rho:.12g}",
                    "bulk_sound_speed_km_s": tokens[-1],
                    "kind": "published_model_calculation",
                }
            )
        assert sum(row["material"] == metal for row in density) == 31
    write_csv("hixson-1992-293K-density-isotherms.csv", density)

    # Tables I/II, checked against rendered pages, including the two LLNL rows
    # explicitly attributed by footnote 4. These are shock states, not isotherms.
    shock_text = """Mo I 723 Ta 4.070 177.80 2.21 7.88 14.19 current
Mo I 693 Ta 4.633 214.73 2.51 8.37 14.59 current
Mo I 727 Pt 4.900 270.63 2.98 8.90 15.35 current
Mo I 726 Ta 5.790 292.86 3.18 9.03 15.75 current
Mo I 694 Ta 6.475 347.23 3.55 9.58 16.23 current
Mo I 695 Pt 6.510 409.00 3.97 10.08 16.85 current
Mo I Mo3 Ta 7.784 457.94 4.29 10.45 17.33 Mitchell1991_ref4
Mo I Mo1 Ta 7.981 476.43 4.40 10.60 17.46 Mitchell1991_ref4
W II 796 Ta 4.514 263.88 2.07 6.61 28.07 current
W II 801 Ta 5.088 315.87 2.34 7.01 28.90 current
W II 797 Ta 6.066 408.33 2.82 7.52 30.80 current
W II 799 Ta 6.877 499.50 3.19 8.13 31.68 current
W II 800 Pt 6.870 603.44 3.60 8.70 32.86 current
W II 803 Pt 7.417 675.93 3.91 8.98 34.10 current"""
    shock_fields = "material table shot impactor impactor_velocity_km_s P_GPa up_km_s us_km_s rho_g_cm3 attribution".split()
    shock = [dict(zip(shock_fields, line.split())) for line in shock_text.splitlines()]
    for row in shock:
        row["kind"] = "published_shock_experiment_state"
        row["uncertainty"] = "not_tabulated_per_row"
    write_csv("hixson-1992-tables1-2-shock.csv", shock)

    pt_text = """TaPt1 21.41 16.65 7.784 8.920 0.134 3.422 0.034 653.4 7.64 34.728 0.535
TaPt3 21.42 16.65 7.829 8.937 0.063 3.444 0.024 659.3 5.16 34.853 0.298
TaPt4 21.42 16.68 6.006 7.682 0.030 2.633 0.011 433.4 1.97 32.597 0.135
TaPt5 21.37 16.68 4.492 6.675 0.035 1.958 0.008 279.3 1.31 30.239 0.121
TaPt6 21.40 16.68 4.324 6.535 0.118 1.885 0.019 263.6 2.98 30.072 0.343
TaPt7 21.40 16.67 6.850 8.291 0.032 3.005 0.014 533.3 2.68 33.570 0.159
TaPt8 21.41 16.68 3.803 6.189 0.035 1.652 0.007 218.9 1.06 29.208 0.110"""
    pt_fields = "shot initial_density_g_cm3 impactor_density_g_cm3 impactor_velocity_km_s us_km_s us_uncertainty_2sigma_km_s up_km_s up_uncertainty_2sigma_km_s P_GPa P_uncertainty_2sigma_GPa rho_g_cm3 rho_uncertainty_2sigma_g_cm3".split()
    pt = [dict(zip(pt_fields, line.split())) for line in pt_text.splitlines()]
    for row in pt:
        row.update(
            material="Pt",
            table="III",
            printed_page=2964,
            kind="published_shock_experiment_state",
            uncertainty_semantics="2sigma_as_printed_not_independent_covariance",
        )
    write_csv("holmes-1989-table3-pt-shock.csv", pt)
    write_json(
        "holmes-1989-300K-parameterization.json",
        {
            "doi": "10.1063/1.344177",
            "locations": "pp.2964-2966; Eq.11 and Table IV",
            "kind": "published_parameterization_of_first_principles_isotherm",
            "T_K": 300,
            "P_range_GPa": [0, 550],
            "V0_atomic_units_cubed": 101.9,
            "volume_normalization": "per_atom",
            "equation": "P=P_T*(1-X)/X**2*exp(eta*(1-X)); X=(V/V0)**(1/3)",
            "printed_coefficients": {
                "P_T_GPa": 798.31,
                "eta": 7.2119,
                "B_T_GPa_rounded": 266,
                "B_T_prime_rounded": 5.81,
                "alpha_T_K_inverse": 0.0000261,
            },
            "coefficient_convention": "Use P_T and eta for the printed curve; B_T and B_T_prime are separately rounded, not algebraically exact equivalents.",
            "provenance": "LMTO zero-temperature pressure plus ion and electron thermal terms; Hugoniot checked against shock experiments; McQueen shock-reduced isotherm plotted separately in Fig.5.",
            "uncertainty": "Authors qualify calculated isotherm accuracy at about 10%; not a pointwise standard deviation.",
            "not_recovered": "Dewaele low-pressure resampling, ultrasonic-constrained refit, weights and original numerical support",
        },
    )
    write_json(
        "hixson-1992-reduction-inputs.json",
        {
            "doi": "10.1063/1.351203",
            "locations": "Sec.IV, Table III, pp.1722-1724",
            "law": "Mie-Gruneisen with rho*gamma constant; Debye caloric model; experimentally fitted linear Us=c0+s*up Hugoniot",
            "materials": {
                "Mo": {
                    "rho0_g_cm3": 10.215,
                    "c0_km_s": 5.122,
                    "s": 1.245,
                    "gamma0": 1.58,
                    "model_3nk_J_g_K": 0.272,
                    "theta0_K_rounded": 411,
                    "theta0_over_293K": 1.4030,
                },
                "W": {
                    "rho0_g_cm3": 19.256,
                    "c0_km_s": 4.022,
                    "s": 1.260,
                    "gamma0": 1.65,
                    "model_3nk_J_g_K": 0.140,
                    "theta0_K_rounded": 340,
                    "theta0_over_293K": 1.1612,
                },
            },
            "reference_state": "293 K; use Hugoniot rho0, not the distinct 25 C crystallographic density also printed in Table III",
            "qualifications": [
                "Adjusted model 3nk differs from mass-based Dulong-Petit value.",
                "Phase changes ignored in reduction.",
                "No explicit electron thermal term.",
                "Gamma uncertainty can grow toward 10% at higher compression; authors estimate resulting isotherms accurate to about 2-3%.",
                "Ta/Pt/Cu primary Hugoniot standards are shared ancestors; Mo/W are not independent absolute standards.",
                "Per-row covariance and Dewaele sampling/weighting not recovered.",
            ],
        },
    )

    assert len(pvt) == 410 and len(density) == 62 and len(shock) == 14 and len(pt) == 7
    for metal in ["Mo", "W"]:
        for temperature in range(100, 1001, 100):
            values = [
                float(row["P_GPa"])
                for row in pvt
                if row["material"] == metal and row["T_K"] == temperature
            ]
            assert all(b > a for a, b in zip(values, values[1:]))
    write_json(
        "upstream-recovery-verification.json",
        {
            "date": "2026-10-23",
            "pdf_sha256": EXPECTED,
            "counts": {
                "published_PVT_calculations": 410,
                "subset_at_300K": 41,
                "published_293K_density_points": 62,
                "Mo_W_shock_states": 14,
                "Pt_shock_states_with_2sigma": 7,
            },
            "checks": [
                "PDF hashes fixed",
                "Tables visually inspected against rendered source pages",
                "Dimensions and pressure monotonicity checked",
                "Source calculated tables kept separate from measured shock states",
                "No Dewaele fit or calibration changed",
            ],
            "OCR_readings_corrected_against_pixels": corrections,
            "scope": "Source-table recovery; no original reduction integration or Dewaele calibration workflow reproduction claimed",
        },
    )
    print(
        json.dumps(
            {
                "PVT": len(pvt),
                "density_293K": len(density),
                "shock_Mo_W": len(shock),
                "shock_Pt": len(pt),
                "out": str(OUT),
            }
        )
    )


if __name__ == "__main__":
    main()
