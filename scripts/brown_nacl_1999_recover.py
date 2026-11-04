#!/usr/bin/env python3
"""Recover source numbers without treating calculated EOS tables as observations.

Run from the repository root. PDF retrieval is deliberately separate from this
offline, deterministic transcription step. Source page images were inspected.
"""

import csv
import hashlib
import json
import re
import subprocess
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
OUT = ROOT / "output/material-candidates/brown-nacl-1999"
SRC = OUT / "sources"
PREFIX = "brown-nacl-1999-"


def write_csv(name, fields, rows):
    path = DATA / (PREFIX + name + ".csv")
    with path.open("w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow(fields)
        writer.writerows(rows)
    return path


def recover():
    # Table 1, pp. 519-520: 19 specimens/runs. Ratios reference 25 C/1 atm.
    # Preserve printed digits, including the 29.50 kbar row on p.519.
    bk = [
        (
            25.387,
            25.354,
            [25, 100, 200, 300, 400, 500],
            [1.06, 3.08, 5.72, 8.38, 10.95, 13.49],
            [0.9958, 0.9961, 0.9969, 0.9978, 0.9989, 1.0000],
        ),
        (25.022, 24.966, [25, 500], [1.73, 14.30], [0.9929, 0.9970]),
        (
            24.775,
            24.690,
            [25, 100, 200, 300, 400, 500],
            [2.63, 4.85, 7.51, 10.12, 12.65, 15.16],
            [0.9890, 0.9893, 0.9900, 0.9910, 0.9921, 0.9931],
        ),
        (
            25.386,
            25.265,
            [25, 100, 200, 300, 400, 500],
            [3.78, 5.78, 8.45, 11.08, 13.64, 16.18],
            [0.9848, 0.9851, 0.9859, 0.9868, 0.9879, 0.9889],
        ),
        (
            24.775,
            24.602,
            [25, 100, 200, 300, 400, 500],
            [5.62, 7.67, 10.39, 13.00, 15.55, 18.05],
            [0.9778, 0.9781, 0.9788, 0.9797, 0.9808, 0.9818],
        ),
        (
            25.411,
            25.189,
            [25, 100, 200, 300, 400, 500],
            [7.25, 9.39, 12.02, 14.59, 17.14, 19.61],
            [0.9721, 0.9724, 0.9731, 0.9740, 0.9751, 0.9762],
        ),
        (
            25.392,
            25.113,
            [25, 100, 200, 300, 400, 500],
            [9.26, 11.39, 14.01, 16.62, 19.19, 21.67],
            [0.9650, 0.9652, 0.9660, 0.9669, 0.9679, 0.9690],
        ),
        (
            25.376,
            25.034,
            [25, 100, 200, 300, 400, 500],
            [11.74, 13.71, 16.34, 18.86, 21.44, 23.92],
            [0.9572, 0.9574, 0.9582, 0.9591, 0.9602, 0.9612],
        ),
        (
            25.370,
            24.933,
            [25, 100, 200, 300, 400, 500],
            [15.55, 17.54, 20.17, 22.82, 25.29, 27.75],
            [0.9453, 0.9456, 0.9463, 0.9472, 0.9483, 0.9493],
        ),
        (
            25.019,
            24.547,
            [25, 400, 500],
            [17.29, 27.08, 29.50],
            [0.9401, 0.9430, 0.9440],
        ),
        (
            25.334,
            24.831,
            [25, 100, 200, 300, 400, 500],
            [18.36, 20.50, 23.11, 25.59, 28.10, 30.50],
            [0.9370, 0.9373, 0.9380, 0.9389, 0.9400, 0.9410],
        ),
        (
            25.324,
            24.788,
            [25, 100, 200, 300, 400, 500],
            [19.95, 21.94, 24.54, 27.03, 29.50, 31.85],
            [0.9329, 0.9332, 0.9339, 0.9348, 0.9358, 0.9369],
        ),
        (
            25.410,
            24.841,
            [25, 100, 200, 300, 400],
            [21.30, 23.36, 25.96, 28.52, 31.09],
            [0.9290, 0.9293, 0.9300, 0.9309, 0.9319],
        ),
        (
            25.409,
            24.788,
            [25, 100, 200, 300],
            [23.59, 25.63, 28.23, 30.82],
            [0.9226, 0.9229, 0.9236, 0.9245],
        ),
        (
            25.408,
            24.742,
            [25, 100, 200],
            [25.69, 27.75, 30.36],
            [0.9172, 0.9174, 0.9181],
        ),
        (25.197, 24.514, [25], [26.97], [0.9142]),
        (25.131, 24.411, [25, 100], [28.83, 30.81], [0.9095, 0.9097]),
        (14.641, 14.213, [25], [29.53], [0.9077]),
        (25.131, 24.369, [25], [31.05], [0.9043]),
    ]
    rows = []
    for run, (ln, lw, ts, ps, vs) in enumerate(bk, 1):
        for t, p, v in zip(ts, ps, vs):
            rows.append(
                [
                    run,
                    519 if run <= 10 else 520,
                    f"{ln:.3f}",
                    f"{lw:.3f}",
                    t,
                    f"{p:.2f}",
                    f"{v:.4f}",
                    "measurement_reduced_with_WC_standard",
                ]
            )
    write_csv(
        "boehler-kennedy-1980-table1",
        [
            "run",
            "source_page",
            "nacl_length_mm",
            "wc_length_mm",
            "temperature_c",
            "pressure_kbar",
            "volume_ratio_25c_1atm",
            "kind",
        ],
        rows,
    )

    # Fritz Table 1 p.202. These are published experimental reductions.
    # Us and standard Us measured; Up, P, rho derived by impedance matching.
    # Two OCR errors (Us=5.59 and 5.94) resolved against the page image.
    shock100 = """4.03 5.67 .33 29 2.358
4.10 5.75 .40 35 2.397
4.12 5.79 .44 39 2.423
4.31 5.94 .58 54 2.502
4.36 5.96 .60 57 2.510
4.59 6.17 .79 78 2.613
4.59 6.19 .81 80 2.625
4.95 6.45 1.03 111 2.735
4.89 6.45 1.04 110 2.748
4.99 4.96 1.07 116 2.755
4.99 6.50 1.09 117 2.766
4.99 6.54 1.12 121 2.789
5.01 6.55 1.13 122 2.793
5.51 6.92 1.44 172 2.932
5.43 6.91 1.44 170 2.949
5.47 6.92 1.45 171 2.942
5.47 6.92 1.45 171 2.941
5.59 5.41 1.53 184 2.977
5.67 7.07 1.57 193 2.992
5.84 5.58 1.70 215 3.052
5.94 7.36 1.82 234 3.121
5.96 7.37 1.84 237 3.129
6.12 7.47 1.92 254 3.151
6.12 7.47 1.92 254 3.152
6.15 7.51 1.95 259 3.166
6.15 7.51 1.95 259 3.167
6.13 7.52 1.96 260 3.185
6.20 7.54 1.97 264 3.170"""
    shock111 = """4.08 5.75 .40 35 2.399
4.32 5.96 .60 56 2.515
4.31 6.00 .63 59 2.538
4.59 6.19 .81 80 2.626
4.59 6.20 .82 81 2.632
4.96 6.44 1.03 110 2.729
4.98 4.96 1.07 115 2.757
5.01 6.50 1.08 117 2.760
4.99 6.50 1.08 117 2.765
4.99 6.53 1.11 119 2.780
5.03 6.55 1.13 123 2.787
5.02 6.56 1.13 123 2.795
5.45 6.92 1.45 171 2.947
5.47 6.92 1.45 171 2.944
5.57 5.41 1.53 184 2.981"""
    rows = []
    for group, block in [(100, shock100), (111, shock111)]:
        for i, line in enumerate(block.splitlines(), 1):
            us, ust, up, p, rho = line.split()
            orient = 110 if (group == 100 and i == 17) else group
            std = (
                "Cu"
                if (group == 100 and i in [10, 18, 20])
                or (group == 111 and i in [7, 15])
                else "2024Al"
            )
            rows.append(
                [
                    len(rows) + 1,
                    orient,
                    std,
                    us,
                    ust,
                    up,
                    p,
                    rho,
                    "published_impedance_match_reduction",
                ]
            )
    write_csv(
        "fritz-1971-table1",
        [
            "row",
            "orientation",
            "standard",
            "shock_velocity_km_s",
            "standard_shock_velocity_km_s",
            "particle_velocity_km_s",
            "pressure_kbar",
            "density_g_cm3",
            "kind",
        ],
        rows,
    )

    # Published coefficients are never refitted in this file.
    coefficients = {
        "enck_dommel_1965": {
            "doi": "10.1063/1.1714229",
            "location": "pp.841-842, Eqs.(3),(5)-(7)",
            "alpha_definition": "DeltaL/(L_room*DeltaT), 50 C intervals; not instantaneous dlnL/dT",
            "temperature_unit": "Celsius",
            "alpha_unit": "1/K",
            "low_ascending_coefficients": [3.836e-5, 2.979e-8, 1.167e-11],
            "high_ascending_coefficients": [6.594e-5, -7.499e-8, 1.100e-10],
            "global_ascending_coefficients": [3.942e-5, 2.012e-8, 2.809e-11],
            "junction_c": 550,
            "measured_interval_c": [50, 765],
            "maximum_estimated_relative_error": 0.015,
            "individual_measurements": 77,
            "heating_curves": 12,
            "samples": 3,
            "individual_rows_published": False,
            "low_polynomial_points": 45,
            "low_fitting_interval_c": [50, 600],
            "high_polynomial_points": 32,
            "high_fitting_interval_c": [450, 765],
        },
        "fritz_1971": {
            "brown_citation_year": 1972,
            "verified_volume_issue_date": "March 1971",
            "source_title": "The Hugoniot Equation of State of Sodium Chloride in the Sodium Chloride Structure",
            "location": "NBS SP326 pp.201-208; Table1 p.202; Table2 and Eq.(1) p.203",
            "linear_us_ascending_coefficients": [3.528, 1.343],
            "linear_printed_parameter_uncertainties": [0.012, 0.009],
            "quadratic_us_ascending_coefficients": [3.403, 1.5422, -0.07345],
            "quadratic_intercept_constraint": "Haussuhl bulk sound speed; not a free fitted intercept",
            "us_fit_rms_km_s": 0.05,
            "rho_initial_g_cm3": 2.1645,
            "molar_mass_g_mol": 58.443,
            "avogadro": 6.02252e23,
            "theta_293k_k": 276,
            "gamma_293k": 1.6044,
            "ks_293k_gpa": 25.07,
            "kt_293k_gpa": 23.74,
            "cp_293k_j_g_k": 0.8627,
            "cv_293k_j_g_k": 0.8168,
            "quadratic_isotherm_max_gpa": 25,
            "linear_isotherm_range_gpa": [20, 32],
            "raw_up_measurement": False,
            "pressure_accuracy_statement": "approximately 1 percent conditional on accepted standard Hugoniots",
        },
        "boehler_kennedy_1980": {
            "doi": "10.1016/0022-3697(80)90183-3",
            "location": "Table1 pp.519-520",
            "runs": 19,
            "pressure_relative_absolute_uncertainty": 0.004,
            "volume_change_relative_absolute_uncertainty": 0.007,
            "temperature_absolute_uncertainty_k": 0.3,
            "length_absolute_accuracy_mm": 0.001,
            "uncertainty_note": "bounds/accuracy statements, not published independent Gaussian row sigmas; volume error applies to compression/change, not entire V/V0",
            "table3_kind": "calculated integrated fitted isotherms; never fitted as observations",
        },
    }
    (DATA / (PREFIX + "published-coefficients.json")).write_text(
        json.dumps(coefficients, indent=2) + "\n"
    )

    # Table 2 p.145, independently measured elastic constants; Cp and alpha
    # reproduced from JANAF (1971) and Enck-Dommel, respectively. Keep the
    # printed V column even though its first value conflicts with p.143 density.
    yam = """294 496.7 131.8 127.8 253.5 240.2 145.2 149.6 117.0 .4701 8.635
338 480.4 133.2 126.8 248.9 233.4 142.1 145.5 121.1 .4726 8.792
381 464.6 133.9 124.9 244.1 226.5 138.4 141.1 125.1 .4751 8.913
419 451.3 134.8 123.7 240.3 220.7 135.5 137.5 128.9 .4774 9.015
449 440.8 135.0 122.6 236.9 215.9 133.2 134.7 131.9 .4793 9.095
479 430.8 135.7 121.6 234.1 211.4 130.8 132.0 135.0 .4812 9.173
511 420.3 136.2 120.5 230.9 206.6 128.3 129.1 138.3 .4833 9.260
544 407.4 134.4 119.3 225.4 199.8 125.7 126.2 141.9 .4856 9.350
563 399.5 134.2 118.4 222.7 196.3 123.7 124.1 143.9 .4869 9.400
571 397.7 133.6 118.3 221.6 194.8 123.4 123.8 144.8 .4874 9.420
596 388.1 132.0 117.3 217.4 189.8 121.4 121.6 147.6 .4892 9.485
637 374.4 131.5 115.7 212.4 183.0 117.9 118.0 152.3 .4922 9.605
675 361.5 130.6 114.1 207.6 176.6 114.6 114.6 156.7 .4951 9.715
702 352.8 130.3 113.0 204.5 172.4 112.3 112.3 159.9 .4972 9.800
745 338.5 129.4 111.2 199.1 165.5 108.5 108.6 165.1 .5008 9.950
766 331.8 129.2 110.2 196.7 162.3 106.5 106.6 167.6 .5025 10.030"""
    write_csv(
        "yamamoto-1987-table2",
        [
            "temperature_k",
            "c11_kbar",
            "c12_kbar",
            "c44_kbar",
            "ks_kbar",
            "kt_kbar",
            "g_reuss_kbar",
            "g_voigt_kbar",
            "alpha_1e6_per_k",
            "printed_specific_volume_cm3_g",
            "cp_1e6_erg_g_k",
        ],
        [line.split() for line in yam.splitlines()],
    )
    write_csv(
        "yamamoto-1987-table2-errors",
        [
            "temperature_k",
            "c11_kbar",
            "c12_kbar",
            "c44_kbar",
            "ks_kbar",
            "kt_kbar",
            "g_reuss_kbar",
            "g_voigt_kbar",
        ],
        [
            [294, 3.7, 4.4, 0.2, 4.1, 3.9, 0.2, 0.2],
            [766, 3.6, 3.9, 0.2, 3.8, 3.1, 0.1, 0.1],
        ],
    )
    # Reduced thermodynamic surfaces calculated from ultrasonic frequency fits;
    # these are validation values, not direct independent PVT observations.
    kvals = [
        [25.209, 22.399, 19.546],
        [26.244, 23.500, 20.718],
        [27.245, 24.562, 21.848],
        [28.210, 25.582, 22.932],
        [29.135, 26.558, 23.965],
        [30.017, 27.485, 24.945],
    ]
    write_csv(
        "spetzler-1972-table4-ks",
        ["pressure_kbar", "temperature_k", "ks_gpa"],
        [
            [p, t, k]
            for p, ks in zip([0, 2, 4, 6, 8, 10], kvals)
            for t, k in zip([300, 550, 800], ks)
        ],
    )
    boehler = [
        [25, 1.90, -0.0467, 0.0008440, -0.000008006, 0.015],
        [200, 3.29, -0.0930, 0.0015439, -0.000010363, 0.019],
        [400, 5.18, -0.1464, 0.0022407, -0.000013077, 0.034],
        [600, 7.49, -0.2299, 0.0035334, -0.000019000, 0.048],
        [800, 10.76, -0.3831, 0.0068994, -0.000047336, 0.056],
    ]
    write_csv(
        "boehler-1981-table1",
        ["temperature_c", "a0", "a1", "a2", "a3", "fit_sd_c_per_kbar"],
        boehler,
    )
    # Selected clearly readable calculated states, not original shock shots.
    write_csv(
        "fritz-1971-table4-checkpoints",
        [
            "pressure_kbar",
            "hugoniot_density_g_cm3",
            "calculated_hugoniot_temperature_k",
            "isotherm293_density_g_cm3",
            "isotherm_volume_ratio",
            "isotherm_lattice_ratio",
        ],
        [
            [200, 3.023, 950, 3.081, 0.7025, 0.8890],
            [225, 3.088, 1087, 3.154, 0.6863, 0.8821],
            [250, 3.148, 1234, 3.223, 0.6716, 0.8757],
            [275, 3.205, 1389, 3.288, 0.6583, 0.8699],
            [300, 3.259, 1553, 3.349, 0.6463, 0.8646],
            [320, 3.299, 1690, 3.396, 0.6374, 0.8606],
        ],
    )
    # Original JANAF 2nd ed., PDF leaf 548, ClNa(cr), table dated Sept30 1964.
    # Thermochemical calories; M=58.448 is printed on this historical table.
    # Rows above 1073.8 K are metastable solid extrapolations, not liquid data.
    write_csv(
        "janaf-1971-nacl-cr-cp",
        [
            "temperature_k",
            "cp_cal_mol_k",
            "cp_j_mol_k",
            "historical_molar_mass_g_mol",
            "status",
        ],
        [
            [
                t,
                f"{cp:.3f}",
                f"{cp * 4.184:.6f}",
                58.448,
                "solid" if t < 1073.8 else "metastable_solid_extension",
            ]
            for t, cp in [
                (0, 0),
                (100, 8.349),
                (200, 11.203),
                (298, 12.072),
                (300, 12.080),
                (400, 12.512),
                (500, 12.892),
                (600, 13.259),
                (700, 13.672),
                (800, 14.176),
                (900, 14.787),
                (1000, 15.503),
                (1100, 16.330),
                (1200, 17.200),
                (1300, 17.960),
                (1400, 18.530),
                (1500, 18.940),
            ]
        ],
    )

    # Independently extract Brown comparison tables from the preserved PDF.
    txt = subprocess.check_output(
        ["pdftotext", "-layout", str(SRC / "brown-1999.pdf"), "-"], text=True
    ).replace("⫺", "-")
    first = txt.split("TABLE I. Properties")[1].split("model 关Eqs.")[0]
    t1 = re.findall(
        r"\b(0\.\d{4})\s+(-?\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d+\.\d{2})\s+(\d{3})\b",
        first,
    )
    assert len(t1) == 21
    write_csv(
        "brown-table1-comparison",
        [
            "specific_volume_cm3_g",
            "static_pressure_gpa",
            "kt_gpa",
            "kprime",
            "gamma",
            "theta_k",
        ],
        t1,
    )
    t2 = []
    for line in txt.split("TABLE II. Equation")[1].split("\f")[0].splitlines():
        f = line.split()
        if (
            len(f) >= 3
            and re.fullmatch(r"0\.\d{4}", f[0])
            and re.fullmatch(r"-?0\.\d{4}", f[1])
        ):
            t2.append(f[:2] + [""] * (12 - len(f)) + f[2:])
    assert len(t2) == 63 and sum(bool(x) for r in t2 for x in r[2:]) == 571
    write_csv(
        "brown-table2-comparison",
        ["specific_volume_cm3_g", "compression"]
        + [f"pressure_{t}k_gpa" for t in range(300, 1201, 100)],
        t2,
    )

    manifest = {"recovery_date": "2026-10-16", "source_files": {}, "datasets": {}}
    for f in sorted(SRC.glob("*.pdf")):
        assert f.read_bytes().startswith(b"%PDF"), f
        manifest["source_files"][f.name] = {
            "sha256": hashlib.sha256(f.read_bytes()).hexdigest(),
            "bytes": f.stat().st_size,
        }
    for f in sorted(DATA.glob(PREFIX + "*")):
        if f.is_file():
            manifest["datasets"][f.name] = {
                "sha256": hashlib.sha256(f.read_bytes()).hexdigest()
            }
    (OUT / "recovery-manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    print(
        json.dumps(
            {
                "bk_rows": sum(len(x[2]) for x in bk),
                "shock_rows": len(rows),
                "brown_t1": len(t1),
                "brown_t2": len(t2),
            }
        )
    )


if __name__ == "__main__":
    recover()
