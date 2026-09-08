"""Normalize Fei's official tables and reproduce published branches or the fp20 LS refit.

Run with --write to regenerate the CSVs from the byte-preserved auxiliary files.
By default, comparisons use published parameters. --refit-fp20-ls reproduces the
separate Peritheos refit from primary observations, never an EOS-generated grid.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
DOI = "10.1029/2007GL030712"
TABLES = {"mg080fe020o": (1, 2), "mg061fe039o": (2, 3), "mg042fe058o": (3, 4)}
FP20_LS_REFIT = "mg080fe020o_fei_2007_ls_b1_bm3_refit"
COLUMNS = (
    "source_order",
    "sample_name",
    "nacl_lattice_a_angstrom",
    "nacl_lattice_error_angstrom",
    "pressure_gpa",
    "sample_lattice_a_angstrom",
    "sample_lattice_error_angstrom",
    "compression_path",
    "nacl_phase",
    "volume_a3_conventional_cell",
    "volume_error_a3",
    "sample_volume_interpretation",
)


def normalized_csv(material: str) -> str:
    """Preserve printed values/order; derive cubic volumes only in the B1 phase."""
    table, auxiliary = TABLES[material]
    raw = DATA / "fei-2007" / f"grl23467-sup-000{auxiliary}-ts0{table}.txt"
    output = io.StringIO(newline="")
    writer = csv.writer(output, lineterminator="\n")
    writer.writerow(COLUMNS)
    for index, line in enumerate(raw.read_text().splitlines()[1:], 1):
        fields = line.split()
        if not fields:
            continue
        name, nacl_a, nacl_error, pressure, a = fields[:5]
        error = fields[5] if len(fields) == 6 else ""
        # Paragraph 8 / Figure 2: the mw39c block is decompression.
        path = "decompression" if name.startswith("mw39c_") else "compression"
        # Use the actual lattice branch, not a rounded pressure threshold:
        # Table S1 has B1 at 29.30 GPa and B2 at 26.86 GPa.
        phase = "B1" if float(nacl_a) > 4 else "B2"
        cubic = material != "mg042fe058o" or float(pressure) < 44
        volume = f"{float(a) ** 3:.12g}" if cubic else ""
        volume_error = (
            f"{3 * float(a) ** 2 * float(error):.12g}" if cubic and error else ""
        )
        writer.writerow(
            (
                index,
                name,
                nacl_a,
                nacl_error,
                pressure,
                a,
                error,
                path,
                phase,
                volume,
                volume_error,
                "B1_cell" if cubic else "rhombohedral_scalar_a_not_full_cell",
            )
        )
    return output.getvalue()


def selected_rows(rows: list[dict], record: dict) -> tuple[list[dict], str]:
    """Explicit validation selection, not an assertion of Fei's exact row mask."""
    low, high = record["validation_selection"]["pressure_gpa"]
    selected = [
        row
        for row in rows
        if row["compression_path"] == "compression"
        and low <= float(row["pressure_gpa"]) <= high
    ]
    return selected, record["validation_selection"]["notes"]


def bm3_pressure(volume, v0, k0, kp):
    ratio = (v0 / np.asarray(volume, dtype=float)) ** (1 / 3)
    return 1.5 * k0 * (ratio**7 - ratio**5) * (1 + 0.75 * (kp - 4) * (ratio**2 - 1))


def refit_fp20_ls() -> dict:
    """Independent fit to the full reported LS interval, not a tuned cutoff."""
    import scipy

    import peritheos
    from peritheos.eos.rt import BM3
    from peritheos.fitting import fit_rt_eos

    rows = list(csv.DictReader(io.StringIO(normalized_csv("mg080fe020o"))))
    selected = [
        row
        for row in rows
        if row["compression_path"] == "compression"
        and 40 <= float(row["pressure_gpa"]) <= 95.48
    ]
    pressure = np.array([float(row["pressure_gpa"]) for row in selected])
    volume = np.array([float(row["volume_a3_conventional_cell"]) for row in selected])
    fit = fit_rt_eos(
        BM3,
        volume=volume,
        pressure=pressure,
        initial={"V0": 74.2, "K0": 170},
        fixed={"K0_prime": 4},
        absolute_sigma=False,
        max_nfev=5000,
    )
    if not fit.success:
        raise RuntimeError(fit.message)
    return {
        "record_identifier": FP20_LS_REFIT,
        "record_kind": "refit",
        "parameters": {name: fit.parameters[name] for name in ("V0", "K0", "K0_prime")},
        "parameter_errors": {
            "V0": fit.standard_errors["V0"],
            "K0": fit.standard_errors["K0"],
            "K0_prime": None,
        },
        "parameter_covariance": {
            "parameter_order": list(fit.free_parameters),
            "matrix": fit.covariance.tolist(),
        },
        "parameter_correlation": {
            "parameter_order": list(fit.free_parameters),
            "matrix": fit.correlation.tolist(),
        },
        "selected_sample_names": [row["sample_name"] for row in selected],
        "selected_source_rows": [int(row["source_order"]) for row in selected],
        "observations": len(selected),
        "pressure_range_gpa": [float(pressure.min()), float(pressure.max())],
        "degrees_of_freedom": fit.degrees_of_freedom,
        "residual_sum_squares_gpa2": float(np.sum(fit.residuals**2)),
        "pressure_rmse_gpa": float(np.sqrt(np.mean(fit.residuals**2))),
        "source_curve_pressure_rmse_gpa": float(
            np.sqrt(np.mean((bm3_pressure(volume, 74.2, 170, 4) - pressure) ** 2))
        ),
        "software": {
            "name": "Peritheos",
            "version": peritheos.__version__,
            "numpy_version": np.__version__,
            "scipy_version": scipy.__version__,
        },
        "solver_success": bool(fit.success),
        "solver_message": fit.message,
    }


def reproduce() -> list[dict]:
    results = []
    for material, (table, _) in TABLES.items():
        rows = list(csv.DictReader(io.StringIO(normalized_csv(material))))
        document = json.loads(
            (ROOT / f"peritheos/data/materials/{material}.eosmat").read_text()
        )
        for record in document["eos_records"]:
            if (
                record["reference"]["doi"] != DOI
                or record["record_kind"] != "published"
            ):
                continue
            selected, selection = selected_rows(rows, record)
            volumes = np.array(
                [float(row["volume_a3_conventional_cell"]) for row in selected]
            )
            pressures = np.array([float(row["pressure_gpa"]) for row in selected])
            pars = record["eos"]["parameters"]
            predicted = bm3_pressure(volumes, pars["V0"], pars["K0"], pars["K0_prime"])
            results.append(
                {
                    "record": record["identifier"],
                    "table": f"S{table}",
                    "kind": "published_parameter_residual_validation_not_refit",
                    "selection": selection,
                    "observations": len(selected),
                    "pressure_rmse_gpa": float(
                        np.sqrt(np.mean((predicted - pressures) ** 2))
                    ),
                }
            )
    return results


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    parser.add_argument("--refit-fp20-ls", action="store_true")
    args = parser.parse_args()
    if args.write:
        for material, (table, _) in TABLES.items():
            (DATA / f"{material}-fei-2007-table-s{table}-pv.csv").write_text(
                normalized_csv(material)
            )
    print(json.dumps(refit_fp20_ls() if args.refit_fp20_ls else reproduce(), indent=2))
