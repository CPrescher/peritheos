"""Generate separate buffer approximations using only published coefficients.

Repository maintenance tool. The original integration tables are read solely
by compare_original(), never by generate(). No EOS coefficient is refitted.
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import io
import json
import math
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
SOURCE = DATA / "campbell-2009-polynomial-source.json"
OUTPUT = DATA / "campbell-2009-independent-polynomial-grid.csv"
ORIGINAL = DATA / "campbell-2009-tables-s4-s5-buffer-results.csv"


def evaluate(pressure_gpa, temperature_k, coefficients):
    """Evaluate the source polynomial; P in GPa, T in K, result log10(fO2)."""
    if not math.isfinite(pressure_gpa) or not math.isfinite(temperature_k):
        raise ValueError("Pressure and temperature must be finite")
    if not 0 <= pressure_gpa <= 100 or not 1000 <= temperature_k <= 2600:
        raise ValueError("Outside the source polynomial grid domain")
    return (
        sum(c * pressure_gpa**k for k, c in enumerate(coefficients["a"]))
        + sum(c * pressure_gpa**k for k, c in enumerate(coefficients["b"]))
        / temperature_k
    )


def generate(source=SOURCE):
    """Compute from formula, coefficients and declared grid; no table input."""
    configuration = json.loads(source.read_bytes())
    buffer = io.StringIO(newline="")
    writer = csv.writer(buffer, lineterminator="\n")
    writer.writerow(
        ["source_polynomial", "temperature_k", "pressure_gpa", "log10_fugacity"]
    )
    grid = configuration["grid"]
    for table, coefficients in configuration["coefficients"].items():
        for temperature in range(
            grid["temperature_k"]["start"],
            grid["temperature_k"]["stop_inclusive"] + 1,
            grid["temperature_k"]["step"],
        ):
            for pressure in range(
                grid["pressure_gpa"]["start"],
                grid["pressure_gpa"]["stop_inclusive"] + 1,
                grid["pressure_gpa"]["step"],
            ):
                writer.writerow(
                    [
                        table,
                        temperature,
                        pressure,
                        format(evaluate(pressure, temperature, coefficients), ".15g"),
                    ]
                )
    return buffer.getvalue().encode()


def compare_original(payload, original=ORIGINAL):
    """Match every P/T/table key; never use original outputs to generate rows."""
    generated = {
        (
            r["source_polynomial"],
            float(r["temperature_k"]),
            float(r["pressure_gpa"]),
        ): float(r["log10_fugacity"])
        for r in csv.DictReader(io.StringIO(payload.decode()))
    }
    differences = {}
    with original.open(newline="", encoding="utf-8") as stream:
        originals = list(csv.DictReader(stream))
    keys = {
        (r["source_table"], float(r["temperature_k"]), float(r["pressure_gpa"]))
        for r in originals
    }
    if keys != set(generated) or len(originals) != len(generated):
        raise ValueError("Original and calculated grid keys differ")
    for row in originals:
        key = (
            row["source_table"],
            float(row["temperature_k"]),
            float(row["pressure_gpa"]),
        )
        differences.setdefault(key[0], []).append(
            (generated[key] - float(row["log10_fugacity"]), key[1], key[2])
        )
    return {
        table: {
            "rows": len(values),
            "rmse_log10_fugacity": math.sqrt(
                sum(d**2 for d, _, _ in values) / len(values)
            ),
            "maximum_absolute_difference_log10_fugacity": max(
                abs(d) for d, _, _ in values
            ),
            "maximum_difference_temperature_k": max(values, key=lambda x: abs(x[0]))[1],
            "maximum_difference_pressure_gpa": max(values, key=lambda x: abs(x[0]))[2],
        }
        for table, values in differences.items()
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = generate()
    if args.check:
        if OUTPUT.read_bytes() != payload:
            raise SystemExit("Calculated Campbell polynomial grid is stale")
    else:
        OUTPUT.write_bytes(payload)
    print(
        json.dumps(
            {
                "sha256": hashlib.sha256(payload).hexdigest(),
                "comparison": compare_original(payload),
            },
            indent=2,
        )
    )


if __name__ == "__main__":
    main()
