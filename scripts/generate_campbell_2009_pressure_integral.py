"""Independently integrate the NNO pressure contribution in Campbell Eq. (6).

Repository maintenance tool. No absolute Gibbs baseline or held output table
is read by generate(). compare_original() is a separate numerical comparison.
The result is a relative log10 fugacity change, NOT a complete Table S5 replay.
"""

from __future__ import annotations

import argparse
import csv
import io
import json
from pathlib import Path

import numpy as np
from scipy.constants import R
from scipy.special import roots_legendre

from scripts.reproduce_campbell_2009_buffers import volume

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos/data/datasets"
SOURCE = DATA / "campbell-2009-pressure-integration-source.json"
OUTPUT = DATA / "campbell-2009-independent-nno-pressure-integral.csv"
ORIGINAL = DATA / "campbell-2009-tables-s4-s5-buffer-results.csv"


def integrate(pressures, temperature, metal, oxide, order=8):
    """Integrate from 0 GPa; cm3/mol times GPa is 1000 J/mol.

    Eight-point Gaussian quadrature is applied to each adjacent pressure
    interval. Source coefficients and the printed factor 0.8686 are retained.
    R uses the current SI value. No pressure-switch or Gibbs law is invented.
    """
    pressures = np.asarray(pressures, dtype=float)
    if (
        pressures.ndim != 1
        or len(pressures) < 2
        or not np.all(np.isfinite(pressures))
        or pressures[0] != 0
        or np.any(np.diff(pressures) <= 0)
        or not np.isfinite(temperature)
        or temperature <= 0
        or order < 2
    ):
        raise ValueError("Require increasing finite pressures from zero and positive T")
    nodes, weights = roots_legendre(order)
    increments = []
    for low, high in zip(pressures[:-1], pressures[1:]):
        midpoint, radius = (high + low) / 2, (high - low) / 2
        differences = [
            volume(p, temperature, oxide) - volume(p, temperature, metal)
            for p in midpoint + radius * nodes
        ]
        increments.append(radius * float(np.dot(weights, differences)))
    return np.r_[0.0, np.cumsum(increments)] * (0.8686 * 1000 / (R * temperature))


def generate(source=SOURCE, order=8):
    """Return independently calculated pressure increments, with no source rows."""
    configuration = json.loads(source.read_bytes())
    p_start, p_stop, p_step = configuration["grid"]["pressure_gpa"]
    t_start, t_stop, t_step = configuration["grid"]["temperature_k"]
    pressures = list(range(p_start, p_stop + 1, p_step))
    stream = io.StringIO(newline="")
    writer = csv.writer(stream, lineterminator="\n")
    writer.writerow(
        ["temperature_k", "pressure_gpa", "delta_log10_fugacity_from_zero_gpa"]
    )
    for temperature in range(t_start, t_stop + 1, t_step):
        values = integrate(
            pressures,
            temperature,
            configuration["metal_parameters"],
            configuration["oxide_parameters"],
            order,
        )
        for pressure, value in zip(pressures, values):
            writer.writerow([temperature, pressure, format(value, ".12g")])
    return stream.getvalue().encode()


def compare_original(payload, original=ORIGINAL):
    """Compare relative differences only; never supply the missing Gibbs baseline."""
    with original.open(newline="", encoding="utf-8") as stream:
        rows = [r for r in csv.DictReader(stream) if r["source_table"] == "Table S5"]
    reference = {
        float(r["temperature_k"]): float(r["log10_fugacity"])
        for r in rows
        if float(r["pressure_gpa"]) == 0
    }
    original_values = {
        (float(r["temperature_k"]), float(r["pressure_gpa"])): float(
            r["log10_fugacity"]
        )
        - reference[float(r["temperature_k"])]
        for r in rows
    }
    calculated = {
        (float(r["temperature_k"]), float(r["pressure_gpa"])): float(
            r["delta_log10_fugacity_from_zero_gpa"]
        )
        for r in csv.DictReader(io.StringIO(payload.decode()))
    }
    if calculated.keys() != original_values.keys():
        raise ValueError("Original and generated pressure grids differ")
    residual = np.array([v - original_values[k] for k, v in calculated.items()])
    return {
        "rows": len(residual),
        "rmse_delta_log10_fugacity": float(np.sqrt(np.mean(residual**2))),
        "maximum_absolute_difference": float(np.max(abs(residual))),
        "scope": "Pressure contribution only; original zero-pressure values used only to difference the comparison, never as generator inputs or distributed generated values.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    payload = generate()
    if args.check:
        # Numerical computation can differ in final floating-point bits across builds.
        actual = list(csv.DictReader(io.StringIO(OUTPUT.read_text())))
        expected = list(csv.DictReader(io.StringIO(payload.decode())))
        if len(actual) != len(expected) or any(
            a["temperature_k"] != e["temperature_k"]
            or a["pressure_gpa"] != e["pressure_gpa"]
            or abs(
                float(a["delta_log10_fugacity_from_zero_gpa"])
                - float(e["delta_log10_fugacity_from_zero_gpa"])
            )
            > 1e-9
            for a, e in zip(actual, expected)
        ):
            raise SystemExit("Independent pressure-integral table is stale")
    else:
        OUTPUT.write_bytes(payload)
    print(json.dumps(compare_original(payload), indent=2))


if __name__ == "__main__":
    main()
