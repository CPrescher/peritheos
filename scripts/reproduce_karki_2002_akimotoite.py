"""Reproduce the Karki and Wentzcovitch (2002) akimotoite EOS audit."""

from __future__ import annotations

import json

from peritheos.eos.rt import BM4

FORMULA_UNITS_PER_CELL = 6.0
PUBLISHED = {
    "static": {
        "temperature_k": 0.0,
        "V0_a3_formula": 43.61,
        "K0_gpa": 210.0,
        "K0_prime": 4.57,
        "K0_double_prime_gpa_inverse": -0.041,
    },
    "300_K": {
        "temperature_k": 300.0,
        "V0_a3_formula": 44.20,
        "K0_gpa": 201.0,
        "K0_prime": 4.64,
        "K0_double_prime_gpa_inverse": -0.042,
    },
    "1000_K": {
        "temperature_k": 1000.0,
        "V0_a3_formula": 45.03,
        "K0_gpa": 182.0,
        "K0_prime": 4.86,
        "K0_double_prime_gpa_inverse": -0.051,
    },
    "2000_K": {
        "temperature_k": 2000.0,
        "V0_a3_formula": 48.59,
        "K0_gpa": 153.0,
        "K0_prime": 5.20,
        "K0_double_prime_gpa_inverse": -0.067,
    },
}


def reproduce() -> dict[str, object]:
    """Check cell-basis conversion and execute each published BM4 curve."""
    records = {}
    for name, row in PUBLISHED.items():
        v0 = row["V0_a3_formula"] * FORMULA_UNITS_PER_CELL
        eos = BM4(
            v0,
            row["K0_gpa"],
            row["K0_prime"],
            row["K0_double_prime_gpa_inverse"],
        )
        records[name] = {
            **row,
            "V0_a3_conventional_z6": v0,
            "pressure_at_v0_gpa": eos.pressure(v0),
            "bulk_modulus_at_v0_gpa": eos.bulk_modulus(v0),
            "pressure_at_0.90_v0_gpa": eos.pressure(0.90 * v0),
        }
    return {
        "formula_units_per_cell": FORMULA_UNITS_PER_CELL,
        "records": records,
        "temperature_trends": {
            "volume_increases_monotonically": all(
                b["V0_a3_formula"] > a["V0_a3_formula"]
                for a, b in zip(PUBLISHED.values(), list(PUBLISHED.values())[1:])
            ),
            "bulk_modulus_decreases_monotonically": all(
                b["K0_gpa"] < a["K0_gpa"]
                for a, b in zip(PUBLISHED.values(), list(PUBLISHED.values())[1:])
            ),
        },
    }


def main() -> None:
    """Print stable JSON output for review and CI."""
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
