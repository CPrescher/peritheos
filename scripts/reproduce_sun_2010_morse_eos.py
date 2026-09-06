"""Verify the Sun et al. (2010) 50-solid Morse EOS parameter table."""

from __future__ import annotations

import csv
from importlib import resources

import numpy as np

from peritheos.eos.rt import Morse3, SunMorse3, SunMorse4

TABLE = "sun-2010-table2-morse-eos-parameters.csv"
MODELS = {
    "mrs3": Morse3,
    "sms3": SunMorse3,
    "sms4": SunMorse4,
}


def load_rows() -> list[dict[str, str]]:
    resource = resources.files("peritheos").joinpath("data", "datasets", TABLE)
    with resource.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def reference_derivative(eos: object) -> float:
    v0 = float(eos.V0)
    step = v0 * 1.0e-5
    plus = float(eos.bulk_modulus(v0 - step))
    minus = float(eos.bulk_modulus(v0 + step))
    pressure_plus = float(eos.pressure(v0 - step))
    pressure_minus = float(eos.pressure(v0 + step))
    return (plus - minus) / (pressure_plus - pressure_minus)


def verify_table() -> dict[str, float | int]:
    rows = load_rows()
    largest_p0 = 0.0
    largest_k0_error = 0.0
    largest_k0_prime_error = 0.0
    for row in rows:
        v0 = float(row["v0_cm3_mol"])
        for prefix, model in MODELS.items():
            k0 = float(row[f"{prefix}_k0_gpa"])
            k0_prime = float(row[f"{prefix}_k0_prime"])
            eos = model(v0, k0, k0_prime)
            largest_p0 = max(largest_p0, abs(float(eos.pressure(v0))))
            largest_k0_error = max(
                largest_k0_error, abs(float(eos.bulk_modulus(v0)) - k0)
            )
            largest_k0_prime_error = max(
                largest_k0_prime_error,
                abs(reference_derivative(eos) - k0_prime),
            )
    return {
        "rows": len(rows),
        "curves": len(rows) * len(MODELS),
        "largest_p0_gpa": largest_p0,
        "largest_k0_error_gpa": largest_k0_error,
        "largest_k0_prime_error": largest_k0_prime_error,
        "sms4_mean_error_all_gpa": float(
            np.mean([float(row["sms4_mean_pressure_error_gpa"]) for row in rows])
        ),
    }


def main() -> None:
    for key, value in verify_table().items():
        print(f"{key}: {value}")


if __name__ == "__main__":
    main()
