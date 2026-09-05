#!/usr/bin/env python3
"""Reproduce Akber-Knutson et al. (2005) static BM3 envelope checks."""

from __future__ import annotations

import json

ROWS = {
    "bridgmanite_akber_knutson_2005_gga_bm3": (167.00, 235.0, 3.84),
    "mgsio3_post_perovskite_akber_knutson_2005_gga_bm3": (167.32, 204.0, 4.18),
    "mg09375al0125si09375o3_bridgmanite_akber_knutson_2005_gga_bm3": (
        167.44,
        232.0,
        3.86,
    ),
    "mg09375al0125si09375o3_post_perovskite_akber_knutson_2005_gga_bm3": (
        167.16,
        207.0,
        4.13,
    ),
    "al2o3_perovskite_akber_knutson_2005_gga_bm3": (171.32, 205.0, 4.03),
    "al2o3_post_perovskite_akber_knutson_2005_gga_bm3": (167.12, 201.0, 4.29),
    "mgal0125si0875o29375_bridgmanite_akber_knutson_2005_gga_bm3": (
        169.44,
        214.0,
        3.96,
    ),
    "mgal0125si0875o29375_post_perovskite_akber_knutson_2005_gga_bm3": (
        169.32,
        194.0,
        4.15,
    ),
    "mgal00625h00625si09375o3_bridgmanite_akber_knutson_2005_gga_bm3": (
        168.80,
        228.0,
        3.85,
    ),
    "mgal00625h00625si09375o3_post_perovskite_akber_knutson_2005_gga_bm3": (
        166.56,
        228.0,
        3.86,
    ),
}


def bm3_pressure(volume: float, v0: float, k0: float, k0_prime: float) -> float:
    """Evaluate the standard Eulerian third-order Birch-Murnaghan EOS."""
    eta = (v0 / volume) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def reproduce() -> dict[str, object]:
    """Return independent endpoint pressures for every Table 1 EOS row."""
    checks = {}
    for identifier, (v0, k0, k0_prime) in ROWS.items():
        checks[identifier] = {
            "V0": v0,
            "K0": k0,
            "K0_prime": k0_prime,
            "pressure_at_1.05_v0_gpa": bm3_pressure(1.05 * v0, v0, k0, k0_prime),
            "pressure_at_0.65_v0_gpa": bm3_pressure(0.65 * v0, v0, k0, k0_prime),
        }
    return {
        "source_statement": "E(V) grids span approximately 105-65% of V0 and reach just above 200 GPa",
        "records": checks,
    }


if __name__ == "__main__":
    print(json.dumps(reproduce(), indent=2, sort_keys=True))
