"""Check the eight Marcondes et al. (2020) static BM3 parameterizations."""

from __future__ import annotations

import json

import numpy as np

TABLE_I = {
    "3fp_hs": (74.19, 170.3, 4.15),
    "3fp_ls": (73.66, 172.1, 4.15),
    "6fp_11nn_hs": (74.38, 171.0, 4.16),
    "6fp_11nn_ls": (73.32, 174.8, 4.17),
    "6fp_11nn_ms": (73.86, 171.3, 4.23),
    "6fp_2nn_hs": (74.15, 172.6, 4.05),
    "6fp_2nn_ls": (73.34, 174.6, 4.17),
    "6fp_2nn_ms": (73.85, 173.0, 4.16),
}


def bm3_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def reproduce() -> dict[str, object]:
    ratios = np.array([1.0, 0.9, 0.8, 0.7])
    return {
        "accepted_records": len(TABLE_I),
        "source_data_status": "parameterization_and_plot_only",
        "records": {
            key: {
                "parameters": {"V0": v0, "K0": k0, "K0_prime": k0_prime},
                "volume_ratios": ratios.tolist(),
                "pressures_gpa": bm3_pressure(ratios * v0, v0, k0, k0_prime).tolist(),
            }
            for key, (v0, k0, k0_prime) in TABLE_I.items()
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
