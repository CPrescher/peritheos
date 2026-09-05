"""Check the complete Leonov et al. (2017) high-spin BM3 records."""

from __future__ import annotations

import json

import numpy as np

BOHR_ANGSTROM = 0.529177210903

TABLE_I = {
    0.000: (144.1, 142.0, 210.0, 73.0),
    0.125: (143.1, 139.0, 205.0, 82.0),
    0.250: (141.3, 137.0, 201.0, 83.0),
    0.375: (139.5, 138.0, 213.0, 77.0),
    0.500: (138.6, 139.0, 200.0, 49.0),
    0.625: (135.5, 142.0, 185.0, 61.0),
    0.750: (133.8, 151.0, 169.0, 52.0),
    0.875: (132.9, 159.0, 158.0, 21.0),
}


def conventional_cell_volume(volume_bohr3_per_formula_unit: float) -> float:
    """Convert the printed B1 formula-unit volume to a conventional Z=4 cell."""
    return 4.0 * volume_bohr3_per_formula_unit * BOHR_ANGSTROM**3


def bm3_pressure(
    volume: np.ndarray, v0: float, k0: float, k0_prime: float
) -> np.ndarray:
    """Evaluate third-order Birch-Murnaghan independently of Peritheos."""
    eta = (v0 / np.asarray(volume, dtype=float)) ** (1.0 / 3.0)
    return (
        1.5 * k0 * (eta**7 - eta**5) * (1.0 + 0.75 * (k0_prime - 4.0) * (eta**2 - 1.0))
    )


def reproduce() -> dict[str, object]:
    """Return source-table values and deterministic curve checkpoints."""
    ratios = np.array([1.0, 0.9, 0.8, 0.7])
    records = {}
    for mg_fraction, (v0_bohr3, hs_k0, ls_k0, transition) in TABLE_I.items():
        v0 = conventional_cell_volume(v0_bohr3)
        records[f"mg_{mg_fraction:.3f}"] = {
            "V0_a3_conventional_cell": v0,
            "K0_hs_gpa": hs_k0,
            "K0_ls_gpa_unexecutable_without_ls_v0": ls_k0,
            "K0_prime_fixed": 4.1,
            "transition_pressure_gpa": transition,
            "volume_ratios": ratios.tolist(),
            "hs_pressures_gpa": bm3_pressure(ratios * v0, v0, hs_k0, 4.1).tolist(),
        }
    return {
        "accepted_high_spin_records": 8,
        "rejected_incomplete_low_spin_rows": 8,
        "records": records,
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
