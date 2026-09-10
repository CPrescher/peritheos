#!/usr/bin/env python3
"""Reconstruct Holmes et al.'s analytical Pt scale and audit its shock evidence."""

from __future__ import annotations

import csv
import json
import math
from pathlib import Path
from typing import Any

import numpy as np

from peritheos import Material, get_material_document

ROOT = Path(__file__).resolve().parents[1]
CHECKPOINT_PATH = (
    ROOT
    / "peritheos"
    / "data"
    / "datasets"
    / "platinum-holmes-1989-equation11-checkpoints.csv"
)
SHOCK_PATH = (
    ROOT / "peritheos" / "data" / "datasets" / "platinum-holmes-1989-table3-shock.csv"
)

# Equation (11) and Table IV. P_T and eta are the coefficients printed for
# the analytical curve; B_T and B_T' in the same table are rounded summaries.
V0_BOHR3_PER_ATOM = 101.9
P_T_GPA = 798.31
ETA = 7.2119
VINET_K0_GPA = P_T_GPA / 3.0
VINET_K0_PRIME = 1.0 + ETA / 1.5

# Equation (12) uses the separately rounded Table IV B_T=266 GPa.
TABLE_IV_B_T_GPA = 266.0
TABLE_IV_ALPHA_PER_K = 0.261e-4
ALPHA_KT_GPA_PER_K = TABLE_IV_ALPHA_PER_K * TABLE_IV_B_T_GPA


def equation_11_pressure(x: float | np.ndarray) -> float | np.ndarray:
    """Evaluate the source's printed 300 K universal-EOS expression."""
    x = np.asarray(x, dtype=float)
    pressure = P_T_GPA * (1.0 - x) / x**2 * np.exp(ETA * (1.0 - x))
    return float(pressure) if pressure.ndim == 0 else pressure


def _rows(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as stream:
        return list(csv.DictReader(stream))


def reproduce() -> dict[str, Any]:
    """Return exact curve-parity and independent shock-row diagnostics."""
    document = get_material_document("platinum")
    identifier = "platinum_holmes_1989_vinet_1"
    source = next(
        record
        for record in document["eos_records"]
        if record["identifier"] == identifier
    )
    eos = Material.from_eosmat(
        document, record_identifiers=[identifier]
    ).get_eos_record(identifier)

    checkpoints = _rows(CHECKPOINT_PATH)
    x = np.asarray([float(row["compression_coordinate_x"]) for row in checkpoints])
    volumes = np.asarray([float(row["unit_cell_volume_a3"]) for row in checkpoints])
    tabulated = np.asarray([float(row["pressure_300k_gpa"]) for row in checkpoints])
    independent = np.asarray(equation_11_pressure(x), dtype=float)
    executable = np.asarray(eos.pressure(volumes, 300.0), dtype=float)

    shock_rows = _rows(SHOCK_PATH)
    momentum_residuals = []
    mass_residuals = []
    for row in shock_rows:
        rho0 = float(row["sample_density_g_cm3"])
        shock_velocity = float(row["shock_velocity_km_s"])
        mass_velocity = float(row["mass_velocity_km_s"])
        pressure = rho0 * shock_velocity * mass_velocity
        density = rho0 / (1.0 - mass_velocity / shock_velocity)
        momentum_residuals.append(pressure - float(row["pressure_gpa"]))
        mass_residuals.append(density - float(row["density_g_cm3"]))

    parameters = source["eos"]["parameters"]
    assert math.isclose(parameters["K0"], VINET_K0_GPA, abs_tol=1e-14)
    assert math.isclose(parameters["K0_prime"], VINET_K0_PRIME, abs_tol=1e-14)
    assert source["thermal"]["parameters"]["alpha_KT"] == ALPHA_KT_GPA_PER_K

    return {
        "source_parameterization": {
            "V0_bohr3_per_atom": V0_BOHR3_PER_ATOM,
            "P_T_gpa": P_T_GPA,
            "eta": ETA,
            "equivalent_vinet_K0_gpa": VINET_K0_GPA,
            "equivalent_vinet_K0_prime": VINET_K0_PRIME,
            "table_iv_rounded_B_T_gpa": TABLE_IV_B_T_GPA,
            "table_iv_rounded_B_T_prime": 5.81,
        },
        "equation_11_checkpoints": {
            "rows": len(checkpoints),
            "pressure_range_gpa": [float(tabulated.min()), float(tabulated.max())],
            "csv_vs_equation_max_abs_gpa": float(
                np.max(np.abs(tabulated - independent))
            ),
            "peritheos_vs_equation_max_abs_gpa": float(
                np.max(np.abs(executable - independent))
            ),
        },
        "equation_12": {
            "alpha_KT_gpa_per_k": ALPHA_KT_GPA_PER_K,
            "thermal_increment_300_to_2000k_gpa": ALPHA_KT_GPA_PER_K * 1700.0,
        },
        "table_iii_shock_qualification": {
            "rows": len(shock_rows),
            "pressure_range_gpa": [
                min(float(row["pressure_gpa"]) for row in shock_rows),
                max(float(row["pressure_gpa"]) for row in shock_rows),
            ],
            "rankine_hugoniot_momentum_max_abs_rounding_gpa": float(
                np.max(np.abs(momentum_residuals))
            ),
            "rankine_hugoniot_mass_max_abs_rounding_g_cm3": float(
                np.max(np.abs(mass_residuals))
            ),
            "fit_role": "qualification_only_not_300k_isotherm_observations",
        },
        "reconstructability": {
            "published_300k_curve": "exact_analytical_reconstruction",
            "underlying_lmto_fit": "not_refittable_no_numerical_e0_or_p0_grid",
            "theoretical_hugoniot": "not_reconstructable_no_numerical_thermal_grid",
        },
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
