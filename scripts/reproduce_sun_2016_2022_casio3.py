#!/usr/bin/env python3
"""Reproduce the Sun 2016 cubic and Sun 2022 tetragonal CaSiO3 fits.

The inputs are direct factual transcriptions of the two papers' Table 1 data.
This script deliberately has no dependency on the removed Sun et al. (2010)
first-principles benchmark.
"""

from __future__ import annotations

import csv
from pathlib import Path

import numpy as np

from peritheos import Material, get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import fit_joint_eos, fit_rt_eos

ROOT = Path(__file__).resolve().parents[1]
DATA = ROOT / "peritheos" / "data" / "datasets"

SUN_2016_DATA = DATA / "ca-perovskite-sun-2016-table1-pvt.csv"
SUN_2022_DATA = DATA / "ca-perovskite-tetragonal-sun-2022-table1-pv.csv"


def _columns(path: Path, *names: str) -> dict[str, np.ndarray]:
    with path.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    return {
        name: np.asarray([float(row[name]) for row in rows], dtype=float)
        for name in names
    }


def fit_sun_2016(*, weighted: bool = False):
    """Fit all 144 cubic Table 1 P-V-T observations to source model 1."""
    values = _columns(
        SUN_2016_DATA,
        "pressure_gpa",
        "pressure_uncertainty_gpa",
        "temperature_k",
        "volume_a3_per_formula_unit",
        "volume_uncertainty_a3",
    )
    document = get_material_document("ca_perovskite")
    executable = Material.from_eosmat(
        document, record_identifiers=["ca_perovskite_sun_2016_bm3_3"]
    ).eos_records[0]
    scale = executable.volume_scale
    return fit_joint_eos(
        MieGruneisenDebye,
        BM3,
        volume=values["volume_a3_per_formula_unit"] * scale,
        temperature=values["temperature_k"],
        pressure=values["pressure_gpa"],
        initial={
            "rt_eos.V0": 45.4 * scale,
            "rt_eos.K0": 249.0,
            "gamma0": 1.8,
            "q": 1.1,
        },
        fixed={
            "rt_eos.K0_prime": 4.0,
            "Tr": 300.0,
            "theta0": 1000.0,
            "n": 5.0,
        },
        configuration={"debye_temperature_law": "integrated_gruneisen"},
        bounds={
            "rt_eos.V0": (20.0 * scale, 80.0 * scale),
            "rt_eos.K0": (10.0, 1000.0),
            "gamma0": (1.0e-9, 10.0),
            "q": (1.0e-9, 10.0),
        },
        pressure_sigma=(values["pressure_uncertainty_gpa"] if weighted else None),
        volume_sigma=(values["volume_uncertainty_a3"] * scale if weighted else None),
        absolute_sigma=weighted,
        max_nfev=5000,
    )


def fit_sun_2022(*, fixed_k0_prime: bool, weighted: bool = False):
    """Fit all 23 tetragonal Table 1 P-V observations."""
    values = _columns(
        SUN_2022_DATA,
        "pressure_gpa",
        "pressure_uncertainty_gpa",
        "volume_a3_conventional_cell",
        "volume_uncertainty_a3",
    )
    initial = {"V0": 182.4, "K0": 229.0}
    fixed = {"K0_prime": 4.0}
    if not fixed_k0_prime:
        initial["K0_prime"] = 4.0
        fixed = {}
    return fit_rt_eos(
        BM3,
        volume=values["volume_a3_conventional_cell"],
        pressure=values["pressure_gpa"],
        initial=initial,
        fixed=fixed,
        bounds={"V0": (100.0, 260.0), "K0": (10.0, 1000.0), "K0_prime": (0.0, 20.0)},
        pressure_sigma=(values["pressure_uncertainty_gpa"] if weighted else None),
        volume_sigma=(values["volume_uncertainty_a3"] if weighted else None),
        absolute_sigma=weighted,
        max_nfev=5000,
    )


def _rmse(result) -> float:
    return float(np.sqrt(np.mean(np.asarray(result.residuals, dtype=float) ** 2)))


def main() -> None:
    cubic = fit_sun_2016()
    cubic_weighted = fit_sun_2016(weighted=True)
    tetragonal_fixed = fit_sun_2022(fixed_k0_prime=True)
    tetragonal_free = fit_sun_2022(fixed_k0_prime=False)
    tetragonal_weighted = fit_sun_2022(fixed_k0_prime=False, weighted=True)
    scale = (
        Material.from_eosmat(
            get_material_document("ca_perovskite"),
            record_identifiers=["ca_perovskite_sun_2016_bm3_3"],
        )
        .eos_records[0]
        .volume_scale
    )

    print("Sun et al. (2016), 144-row cubic BM3-MGD, unweighted P residuals")
    for name in cubic.free_parameters:
        value = cubic.parameters[name]
        error = cubic.standard_errors[name]
        if name == "rt_eos.V0":
            value /= scale
            error /= scale
        print(f"  {name}: {value:.12g} +/- {error:.12g}")
    print(f"  pressure RMSE: {_rmse(cubic):.12g} GPa")
    print(
        "  weighted discriminator:",
        {
            name: cubic_weighted.parameters[name]
            / (scale if name == "rt_eos.V0" else 1.0)
            for name in cubic_weighted.free_parameters
        },
    )

    for label, result in (
        ("fixed K0'=4", tetragonal_fixed),
        ("free K0'", tetragonal_free),
    ):
        print(f"Sun et al. (2022), 23-row tetragonal BM3, {label}")
        for name in result.free_parameters:
            value = result.parameters[name] / (4.0 if name == "V0" else 1.0)
            error = result.standard_errors[name] / (4.0 if name == "V0" else 1.0)
            print(f"  {name}: {value:.12g} +/- {error:.12g}")
        print(f"  pressure RMSE: {_rmse(result):.12g} GPa")
    print(
        "  weighted free-fit discriminator:",
        {
            name: tetragonal_weighted.parameters[name] / (4.0 if name == "V0" else 1.0)
            for name in tetragonal_weighted.free_parameters
        },
    )


if __name__ == "__main__":
    main()
