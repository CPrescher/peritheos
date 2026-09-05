#!/usr/bin/env python3
"""Reproduce the Datchi et al. (2007) cubic-BN thermal EOS.

Primary source: Table IV, Table V, and Equations (2)-(4),
doi:10.1103/PhysRevB.75.214104.
"""

from __future__ import annotations

import csv
from importlib import resources

import numpy as np

from peritheos import get_material_document
from peritheos.eos.rt import Vinet
from peritheos.eos.thermal import MieGruneisenDebye
from peritheos.fitting import FitResult, fit_thermal_eos
from peritheos.units import cell_volume_to_molar_volume

DATASET_ID = "cubic_boron_nitride_datchi_2007_table4_pvt"
CELL_TO_MOLAR = float(cell_volume_to_molar_volume(1.0, 4.0))
REFERENCE = Vinet(47.2208 * CELL_TO_MOLAR, 397.0, 3.62)
FIXED = {"Tr": 295.0, "theta0": 1700.0, "gamma0": 1.04, "n": 2.0}
CONFIGURATION = {
    "debye_temperature_law": "integrated_gruneisen",
    "thermal_pressure_reference": "absolute_zero",
}


def load_data() -> np.ndarray:
    """Load all printed Table IV rows and reconstruct cell volume as a^3."""
    document = get_material_document("boron_nitride")
    dataset = next(
        item for item in document["datasets"] if item["identifier"] == DATASET_ID
    )
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    return np.array(
        [
            (
                float(row["pressure_gpa"]),
                float(row["temperature_k"]),
                float(row["lattice_a_angstrom"]),
            )
            for row in rows
        ],
        dtype=[("pressure", float), ("temperature", float), ("a", float)],
    )


def model(q: float = 4.0) -> MieGruneisenDebye:
    """Build the source's 0 K Vinet plus absolute Debye-pressure EOS."""
    return MieGruneisenDebye(
        REFERENCE,
        q=q,
        **FIXED,
        **CONFIGURATION,
    )


def published_diagnostics(data: np.ndarray) -> tuple[float, float]:
    """Return pressure RMSE and the independently reported 300 K P=0 volume."""
    cell_volume = data["a"] ** 3
    predicted = model().pressure(
        cell_volume * CELL_TO_MOLAR,
        data["temperature"],
    )
    rmse = float(np.sqrt(np.mean((predicted - data["pressure"]) ** 2)))
    volume_300 = float(model().volume(0.0, 300.0) / CELL_TO_MOLAR / 8.0)
    return rmse, volume_300


def refit_q(data: np.ndarray) -> FitResult:
    """Fit q alone to all printed states with unweighted pressure residuals."""
    return fit_thermal_eos(
        MieGruneisenDebye,
        REFERENCE,
        data["a"] ** 3 * CELL_TO_MOLAR,
        data["temperature"],
        data["pressure"],
        initial={"q": 4.0},
        fixed=FIXED,
        configuration=CONFIGURATION,
        bounds={"q": (-20.0, 20.0)},
    )


def main() -> None:
    data = load_data()
    rmse, volume_300 = published_diagnostics(data)
    result = refit_q(data)
    print(f"Loaded {len(data)} Table IV rows")
    print(f"Published q=4 curve: pressure RMSE={rmse:.6f} GPa")
    print(
        "Published-curve P=0 volume at 300 K: "
        f"{volume_300:.7f} A^3/atom (Table VI: 5.9055)"
    )
    print(
        f"Unweighted q-only refit: q={result.parameters['q']:.6f} +/- "
        f"{result.standard_errors['q']:.6f}; "
        f"RMSE={np.sqrt(np.mean(result.residuals**2)):.6f} GPa; "
        f"dof={result.degrees_of_freedom}"
    )


if __name__ == "__main__":
    main()
