#!/usr/bin/env python3
"""Reproduce the shipped Noguchi et al. (2013) CaSiO3 model-1 checks.

The complete Table 1 transcription used for the independent refit is not
redistributed because doi:10.1007/s00269-012-0549-1 states no reusable data
license. Refit results and source-row selection are preserved in the material
record and docs/literature-reproductions.md.
"""

from __future__ import annotations

from peritheos import get_material_document
from peritheos.materials import Material

RECORD_ID = "ca_perovskite_noguchi_2013_bm2_mgd_1"


def isothermal_bulk_modulus(record, volume: float, temperature: float) -> float:
    """Evaluate -V(dP/dV)_T in the record's public cell-volume unit."""
    step = volume * 1.0e-6
    derivative = (
        record.pressure(volume + step, temperature)
        - record.pressure(volume - step, temperature)
    ) / (2.0 * step)
    return -volume * derivative


def main() -> None:
    document = get_material_document("ca_perovskite")
    source = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    record = Material.from_eosmat(document, record_identifiers=[RECORD_ID]).eos_records[
        0
    ]

    reference_pressure = record.pressure(46.5, 700.0)
    extrapolated_volume = record.volume(0.0, 300.0)
    extrapolated_modulus = isothermal_bulk_modulus(record, extrapolated_volume, 300.0)

    print(f"record={RECORD_ID}")
    print(f"reference P(46.5 A^3, 700 K)={reference_pressure:.12g} GPa")
    print(
        "300 K zero-pressure extrapolation: "
        f"V={extrapolated_volume:.8f} A^3, K={extrapolated_modulus:.8f} GPa "
        "(paper: 45.8 A^3, 225 GPa)"
    )
    print("Figure 4 model-1 isotherm checks at V=34 A^3:")
    for temperature in (700.0, 1600.0, 2100.0):
        print(
            f"  T={temperature:.0f} K: P={record.pressure(34.0, temperature):.8f} GPa"
        )

    refit = source["scientific_validation"]["independent_refit"]
    print(f"independent refit status={refit['result']}")
    print(refit["finding"])


if __name__ == "__main__":
    main()
