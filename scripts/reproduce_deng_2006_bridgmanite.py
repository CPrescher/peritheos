#!/usr/bin/env python3
"""Reproduce the complete Deng et al. (2006) third-order Birch branch."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

IDENTIFIER = "bridgmanite_deng_2006_lda_bm3"


def reproduce() -> dict[str, object]:
    document = get_material_document("bridgmanite")
    record = next(
        row for row in document["eos_records"] if row["identifier"] == IDENTIFIER
    )
    assert record["eos"]["parameters"] == {"V0": 161.658, "K0": 255.0, "K0_prime": 4.0}
    eos = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).get_eos_record(IDENTIFIER)
    pressure = float(eos.pressure(0.8 * 161.658))
    recovered = float(eos.volume(pressure))
    assert isclose(recovered, 0.8 * 161.658, rel_tol=2e-10)
    return {
        "record": IDENTIFIER,
        "pressure_at_0.8_v0_gpa": pressure,
        "roundtrip_volume_a3_cell": recovered,
    }


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
