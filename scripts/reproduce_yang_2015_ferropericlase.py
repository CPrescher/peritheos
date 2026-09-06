#!/usr/bin/env python3
"""Execute the defensible Yang et al. (2015) high-spin BM3 branch."""

from __future__ import annotations

import json
from math import isclose

from peritheos import Material, get_material_document

MATERIAL = "mg092fe008o"
RECORD = "mg092fe008o_yang_2015_hs_bm3_reference"


def reproduce() -> dict[str, object]:
    document = get_material_document(MATERIAL)
    record = next(row for row in document["eos_records"] if row["identifier"] == RECORD)
    expected_v0 = 4.1996**3
    assert isclose(record["eos"]["parameters"]["V0"], expected_v0, rel_tol=1e-14)
    assert record["eos"]["parameters"] == {
        "V0": 74.06683401593601,
        "K0": 152.5,
        "K0_prime": 4.1,
    }
    eos = Material.from_eosmat(document, record_identifiers=[RECORD]).get_eos_record(
        RECORD
    )
    checkpoints = {}
    for fraction in (0.98, 0.90, 0.82):
        target_volume = fraction * expected_v0
        pressure = float(eos.pressure(target_volume))
        recovered = float(eos.volume(pressure))
        assert isclose(recovered, target_volume, rel_tol=2e-10)
        checkpoints[str(fraction)] = {
            "volume_a3": target_volume,
            "pressure_gpa": pressure,
            "recovered_volume_a3": recovered,
        }
    return {"record": RECORD, "checkpoints": checkpoints}


def main() -> None:
    print(json.dumps(reproduce(), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
