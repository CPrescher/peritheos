"""Reproduce the Dewaele et al. (2000) MgO Table 3 EOS checks."""

from __future__ import annotations

from peritheos import get_material_document
from peritheos.materials import Material

DOI = "10.1029/1999JB900364"
RATIO = 0.667


def main() -> None:
    document = get_material_document("mgo")
    source_records = [
        record
        for record in document["eos_records"]
        if record["reference"].get("doi", "").lower() == DOI.lower()
    ]
    identifiers = [record["identifier"] for record in source_records]
    material = Material.from_eosmat(document, record_identifiers=identifiers)
    by_id = {record.identifier: record for record in material.eos_records}

    for source_record in source_records:
        identifier = source_record["identifier"]
        volume = source_record["eos"]["parameters"]["V0"] * RATIO
        pressure = by_id[identifier].pressure(volume, 300.0, check_validity=False)
        print(f"{identifier}: P(V/V0={RATIO})={pressure:.8f} GPa")


if __name__ == "__main__":
    main()
