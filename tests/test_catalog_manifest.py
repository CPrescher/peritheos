"""The manifest counts the catalog exactly; other tests take their totals from it."""

from collections import Counter

from peritheos import get_material_document, list_material_documents
from tests.catalog_state import manifest


def test_manifest_counts_every_material_record_and_status():
    documents = [get_material_document(m) for m in list_material_documents()]
    records = [r for d in documents for r in d["eos_records"]]
    data = manifest()
    assert data["materials"] == len(documents)
    assert data["eos_records"] == len(records)
    assert len({r["identifier"] for r in records}) == len(records)
    assert data["scientific_validation"]["counts"] == dict(
        sorted(Counter(r["scientific_validation"]["status"] for r in records).items())
    )
    calibration = data["pressure_calibration"]
    assert calibration["status_counts"] == {
        status: sum(r["pressure_calibration"]["status"] == status for r in records)
        for status in calibration["status_counts"]
    }
    assert sum(calibration["status_counts"].values()) == len(records)
    assert calibration["recalculation_counts"] == dict(
        sorted(
            Counter(
                r["pressure_calibration"]["recalculation"]["status"] for r in records
            ).items()
        )
    )
