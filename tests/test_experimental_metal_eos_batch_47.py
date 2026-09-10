import csv
from collections import Counter
from importlib import resources

import numpy as np
import pytest

from peritheos import Material, get_material_document, list_material_documents
from peritheos.eos.rt import Baonza


def _rows():
    resource = resources.files("peritheos").joinpath(
        "data", "datasets", "experimental-metal-eos-batch-47.csv"
    )
    with resource.open(encoding="utf-8", newline="") as stream:
        return list(csv.DictReader(stream))


def test_batch_has_exact_expected_provenance_split_and_unique_records():
    rows = _rows()
    assert len(rows) == 47
    assert len({row["record_identifier"] for row in rows}) == 47
    assert Counter(row["family"] for row in rows) == {
        "experimental_fit": 44,
        "thermal_mixed": 3,
    }
    assert Counter(row["reference_doi"] for row in rows) == {
        "10.3390/min9110684": 30,
        "10.1103/PhysRevB.70.094112": 8,
        "10.1029/2012JB009292": 6,
        "10.3390/cryst15030221": 3,
    }


def test_every_batch_record_is_primary_validated_and_executable():
    identifiers = {row["record_identifier"] for row in _rows()}
    found = set()
    for material_identifier in list_material_documents():
        document = get_material_document(material_identifier)
        selected = [
            record["identifier"]
            for record in document["eos_records"]
            if record["identifier"] in identifiers
        ]
        if not selected:
            continue
        material = Material.from_eosmat(document, record_identifiers=selected)
        for record in material.eos_records:
            raw = next(
                item
                for item in document["eos_records"]
                if item["identifier"] == record.identifier
            )
            v0 = raw["eos"]["parameters"]["V0"]
            assert record.pressure(v0, record.reference_temperature) == pytest.approx(
                0.0, abs=1.0e-8
            )
            found.add(record.identifier)
    assert found == identifiers


def test_baonza_reference_state_derivative_and_vectorization():
    eos = Baonza(10.0, 140.0, 4.5)
    assert eos.pressure(10.0) == pytest.approx(0.0, abs=1.0e-13)
    assert eos.bulk_modulus(10.0) == pytest.approx(140.0)
    step = 1.0e-5
    derivative = (eos.bulk_modulus(10.0 - step) - eos.bulk_modulus(10.0 + step)) / (
        eos.pressure(10.0 - step) - eos.pressure(10.0 + step)
    )
    assert derivative == pytest.approx(4.5, rel=1.0e-8)
    pressure = eos.pressure(np.array([10.0, 9.5, 9.0]))
    assert pressure.shape == (3,)
    assert np.all(np.diff(pressure) > 0.0)


def test_thermal_hcp_iron_records_add_temperature_dependence():
    document = get_material_document("iron")
    identifiers = [
        record["identifier"]
        for record in document["eos_records"]
        if "iron_zhang_2025_fit" in record["identifier"]
    ]
    assert len(identifiers) == 3
    material = Material.from_eosmat(document, record_identifiers=identifiers)
    for record in material.eos_records:
        raw = next(
            item
            for item in document["eos_records"]
            if item["identifier"] == record.identifier
        )
        v0 = raw["eos"]["parameters"]["V0"]
        assert record.pressure(v0, 300.0) == pytest.approx(0.0, abs=1.0e-9)
        assert record.pressure(v0, 2000.0) > 0.0
