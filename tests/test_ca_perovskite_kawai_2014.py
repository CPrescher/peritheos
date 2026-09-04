import csv
import io
from importlib import resources

import numpy as np
import pytest

from peritheos import get_material_document
from peritheos.materials import Material

IDENTIFIER = "ca_perovskite_kawai_2014_vinet_mgd_3"


def _source_and_record():
    document = get_material_document("ca_perovskite")
    source = next(
        record for record in document["eos_records"] if record["identifier"] == IDENTIFIER
    )
    record = Material.from_eosmat(
        document, record_identifiers=[IDENTIFIER]
    ).eos_records[0]
    return document, source, record


def test_kawai_2014_vinet_mgd_parameters_and_reference_state():
    _, source, record = _source_and_record()

    assert source["eos"]["type"] == "Vinet"
    assert source["eos"]["parameters"] == {
        "V0": pytest.approx(46.17),
        "K0": pytest.approx(203.95),
        "K0_prime": pytest.approx(4.76),
    }
    assert source["parameter_errors"] == {
        "V0": None,
        "K0": None,
        "K0_prime": None,
    }
    assert source["thermal"]["parameters"] == {
        "Tr": pytest.approx(1000.0),
        "theta0": pytest.approx(1100.0),
        "gamma0": pytest.approx(1.576),
        "q": pytest.approx(0.96),
        "n": pytest.approx(5.0),
    }
    assert source["thermal"]["debye_temperature_law"] == "integrated_gruneisen"
    assert source["thermal"]["fixed_parameters"] == ["Tr", "theta0", "n"]
    assert record.reference_volume == pytest.approx(46.17)
    assert record.pressure(46.17, 1000.0) == pytest.approx(0.0, abs=1.0e-12)
    assert record.thermal_pressure_increment(46.17, 1000.0) == pytest.approx(
        0.0, abs=1.0e-12
    )


def test_kawai_2014_table1_isochors_and_round_trip():
    document, source, record = _source_and_record()
    dataset = next(
        item
        for item in document["datasets"]
        if item["identifier"] == "ca_perovskite_kawai_2014_table1_isochors"
    )
    payload = (
        resources.files("peritheos.data")
        .joinpath(dataset["resource"]["path"])
        .read_text(encoding="utf-8")
    )
    rows = list(csv.DictReader(io.StringIO(payload)))

    assert len(rows) == 60
    assert rows[0] == {
        "compression_1_minus_v_over_v0": "0.00",
        "volume_ratio": "1.00",
        "volume_a3_per_formula_unit": "46.1700",
        "temperature_k": "1000",
        "pressure_gpa": "0.0",
    }
    assert rows[-1]["pressure_gpa"] == "184.5"
    residuals = [
        record.pressure(
            float(row["volume_a3_per_formula_unit"]),
            float(row["temperature_k"]),
        )
        - float(row["pressure_gpa"])
        for row in rows
    ]
    assert np.sqrt(np.mean(np.square(residuals))) == pytest.approx(
        0.0373438719, abs=2.0e-8
    )
    assert np.max(np.abs(residuals)) < 0.1

    volume = record.volume(100.0, 2500.0, check_validity=True)
    assert record.pressure(volume, 2500.0, check_validity=True) == pytest.approx(
        100.0
    )
    assert source["fit_datasets"] == []
    assert (
        source["scientific_validation"]["primary_data_check"]["status"]
        == "theoretical_parameterization_only"
    )
