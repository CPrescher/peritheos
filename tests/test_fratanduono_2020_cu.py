"""Published Cu reference state, schema, exports and same-run reuse."""

import copy

import jsonschema
import numpy as np
import pytest

from peritheos import (
    Material,
    eosmat_schema,
    get_eos_record,
    get_eos_record_document,
    get_material_document,
    validate_pressure_calibration_references,
)
from peritheos.eos.rt import Vinet3
from scripts.validate_primary_eos_refits import (
    SHEN_CU_RECORD,
    SHEN_SELECTIONS,
    _shen_series,
    shen_cu_pressure_gpa,
)


def test_primary_record_and_material_roundtrip():
    document = get_material_document("copper")
    jsonschema.validate(document, eosmat_schema())
    record = get_eos_record_document(SHEN_CU_RECORD)
    assert record["eos"]["parameters"] == dict(
        V0=pytest.approx(4 * 63.546 / (8.939 * 0.602214076)),
        K0=133.6,
        eta=6.29,
        beta=2.06,
        psi=1.65,
    )
    assert record["parameter_errors"] == dict(
        V0=None, K0=0.8, eta=0.8, beta=0.4, psi=0.6
    )
    assert record["fixed_parameters"] == ["V0"]
    assert record["parameter_error_confidence"] is None
    assert record["parameter_covariance"] is None
    assert record["temperature_ref"] == 298
    assert record.get("experimental_pressure_range_gpa") is None
    assert "fit_datasets" not in record
    check = record["scientific_validation"]
    assert check["status"] == "primary_source_validated"
    assert check["primary_data_check"]["dataset_identifiers"] == []
    assert "reduced-isentrope" in check["primary_data_check"]["finding"]
    material = Material.from_dict(document)
    exported = material.to_dict()
    jsonschema.validate(exported, eosmat_schema())
    recovered = Material.from_dict(exported).get_eos_record(SHEN_CU_RECORD)
    assert isinstance(recovered.eos, Vinet3)
    assert recovered.eos.parameter_values() == record["eos"]["parameters"]
    assert recovered.pressure(45.94823) == pytest.approx(3.908992349216)
    assert recovered.volume(3.908992349216) == pytest.approx(45.94823)


@pytest.mark.parametrize("change", ["missing_psi", "ordinary_model", "extra_kprime"])
def test_schema_rejects_incomplete_or_mislabeled_model(change):
    document = copy.deepcopy(get_material_document("copper"))
    eos = next(
        r["eos"] for r in document["eos_records"] if r["identifier"] == SHEN_CU_RECORD
    )
    if change == "missing_psi":
        del eos["parameters"]["psi"]
    elif change == "ordinary_model":
        eos["model"] = "vinet"
    else:
        eos["parameters"]["K0_prime"] = 5.1933333333
    with pytest.raises(jsonschema.ValidationError):
        jsonschema.validate(document, eosmat_schema())


def test_shen_uses_catalog_record_and_preserves_pairing(monkeypatch):
    validate_pressure_calibration_references()
    eos = get_eos_record(SHEN_CU_RECORD).eos
    volumes = np.array([[45.94823, 40], [30, eos.V0]])
    np.testing.assert_array_equal(shen_cu_pressure_gpa(volumes), eos.pressure(volumes))
    for identifier in SHEN_SELECTIONS:
        record = get_eos_record_document(identifier)
        assert (
            record["pressure_calibration"]["methods"][0]["reference_eos_record"]
            == SHEN_CU_RECORD
        )
    # A tiny same-run fixture exercises pairing, Pt first/last preservation,
    # and exclusion of the wrong experiment independently of the saved ledger.
    rows = [
        dict(experiment=exp, run_number=run, phase=phase, unit_cell_volume_a3=str(v))
        for exp, run, phase, v in [
            ("DAC-2", "2", "Pt", 55),
            ("DAC-2", "1", "Cu", 45.94823),
            ("DAC-2", "2", "Cu", 40),
            ("DAC-2", "1", "Pt", 59),
            ("DAC-2", "1", "Pt", 58),
            ("DAC-1", "1", "Pt", 60),
        ]
    ]
    monkeypatch.setattr(
        "scripts.validate_primary_eos_refits._load_rows", lambda _: rows
    )
    result = _shen_series(
        {"identifier": "platinum_shen_2026_vinet_2"}, {"identifier": "fixture"}
    )
    np.testing.assert_array_equal(result.volume, [55, 59, 58])
    np.testing.assert_allclose(
        result.pressure, eos.pressure(np.array([40, 45.94823, 45.94823]))
    )
    assert result.pressure_sigma is None
    assert result.volume_sigma is None
