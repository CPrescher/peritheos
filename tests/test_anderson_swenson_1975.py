"""Independent published coefficients and source-table reproduction checks."""

import json

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.eos.rt import OddInversePower
from peritheos.errors import EosValidationError
from scripts.reproduce_anderson_swenson_1975 import FACTOR, REPORT, reproduce, rows


def test_primary_source_diagnostics(assert_audit_close):
    result = reproduce()
    assert_audit_close(result, json.loads(REPORT.read_text()))
    assert (result["raw_table_cells"], result["usable_observations"]) == (394, 379)
    assert result["independent_refit_status"] == "not_reproduced"
    for check in result["isotherms"]:
        assert check["python_max_difference_gpa"] < 4e-14
        assert check["inverse_max_difference_cm3_mol"] < 1e-9
        assert (
            abs(
                check["calculated_bulk_at_zero_kbar"]
                - check["published_bulk_at_zero_kbar"]
            )
            < 0.051
        )
    assert result["eq9_vs_table1_max_difference_gpa"] < 0.0006
    # The data do not warrant an exact raw-data-fit reproduction claim.
    assert result["isotherms"][0]["literal_table11_corrected_volume_rmse_cm3_mol"] > 0.1


def test_source_sentinels_temperatures_and_normalization():
    data = rows("argon-anderson-1972-appendix-a.csv")
    assert sum(r["usable"] == "0" for r in data) == 15
    for r in data:
        if r["usable"] == "0":
            assert float(r["reported_volume_cm3_mol"]) == 0
            assert r["volume_a3"] == ""
        else:
            assert float(r["volume_a3"]) == pytest.approx(
                float(r["reported_volume_cm3_mol"]) * FACTOR
            )
    assert any(float(r["temperature_k"]) == 67.9 for r in data)
    assert max(float(r["pressure_kbar"]) for r in data) == 20


def test_analytic_bulk_modulus_and_material_roundtrip():
    document = get_material_document("argon_fcc")
    material = Material.from_eosmat(document)
    restored = Material.from_eosmat(material.to_eosmat())
    for record in material.eos_records:
        if "anderson_swenson" not in record.identifier:
            continue
        model = record.eos
        assert isinstance(model, OddInversePower)
        v = model.V0 * np.array([0.8, 0.9, 1.0])
        h = v * 1e-5
        numerical = -v * (model.pressure(v + h) - model.pressure(v - h)) / (2 * h)
        assert model.bulk_modulus(v) == pytest.approx(numerical, rel=2e-8)
        assert restored.get_eos_record(record.identifier).pressure(v) == pytest.approx(
            record.pressure(v)
        )
        assert record.reference_temperature in [4.2, 20, 40, 60, 77]
        assert not record.eosmat_metadata.get("fit_datasets")
        assert record.eosmat_metadata["comparison_datasets"] == [
            "argon_anderson_1972_appendix_a"
        ]


@pytest.mark.parametrize("v", [0, -1, float("nan")])
def test_invalid_volume_rejected(v):
    with pytest.raises(EosValidationError):
        OddInversePower(10, -1, 2, -0.5, 0).pressure(v)
