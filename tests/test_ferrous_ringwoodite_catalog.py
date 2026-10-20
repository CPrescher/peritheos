"""The accepted static Fe endmember is public; unresolved thermal work stays staged."""

import numpy as np
import pytest

from peritheos import Material, get_material_document
from peritheos.errors import MaterialError
from scripts.ferrous_ringwoodite_static import RECORD, ledger_outcome


def test_public_material_loads_static_record_and_retains_primary_observations():
    document = get_material_document("ferrous_ringwoodite")
    assert [r["identifier"] for r in document["eos_records"]] == [RECORD]
    material = Material.from_eosmat(document)
    record = material.default_record()
    assert record.identifier == RECORD
    assert record.reference_temperature == 296.15
    assert record.pressure(553.5, 296.15) == pytest.approx(1.867379040453365)
    assert record.pressure(547.7, 296.15) == pytest.approx(4.0556076738036)
    p = np.array([0.0, 2.0, 4.0])
    np.testing.assert_allclose(
        record.pressure(record.volume(p, 296.15), 296.15), p, atol=1e-8
    )
    with pytest.raises(MaterialError, match="isothermal"):
        record.volume(1, 300, check_validity=True)
    with pytest.raises(MaterialError, match="outside"):
        record.volume(5, 296.15, check_validity=True)
    dataset = material.get_dataset("ferrous_ringwoodite_finger_1979_table2")
    assert len(dataset) == 3
    assert dataset.checksum_verified
    assert document["eos_records"][0]["record_kind"] == "derived"
    assert document["eos_records"][0]["eos"]["parameters"] == {"V0": 558.7, "K0": 196.0}


def test_static_reproduction_and_plot_selection_agree():
    from scripts.reconstruction_plot_data import curve_comparison

    document = get_material_document("ferrous_ringwoodite")
    record = document["eos_records"][0]
    result = ledger_outcome(record)
    assert result["parameters"][0]["refit"] == pytest.approx(196.2175873043)
    assert result["observations"] == 2
    comparison = curve_comparison(document, record)
    assert comparison is not None
    assert result["rmse_gpa"] < result["published_rmse_gpa"]
