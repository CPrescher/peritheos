import copy
import csv
import hashlib
import json

import numpy as np
import pytest

from peritheos import (
    Material,
    convert_eosmat,
    create_eosmat,
    get_eos_record,
    get_material,
    get_material_document,
    load_eosmat,
    validate_eosmat_document,
)
from peritheos.eos.rt import BM2, BM2Reference
from peritheos.eos.thermal import ThermalReferenceStateEOS
from peritheos.errors import EosValidationError
from peritheos.fitting import fit_rt_eos
from scripts import reproduce_bykova_2018_experimental as audit


def test_reference_identity_independent_pressure_and_inversion():
    eos = BM2Reference(375.04, 254, 32.7)
    assert eos.pressure(375.04) == 32.7
    assert eos.bulk_modulus(375.04) == 254
    assert eos.pressure(375.04 / 8) == pytest.approx(36608.7, abs=1e-8)
    assert eos.bulk_modulus_derivative(375.04) == pytest.approx(4, abs=1e-7)
    pressures = np.array([[20, 32.7], [50, 70.4]])
    np.testing.assert_allclose(
        eos.pressure(eos.volume(pressures)), pressures, atol=1e-8
    )
    assert eos.volume(32.7) == 375.04
    volumes = np.array([330, 375.04, 390])
    np.testing.assert_array_equal(
        BM2Reference(375.04, 254, 0).pressure(volumes),
        BM2(375.04, 254).pressure(volumes),
    )
    np.testing.assert_array_equal(
        eos.bulk_modulus(volumes), BM2(375.04, 254).bulk_modulus(volumes)
    )


@pytest.mark.parametrize(
    "parameters",
    [(0, 254, 32.7), (375, -1, 32.7), (375, 254, np.nan), (np.inf, 254, 32.7)],
)
def test_invalid_reference_parameters(parameters):
    with pytest.raises(EosValidationError):
        BM2Reference(*parameters)


def test_native_fitting_and_thermal_reference_change_preserve_offset():
    volumes = np.linspace(335, 385, 20)
    eta = (375.04 / volumes) ** (1 / 3)
    pressures = 32.7 + 1.5 * 254 * (eta**7 - eta**5)
    result = fit_rt_eos(
        BM2Reference,
        volumes,
        pressures,
        initial={"V0": 374, "K0": 250},
        fixed={"P0": 32.7},
    )
    assert result.success
    assert result.model.P0 == 32.7
    assert result.parameters["V0"] == pytest.approx(375.04, abs=1e-6)
    assert result.parameters["K0"] == pytest.approx(254, abs=1e-5)
    thermal = ThermalReferenceStateEOS(
        BM2Reference(375.04, 254, 32.7), Tr=298, alpha0=1e-5, dK_dT=-0.02
    )
    shifted_volume = 375.04 * np.exp(1e-5 * (1000 - 298))
    assert thermal.pressure(shifted_volume, 1000) == pytest.approx(32.7, abs=1e-10)


def test_eosmat_native_and_v3_round_trip_and_standalone_example():
    doc = create_eosmat(
        "reference_demo", "SiO2", "BM2Reference", audit.PARAMETERS, temperature_ref=298
    )
    for d in [
        doc,
        convert_eosmat(doc, version=3),
        convert_eosmat(convert_eosmat(doc, version=3), version=4),
    ]:
        validate_eosmat_document(d)
        material = Material.from_eosmat(d, require_primary_validation=False)
        assert material.eos_records[0].pressure(375.04) == 32.7
        assert material.eos_records[0].eos.parameter_values() == audit.PARAMETERS
    bad = copy.deepcopy(doc)
    del bad["eos_records"][0]["parameters"]["P0"]
    with pytest.raises(ValueError):
        Material.from_eosmat(bad, require_primary_validation=False)
    example = load_eosmat(audit.ROOT / "docs/examples/bm2-reference.eosmat")
    assert (
        Material.from_eosmat(example, require_primary_validation=False)
        .eos_records[0]
        .pressure(375.04)
        == 32.7
    )


def test_published_record_preserves_coefficients_scope_data_and_errors():
    document = get_material_document("coesite_iv")
    record = next(
        r for r in document["eos_records"] if r["identifier"] == audit.RECORD_ID
    )
    assert record["eos"]["parameters"] == audit.PARAMETERS
    assert record["record_kind"] == "published"
    assert record["phase_scope"] == ["coesite_iv", "coesite_v"]
    assert record["parameter_errors"] == {"V0": 0.48, "K0": 9, "P0": None}
    assert record["parameter_error_confidence"] is None
    assert record["scientific_validation"]["usage_recommendation"] == "not_recommended"
    assert (
        get_material("coesite_iv").default_equilibrium_record().identifier
        != audit.RECORD_ID
    )
    assert get_eos_record(audit.RECORD_ID).pressure(375.04) == 32.7
    assert (
        hashlib.sha256(audit.TABLE.read_bytes()).hexdigest()
        == "b750673a49ef8c2d56241bd810a4d8f0b0ada801c65df850573ae579089f57de"
    )
    with audit.TABLE.open(newline="") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 16
    assert rows[12]["pressure_token"] == "44.5(5)"
    assert rows[12]["volume_token"] == "346.4(8)"
    assert rows[0]["sample_id"] == "DAC#2 c1"
    assert rows[11]["reported_volume_a3_per_formula_unit"] == "22.022"
    material_data = get_material("coesite_iv").get_dataset(audit.TABLE_ID)
    assert len(material_data.values("pressure_gpa")) == 16


def test_vector_curve_and_observations_have_separate_validation_claims():
    result = audit.reproduce()
    assert result["curve_check"]["observations"] == 23
    assert result["curve_check"]["all_within_graphical_bounds"] is True
    assert result["curve_check"]["max_absolute_residual_gpa"] == pytest.approx(
        0.04263156, abs=1e-7
    )
    obs = result["observations_check"]
    assert obs["observations"] == 16
    assert obs["coesite_iv_rows"] == 12
    assert obs["rmse_gpa"] == pytest.approx(4.04658696, abs=1e-7)
    assert obs["max_absolute_residual_gpa"] == pytest.approx(11.85456532, abs=1e-7)
    assert obs["largest_residual_source_row"] == 13
    assert result["author_workflow_reproduced"] is False
    assert result["diagnostic_refit"]["adopted_as_replacement"] is False
    from scripts.bykova_replay_checks import check_report

    check_report(json.loads(audit.RESULT.read_text()), result)
