import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np
import pytest

from peritheos import Material, get_eos_record, get_material_document
from scripts.reproduce_lv_2016_qandilite import load_data, pressure, reproduce, slope

ROOT = Path(__file__).resolve().parents[1]


def test_table1_transcription_and_crystallographic_identity():
    document = get_material_document("qandilite")
    data = load_data()
    assert len(data["pressure_gpa"]) == 18
    assert np.isnan(data["pressure_sigma_gpa"][0])
    np.testing.assert_array_equal(data["pressure_sigma_gpa"][1:], 0.1)
    assert data["pressure_gpa"][[0, 1, -1]].tolist() == [0.0001, 0.4, 14.9]
    assert data["volume_a3"][[0, 10, -1]].tolist() == [602.58, 576.97, 559.91]
    assert data["volume_sigma_a3"][[0, 15, -1]].tolist() == [0.49, 0.31, 0.34]
    assert np.max(np.abs(data["a_a"] ** 3 - data["volume_a3"])) < 0.102
    resource = document["datasets"][0]["resource"]
    assert (
        hashlib.sha256(
            (ROOT / "peritheos/data" / resource["path"]).read_bytes()
        ).hexdigest()
        == resource["sha256"]
    )
    assert document["composition"]["reported_formula"] == "Mg2.00(1)Ti1.00(1)O4"
    assert document["space_group_number"] == 227
    assert document["formula_units_per_cell"] == 8
    contents = defaultdict(float)
    for site in document["atom_sites"]:
        contents[site["element"]] += site["multiplicity"] * site["occupancy"]
    assert dict(contents) == {"Mg": 16, "Ti": 8, "O": 32}
    assert document["lattice"]["a"] == 8.44192
    assert document["atom_sites"][-1]["x"] == 0.2594
    assert document["eos_records"][0]["default"]
    assert document["source"]["citation"]["doi"] == "10.2138/am-2003-5-615"


@pytest.mark.parametrize("order", [2, 3])
def test_published_observations_equation_inversion_and_interchange(order):
    document = get_material_document("qandilite")
    record = get_eos_record(f"qandilite_lv_2016_bm{order}")
    params = [603.21, 172.0] if order == 2 else [603.10, 175.0, 3.5]
    volumes = load_data()["volume_a3"]
    np.testing.assert_allclose(
        record.pressure(volumes), pressure(volumes, params), atol=2e-12
    )
    # An independent measured high-P state, with its printed one-sigma errors.
    assert record.volume(14.9) == pytest.approx(559.91, abs=0.34)
    assert record.pressure(559.91) == pytest.approx(14.9, abs=0.18)
    assert record.eos.pressure(params[0]) == pytest.approx(0, abs=1e-12)
    assert record.eos.bulk_modulus(params[0]) == pytest.approx(params[1])
    assert record.eos.bulk_modulus(570) == pytest.approx(-570 * slope(570, params))
    assert record.within_validity(record.volume(10.0), 300.0)
    with pytest.raises(ValueError, match="outside the published calibration"):
        record.pressure(550.0, 300.0, check_validity=True)
    with pytest.raises(ValueError, match="isothermal 300 K"):
        record.pressure(580.0, 301.0, check_validity=True)
    pressures = np.linspace(0.0001, 14.9, 19)
    np.testing.assert_allclose(
        record.pressure(record.volume(pressures)), pressures, atol=1e-9
    )
    restored = Material.from_eosmat(
        json.loads(json.dumps(document)), record_identifiers=[record.identifier]
    )
    assert restored.eos_records[0].pressure(570) == pytest.approx(record.pressure(570))
    assert (
        record.pressure_calibration["methods"][0]["reference_calibration_record"]
        == "ruby_mao_1978"
    )


def test_independent_refits_preserve_published_parameters():
    results = reproduce()
    saved = json.loads((ROOT / "docs/data/lv-2016-qandilite.json").read_text())[
        "results"
    ]
    for order in [2, 3]:
        key = f"bm{order}"
        result = results[key]
        record = get_material_document("qandilite")["eos_records"][order - 2]
        assert result["implementation_max_error_gpa"] < 2e-12
        assert result["published_pressure_rmse_gpa"] < 0.13
        for name, value in result["all_rows_unweighted_pressure"]["parameters"].items():
            assert (
                abs(value - record["eos"]["parameters"][name])
                < record["parameter_errors"][name]
            )
        for mode in [
            "all_rows_unweighted_pressure",
            "all_rows_weighted_volume",
            "all_rows_effective_variance",
            "dac_only_unweighted_pressure",
        ]:
            assert result[mode]["parameters"] == pytest.approx(
                saved[key][mode]["parameters"], rel=1e-5
            )
        assert record["record_kind"] == "published"
