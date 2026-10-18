"""Native reading order preserves complete scientific content and legacy layout."""

import copy
import json
from pathlib import Path

from peritheos import convert_eosmat, load_eosmat, save_eosmat
from peritheos._eosmat_format import format_eosmat

ROOT = Path(__file__).resolve().parents[1]


def canonical(document):
    # Retain number spelling, signed zero, nulls and every ordered array.
    return json.dumps(document, sort_keys=True, ensure_ascii=False, allow_nan=False)


def test_native_practical_information_precedes_long_citations_and_evidence(tmp_path):
    document = load_eosmat(ROOT / "peritheos/data/materials/titanium_alpha.eosmat")
    record = document["eos_records"][0]
    # Exercise late-added extensions, not just fields in already-formatted files.
    record["source_specific_limit"] = {"notes": "Synthetic scope qualification."}
    record["extension"] = {
        "thermal": {"parameters": {"b": 2, "a": 1}, "model": "literal metadata"},
        "evidence": {"z": 0, "a": -0.0},
        "ordered": [None, -0.0, 1.2345678901234567, 1e-30, 10**40],
    }
    document["eos_records"][0] = dict(reversed(list(record.items())))
    before = copy.deepcopy(document)
    path = tmp_path / "native.eosmat"
    save_eosmat(path, document)
    result = json.loads(path.read_text())
    keys = list(result["eos_records"][0])
    assert keys[:5] == ["identifier", "label", "model", "temperature_ref", "parameters"]
    for key in (
        "parameter_error_confidence",
        "experimental_pressure_range_gpa",
        "experimental_temperature_range_k",
        "notes",
        "source_specific_limit",
    ):
        assert keys.index(key) < keys.index("reference") < keys.index("evidence")
    assert keys[-1] == "evidence"
    extension = result["eos_records"][0]["extension"]
    assert list(extension["thermal"]) == ["parameters", "model"]
    assert list(extension["thermal"]["parameters"]) == ["b", "a"]
    assert list(extension["evidence"]) == ["z", "a"]
    assert canonical(before) == canonical(result) == canonical(document)
    assert canonical(convert_eosmat(before, version=3)) == canonical(
        convert_eosmat(result, version=3)
    )
    assert path.read_text() == format_eosmat(result)


def test_hold_selection_scope_and_covariance_are_visible_without_changing_meaning():
    document = load_eosmat(ROOT / "docs/examples/simple-material.eosmat")
    document["recommendation_hold"] = {"status": "on_hold", "reason": "Synthetic test."}
    record = document["eos_records"][0]
    record.update(
        catalog_access="explicit_selection",
        default=False,
        validity={"pressure_gpa": [0, 5]},
        experimental_pressure_range_gpa=[0, 10],
        parameter_covariance={
            "parameters": ["K0_prime", "K0"],
            "matrix": [[2.0, 0.25], [0.25, 3.0]],
        },
    )
    rendered = format_eosmat(document)
    result = json.loads(rendered)
    assert rendered.index('"recommendation_hold"') < rendered.index('"eos_records"')
    record = result["eos_records"][0]
    keys = list(record)
    assert keys.index("validity") < keys.index("parameter_covariance")
    assert keys.index("catalog_access") < keys.index("parameter_covariance")
    assert (
        record["experimental_pressure_range_gpa"] != record["validity"]["pressure_gpa"]
    )
    assert canonical(result) == canonical(document)


def test_every_native_thermal_component_has_model_parameters_then_laws():
    count = 0
    for path in sorted((ROOT / "peritheos/data/materials").glob("*.eosmat")):
        document = load_eosmat(path)
        before = canonical(document)
        rendered = format_eosmat(document)
        result = json.loads(rendered)
        assert canonical(result) == before, path.name
        assert canonical(document) == before, path.name
        assert format_eosmat(result) == rendered, path.name
        for original, record in zip(document["eos_records"], result["eos_records"]):
            assert list(record)[-1] == "evidence", path.name
            if "thermal" not in record:
                continue
            count += 1
            assert list(record["thermal"])[:2] == ["model", "parameters"], path.name
            assert list(record["thermal"]["parameters"]) == list(
                original["thermal"]["parameters"]
            ), path.name
    assert count > 0


def test_legacy_thermal_order_stays_unchanged():
    document = convert_eosmat(
        load_eosmat(ROOT / "docs/examples/thermal-material.eosmat"), version=3
    )
    thermal = document["eos_records"][0]["thermal"]
    thermal = dict(reversed(list(thermal.items())))
    document["eos_records"][0]["thermal"] = thermal
    result = json.loads(format_eosmat(document))
    assert list(result["eos_records"][0]["thermal"]) == list(thermal)
    assert canonical(result) == canonical(document)
