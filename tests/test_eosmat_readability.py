"""Readable EOSMAT remains lossless, compatible, and honest about source status."""

import copy
import json
import math
import subprocess
import sys
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import Material, convert_eosmat, eosmat_schema, load_eosmat, save_eosmat
from peritheos._eosmat_format import format_eosmat
from peritheos.errors import MaterialError
from scripts.format_eosmat import format_file
from scripts.generate_material_summaries import (
    availability,
    coefficient,
    parameters,
    render_material,
)

ROOT = Path(__file__).resolve().parents[1]
TEMPLATE = ROOT / "docs/examples/simple-material.eosmat"


def test_template_is_valid_but_requires_explicit_unaudited_loading():
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    Draft202012Validator(eosmat_schema()).validate(document)
    with pytest.raises(MaterialError, match="pending_primary_source_check"):
        Material.from_eosmat(document)
    material = Material.from_eosmat(document, require_primary_validation=False)
    record = material.eos_records[0]
    assert float(record.pressure(10.0)) == 0.0
    assert record.scientific_validation_status == "pending_primary_source_check"


def test_readable_writer_preserves_extensions_array_order_and_numerical_precision(
    tmp_path,
):
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    document["extension"] = {
        "text": 'Å³, quotes " and backslash \\ and line\nbreak',
        "values": [None, -0.0, 1.2345678901234567, 1e-30, 10**40],
        "matrix": [[3, 2, 1], [6, 5, 4]],
        "empty_object": {},
        "empty_array": [],
    }
    sibling = copy.deepcopy(document["eos_records"][0])
    sibling["identifier"] = "another_example"
    document["eos_records"].insert(0, sibling)
    before = copy.deepcopy(document)
    path = tmp_path / "material.eosmat"
    save_eosmat(path, document)
    loaded = load_eosmat(path)
    assert loaded == before == document
    assert math.copysign(1, loaded["extension"]["values"][1]) == -1
    payload = path.read_text(encoding="utf-8")
    assert payload == format_eosmat(loaded)
    assert path.read_bytes() == payload.encode("utf-8")
    assert payload.index('"identifier": "example_material"') < payload.index(
        '"eos_records"'
    )
    assert '"parameters": {"V0": 10.0, "K0": 160.0, "K0_prime": 4.0}' in payload
    assert format_file(path) is False


def test_failed_save_does_not_truncate_an_existing_file(tmp_path):
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    document["extension"] = {"nonfinite": float("nan")}
    path = tmp_path / "existing.eosmat"
    path.write_text("previous contents", encoding="utf-8")
    with pytest.raises(ValueError):
        save_eosmat(path, document)
    assert path.read_text(encoding="utf-8") == "previous contents"


def test_format_check_is_read_only_and_rewrite_is_lossless(tmp_path):
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    path = tmp_path / "material.eosmat"
    original = json.dumps(document).encode()
    path.write_bytes(original)
    command = [sys.executable, "-m", "scripts.format_eosmat"]
    result = subprocess.run(
        command + ["--check", str(path)], cwd=ROOT, capture_output=True
    )
    assert result.returncode == 1
    assert path.read_bytes() == original
    subprocess.run(command + [str(path)], cwd=ROOT, check=True, capture_output=True)
    assert load_eosmat(path) == document
    assert (
        subprocess.run(
            command + ["--check", str(path)], cwd=ROOT, capture_output=True
        ).returncode
        == 0
    )


def test_summary_distinguishes_holds_zero_kelvin_and_restricted_records():
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    record = document["eos_records"][0]
    record["temperature_ref"] = 0.0
    record["scientific_validation"]["status"] = "not_reproduced"
    record["catalog_access"] = "explicit_selection"
    document["recommendation_hold"] = {
        "status": "on_hold",
        "scope": "equilibrium_default",
        "reason": "Synthetic test hold; no replacement default selected.",
    }
    executable = {r["identifier"] for r in document["eos_records"]}
    sheet = render_material(document, executable=executable, defaults={})
    assert "**Equilibrium:** None" in sheet
    assert "Default selection hold" in sheet
    assert "Explicit diagnostic selection; not recommended" in sheet
    assert "| 0.0 |" in sheet
    assert record["identifier"] in sheet
    for record in document["eos_records"]:
        assert f'<a id="record-{record["identifier"]}"></a>' in sheet


def test_summary_keeps_fixed_values_and_null_reference_conventions_distinct():
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    record = document["eos_records"][0]
    record["thermal"] = {
        "parameters": {"Tr": None},
        "parameter_errors": {"Tr": None},
        "fixed_parameters": [],
    }
    rendered = "\n".join(parameters(record))
    assert "4.0 | Not reported | Yes" in rendered
    assert "null (model convention)" in rendered
    assert "standard deviation" not in rendered
    assert coefficient(record, "K0_prime") == "4.0 (fixed)"
    record["parameter_errors"]["K0"] = 0.5
    assert coefficient(record, "K0") == "160.0 ± 0.5"
    assert coefficient(record, "not_a_parameter") == "—"
    record["eos"]["fixed_parameters"] = ["V0"]
    assert coefficient(record, "V0") == "10.0 (fixed)"
    assert "10.0 | Not reported | Yes" in "\n".join(parameters(record))


def test_summary_includes_document_only_and_reference_variant_records():
    document = load_eosmat(ROOT / "peritheos/data/materials/wadsleyite.eosmat")
    assert "Document only" in render_material(document, executable=set(), defaults={})
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    parent = document["eos_records"][0]
    variant = copy.deepcopy(parent)
    variant["identifier"] = "example_scale_variant"
    variant["pressure_scale_variant_of"] = parent["identifier"]
    document["eos_records"].append(variant)
    sheet = render_material(document, executable={variant["identifier"]}, defaults={})
    assert "Pressure-scale reference variant" in sheet
    assert f"(#record-{variant['pressure_scale_variant_of']})" in sheet


def test_source_checked_record_can_still_be_not_recommended():
    document = convert_eosmat(load_eosmat(TEMPLATE), version=3)
    record = document["eos_records"][0]
    record["scientific_validation"] = {
        "status": "primary_source_validated",
        "usage_recommendation": "not_recommended_for_quantitative_use",
        "reproduction_status": "not_achieved",
    }
    executable = {record["identifier"]}
    assert (
        availability(record, executable)
        == "Executable; not recommended for quantitative use"
    )
    sheet = render_material(document, executable=executable, defaults={})
    assert sheet.index("**Usage qualification:**") < sheet.index("<details")
    assert "not\\_achieved" in sheet


def test_generated_material_summaries_and_formatting_are_current():
    for module in ("scripts.format_eosmat", "scripts.generate_material_summaries"):
        subprocess.run([sys.executable, "-m", module, "--check"], cwd=ROOT, check=True)
