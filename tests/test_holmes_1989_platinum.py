import importlib.util
import json
from pathlib import Path

import pytest
from jsonschema import Draft202012Validator

from peritheos import get_material_document

ROOT = Path(__file__).parents[1]


def _load_reproduction():
    path = ROOT / "scripts" / "reproduce_holmes_1989_platinum.py"
    spec = importlib.util.spec_from_file_location(
        "reproduce_holmes_1989_platinum", path
    )
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_holmes_equation_11_is_exactly_executable_but_not_refitted_to_shock_rows():
    result = _load_reproduction().reproduce()

    curve = result["equation_11_checkpoints"]
    assert curve["rows"] == 8
    assert curve["pressure_range_gpa"] == pytest.approx([0.0, 550.0], abs=2e-10)
    assert curve["csv_vs_equation_max_abs_gpa"] < 5e-12
    assert curve["peritheos_vs_equation_max_abs_gpa"] < 5e-10

    shock = result["table_iii_shock_qualification"]
    assert shock["rows"] == 7
    assert shock["pressure_range_gpa"] == [218.9, 659.3]
    assert shock["rankine_hugoniot_momentum_max_abs_rounding_gpa"] < 0.15
    assert shock["rankine_hugoniot_mass_max_abs_rounding_g_cm3"] < 0.008
    assert shock["fit_role"] == "qualification_only_not_300k_isotherm_observations"

    assert result["reconstructability"] == {
        "published_300k_curve": "exact_analytical_reconstruction",
        "underlying_lmto_fit": "not_refittable_no_numerical_e0_or_p0_grid",
        "theoretical_hugoniot": "not_reconstructable_no_numerical_thermal_grid",
    }


def test_holmes_record_preserves_raw_source_coefficients_and_dataset_roles():
    document = get_material_document("platinum")
    record = next(
        row
        for row in document["eos_records"]
        if row["identifier"] == "platinum_holmes_1989_vinet_1"
    )
    assert record["source_parameterization"] == {
        "equation": "Equation (11)",
        "V0_bohr3_per_atom": 101.9,
        "P_T_gpa": 798.31,
        "eta": 7.2119,
        "equivalent_parameter_mapping": {
            "K0": "P_T/3",
            "K0_prime": "1+eta/1.5",
        },
        "table_iv_rounded_parameters": {"K0_gpa": 266.0, "K0_prime": 5.81},
    }
    assert record["eos"]["parameters"]["K0"] == pytest.approx(798.31 / 3.0)
    assert record["eos"]["parameters"]["K0_prime"] == pytest.approx(1.0 + 7.2119 / 1.5)
    assert record["scientific_validation"]["primary_data_check"][
        "dataset_identifiers"
    ] == [
        "platinum_holmes_1989_table3_shock",
        "platinum_holmes_1989_equation11_checkpoints",
    ]
    datasets = {row["identifier"]: row for row in document["datasets"]}
    assert datasets["platinum_holmes_1989_table3_shock"]["kind"] == "shock_hugoniot"
    checkpoints = datasets["platinum_holmes_1989_equation11_checkpoints"]
    assert checkpoints["kind"] == "source_derived_eos_checkpoints"
    assert "not independent fit observations" in checkpoints["notes"]
    assert "fit_datasets" not in record


def test_holmes_platinum_document_matches_normative_schema():
    schema = json.loads(
        (ROOT / "peritheos/data/eosmat-v3.schema.json").read_text(encoding="utf-8")
    )
    document = get_material_document("platinum")
    assert list(Draft202012Validator(schema).iter_errors(document)) == []


def test_holmes_audit_documents_the_missing_lmto_inputs():
    text = (ROOT / "docs/literature-reproductions/holmes-1989-platinum.md").read_text(
        encoding="utf-8"
    )
    for phrase in (
        "E0(V)",
        "P0(V)",
        "density-of-states values",
        "not independent fit observations",
        "never fitted as a 300 K isotherm",
    ):
        assert phrase in text
