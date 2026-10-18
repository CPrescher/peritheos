"""Standalone authoring stays small, explicit, and compatible with full EOSMAT."""

import copy
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest
from jsonschema import Draft202012Validator

from peritheos import (
    Material,
    convert_eosmat,
    create_eosmat,
    eosmat_schema,
    load_eosmat,
    save_eosmat,
)
from peritheos.__main__ import main
from peritheos.eos.rt import BM3
from peritheos.errors import MaterialError

ROOT = Path(__file__).resolve().parents[1]
COEFFICIENTS = {"V0": 74.7, "K0": 160.0, "K0_prime": 4.2}


def create(**options):
    arguments = {"temperature_ref": 300, **options}
    return create_eosmat("my_material", "MgO", "BM3", COEFFICIENTS, **arguments)


@pytest.mark.parametrize("model", ["BM3", "birch_murnaghan_3"])
def test_standalone_file_is_small_valid_and_numerically_usable(tmp_path, model):
    document = create_eosmat(
        "my_material", "MgO", model, COEFFICIENTS, temperature_ref=300
    )
    Draft202012Validator(eosmat_schema(4)).validate(document)
    path = tmp_path / "my_material.eosmat"
    save_eosmat(path, document)
    assert len(path.read_text().splitlines()) <= 25
    assert load_eosmat(path) == document
    assert "datasets" not in document and "source" not in document
    assert "formula_units_per_cell" not in document
    record = document["eos_records"][0]
    assert record["model"] == "BM3"
    assert record["parameters"] == COEFFICIENTS
    assert "eos" not in record and "parameter_errors" not in record
    assert "evidence" not in record
    with pytest.raises(MaterialError, match="pending_primary_source_check"):
        Material.from_eosmat(document)
    material = Material.from_eosmat(document, require_primary_validation=False)
    volumes = np.array([74.7, 70, 65])
    np.testing.assert_allclose(
        material.eos_records[0].pressure(volumes), BM3(**COEFFICIENTS).pressure(volumes)
    )
    for pressure in (0, 10, 50):
        record = material.eos_records[0]
        assert record.pressure(record.volume(pressure)) == pytest.approx(
            pressure, abs=1e-8
        )


def test_authoring_preserves_supplied_errors_precision_and_source_without_mutation():
    reference = {"authors": ["Example"], "year": 2026, "location": "Table 2"}
    errors = {"V0": None, "K0": 0.0, "K0_prime": 0.013456789012345678}
    fixed = ["K0_prime"]
    before = copy.deepcopy((reference, errors, fixed, COEFFICIENTS))
    document = create(
        reference=reference, parameter_errors=errors, fixed_parameters=fixed
    )
    record = document["eos_records"][0]
    legacy = convert_eosmat(document, version=3)["eos_records"][0]
    assert legacy["parameter_errors"] == errors
    assert legacy["fixed_parameters"] == fixed
    assert record["reference"] == reference
    record["reference"]["authors"].append("Other")
    record["parameters"]["K0"]["error"] = 99
    record["parameters"]["K0_prime"]["fixed"] = False
    record["parameters"]["K0"]["value"] = 99
    assert (reference, errors, fixed, COEFFICIENTS) == before


@pytest.mark.parametrize("model", ["MieGruneisenDebye", "mie_gruneisen_debye"])
@pytest.mark.parametrize("law", ["integrated_gruneisen", "variable_exponent"])
def test_thermal_component_keeps_cell_basis_options_errors_and_missing_uncertainty(
    model,
    law,
):
    thermal = {
        "model": model,
        "parameters": {"Tr": 300, "theta0": 760, "gamma0": 1.5, "q": 1, "n": 2},
        "debye_temperature_law": law,
        "parameter_errors": {"theta0": None, "q": 0.1},
        "fixed_parameters": ["n"],
        "extension": {"original_tokens": ["760", None]},
    }
    before = copy.deepcopy(thermal)
    with pytest.raises(MaterialError, match="formula_units_per_cell"):
        create(thermal=thermal)
    document = create(thermal=thermal, formula_units_per_cell=4)
    Draft202012Validator(eosmat_schema(4)).validate(document)
    component = convert_eosmat(document, version=3)["eos_records"][0]["thermal"]
    assert component["model"] == "mie_gruneisen_debye"
    assert component["type"] == "MieGruneisenDebye"
    assert {k: v for k, v in component.items() if k not in {"type", "model"}} == {
        k: v for k, v in thermal.items() if k != "model"
    }
    assert thermal == before
    record = Material.from_eosmat(
        document, require_primary_validation=False
    ).eos_records[0]
    assert record.eos.rt_eos.V0 == pytest.approx(74.7 * 6.02214076e-2 / 4)
    assert record.eos.debye_temperature_law == law
    assert record.pressure(74.7, 1000) > record.pressure(74.7, 300)


def test_zero_kelvin_isotherm_remains_explicit():
    document = create(temperature_ref=0)
    assert document["eos_records"][0]["temperature_ref"] == 0
    record = Material.from_eosmat(
        document, require_primary_validation=False
    ).eos_records[0]
    assert record.pressure(74.7) == 0


@pytest.mark.parametrize(
    "options,match",
    [
        ({"fixed_parameters": "K0"}, "sequence"),
        ({"fixed_parameters": ["missing"]}, "unique supplied"),
        ({"fixed_parameters": ["K0", "K0"]}, "unique supplied"),
        ({"parameter_errors": {"missing": 1}}, "supplied EOS parameters"),
        ({"parameter_errors": {"K0": float("nan")}}, "finite"),
        ({"parameter_errors": {"K0": -1}}, "nonnegative"),
        ({"temperature_ref": float("nan")}, "finite"),
        ({"temperature_ref": True}, "number"),
        ({"name": ""}, "nonempty"),
        ({"reference": ""}, "nonempty"),
        ({"thermal": {}}, "type or model"),
        ({"thermal": {"model": "unknown"}}, "Unknown EOS model"),
        (
            {"thermal": {"model": "multi_oscillator_gruneisen_thermal_pressure"}},
            "explicit formulation",
        ),
        (
            {
                "thermal": {
                    "type": "MieGruneisenDebye",
                    "model": "mie_gruneisen_einstein",
                }
            },
            "does not match",
        ),
    ],
)
def test_invalid_or_ambiguous_input_is_not_silently_repaired(options, match):
    with pytest.raises(ValueError, match=match):
        create(**options)


@pytest.mark.parametrize(
    "identifier,formula,model,parameters",
    [
        ("My material", "MgO", "BM3", COEFFICIENTS),
        ("my_material", "", "BM3", COEFFICIENTS),
        ("my_material", "MgO", "unknown", COEFFICIENTS),
        ("my_material", "MgO", "BM3", {"V0": 74.7, "K0": 160}),
        ("my_material", "MgO", "BM3", {**COEFFICIENTS, "K0": -1}),
        ("my_material", "MgO", "BM3", {**COEFFICIENTS, "K0": True}),
        ("my_material", "MgO", "LinearUsUpHugoniot", {}),
    ],
)
def test_bad_identity_coefficients_and_missing_hugoniot_state_fail(
    identifier, formula, model, parameters
):
    with pytest.raises(ValueError):
        create_eosmat(identifier, formula, model, parameters, temperature_ref=300)


def command(path):
    return [
        "create",
        str(path),
        "--formula",
        "MgO",
        "--model",
        "BM3",
        "--temperature",
        "300",
        "--parameters",
        "V0=74.7",
        "K0=160",
        "K0_prime=4.2",
    ]


def test_module_command_creates_an_independent_file(tmp_path):
    path = tmp_path / "my_material.eosmat"
    result = subprocess.run(
        [sys.executable, "-m", "peritheos", *command(path)],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(ROOT)},
        text=True,
        capture_output=True,
    )
    assert result.returncode == 0, result.stderr
    assert "EOSMAT v4" in result.stdout
    assert load_eosmat(path) == create()


def test_command_refuses_to_overwrite_existing_files(tmp_path, capsys):
    path = tmp_path / "my_material.eosmat"
    path.write_text("original file")
    with pytest.raises(SystemExit) as error:
        main(command(path))
    assert error.value.code == 2
    assert path.read_text() == "original file"
    assert "File exists" in capsys.readouterr().err


def test_command_names_missing_coefficients_without_a_traceback(tmp_path, capsys):
    path = tmp_path / "my_material.eosmat"
    with pytest.raises(SystemExit) as error:
        main(command(path) + ["--parameters", "V0=74.7", "K0=160"])
    assert error.value.code == 2
    message = capsys.readouterr().err
    assert "K0_prime" in message
    assert "Traceback" not in message
    assert not path.exists()


def test_command_supports_thermal_parameters_and_explicit_cell_contents(tmp_path):
    path = tmp_path / "my_material.eosmat"
    assert (
        main(
            command(path)
            + [
                "--z",
                "4",
                "--errors",
                "K0=2",
                "--fixed",
                "K0_prime",
                "--thermal-model",
                "MieGruneisenDebye",
                "--thermal-parameters",
                "Tr=300",
                "theta0=760",
                "gamma0=1.5",
                "q=1",
                "n=2",
            ]
        )
        == 0
    )
    document = load_eosmat(path)
    assert document["formula_units_per_cell"] == 4
    record = document["eos_records"][0]
    assert record["thermal"]["model"] == "MieGruneisenDebye"
    assert record["thermal"]["parameters"]["n"] == 2
    assert record["parameters"]["K0"]["error"] == 2
    assert record["parameters"]["K0_prime"]["fixed"] is True


@pytest.mark.parametrize(
    "extra",
    [
        ["--parameters", "V0=nan", "K0=160", "K0_prime=4"],
        ["--parameters", "V0=74", "V0=75"],
        ["--parameters", "V0:74"],
        ["--parameters", "V0=word"],
        ["--thermal-model", "MieGruneisenDebye"],
        ["--thermal-parameters", "Tr=300"],
    ],
)
def test_command_reports_errors_without_creating_files(tmp_path, extra):
    path = tmp_path / "my_material.eosmat"
    with pytest.raises(SystemExit) as error:
        main(command(path) + extra)
    assert error.value.code == 2
    assert not path.exists()


def test_full_v3_metadata_can_be_added_without_a_second_format(tmp_path):
    document = create()
    record = document["eos_records"][0]
    document["extension"] = {"original": [None, -0.0, 1.2345678901234567]}
    document["datasets"] = [
        {
            "identifier": "my_observations",
            "kind": "pressure_volume",
            "reference": "Private experiment",
            "source_location": "Run 1",
            "used_by_eos_records": [record["identifier"]],
            "columns": [
                {"name": "p", "quantity": "pressure", "unit": "GPa", "role": "value"},
                {
                    "name": "v",
                    "quantity": "unit_cell_volume",
                    "unit": "angstrom^3",
                    "role": "value",
                },
            ],
            "rows": [[0, 74.7], [10, None]],
        }
    ]
    record["evidence"] = {"parameter_provenance": {"K0": "Private unweighted fit"}}
    record["parameter_covariance"] = {
        "parameter_order": ["K0", "K0_prime"],
        "matrix": [[1, 0.01], [0.01, 0.01]],
    }
    before = copy.deepcopy(document)
    path = tmp_path / "full.eosmat"
    save_eosmat(path, document)
    assert load_eosmat(path) == document == before
    assert json.loads(path.read_text())["datasets"][0]["rows"][1][1] is None
