"""Equation-level regression cases traceable to the cited literature."""

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos.eos.rt import (
    BM2,
    BM3,
    BM4,
    Holzapfel,
    ModifiedTait,
    Murnaghan,
    NaturalStrain3,
    Vinet,
)
from peritheos.eos.thermal import (
    MieGruneisenDebye,
    MieGruneisenEinstein,
    ThermalModifiedTait,
    ThermalReferenceStateEOS,
)

REFERENCE_FILE = (
    Path(__file__).parent.parent
    / "crates"
    / "peritheos"
    / "tests"
    / "data"
    / "literature_reference_cases.json"
)
CASES = json.loads(REFERENCE_FILE.read_text())


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["id"])
def test_literature_equation_reference_case(case):
    room_temperature_models = {
        "BM2": BM2,
        "BM3": BM3,
        "BM4": BM4,
        "Holzapfel": Holzapfel,
        "ModifiedTait": ModifiedTait,
        "Murnaghan": Murnaghan,
        "NaturalStrain3": NaturalStrain3,
        "Vinet": Vinet,
    }
    if case["model"] in room_temperature_models:
        model = room_temperature_models[case["model"]](**case["parameters"])
    elif case["model"] == "ThermalModifiedTait":
        reference = ModifiedTait(**case["reference_parameters"])
        model = ThermalModifiedTait(reference, **case["parameters"])
    elif case["model"] in {"MieGruneisenDebye", "MieGruneisenEinstein"}:
        reference = BM3(**case["reference_parameters"])
        thermal_models = {
            "MieGruneisenDebye": MieGruneisenDebye,
            "MieGruneisenEinstein": MieGruneisenEinstein,
        }
        model = thermal_models[case["model"]](reference, **case["parameters"])
    elif case["model"] == "ThermalReferenceStateEOS":
        reference = BM3(**case["reference_parameters"])
        model = ThermalReferenceStateEOS(
            reference,
            **case["parameters"],
            **case["configuration"],
        )
    else:
        raise AssertionError(f"Unknown reference-case model: {case['model']}")

    arguments = [case["volume"]]
    if "temperature" in case:
        arguments.append(case["temperature"])
    value = getattr(model, case["quantity"])(*arguments)

    assert case["doi"]
    # SciPy's CODATA gas constant changed slightly across supported releases.
    assert np.isclose(value, case["expected"], rtol=1.0e-10, atol=1.0e-10)
