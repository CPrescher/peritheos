"""Equation-level regression cases traceable to the cited literature."""

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos.eos.rt import BM3, ModifiedTait, Murnaghan, NaturalStrain3
from peritheos.eos.thermal import MieGruneisenDebye, ThermalModifiedTait

REFERENCE_FILE = Path(__file__).parent / "data" / "literature_reference_cases.json"
CASES = json.loads(REFERENCE_FILE.read_text())


@pytest.mark.parametrize("case", CASES, ids=lambda case: case["id"])
def test_literature_equation_reference_case(case):
    if case["model"] == "Murnaghan":
        model = Murnaghan(**case["parameters"])
    elif case["model"] == "NaturalStrain3":
        model = NaturalStrain3(**case["parameters"])
    elif case["model"] == "ModifiedTait":
        model = ModifiedTait(**case["parameters"])
    elif case["model"] == "ThermalModifiedTait":
        reference = ModifiedTait(**case["reference_parameters"])
        model = ThermalModifiedTait(reference, **case["parameters"])
    elif case["model"] == "MieGruneisenDebye":
        reference = BM3(**case["reference_parameters"])
        model = MieGruneisenDebye(reference, **case["parameters"])
    else:
        raise AssertionError(f"Unknown reference-case model: {case['model']}")

    arguments = [case["volume"]]
    if "temperature" in case:
        arguments.append(case["temperature"])
    value = getattr(model, case["quantity"])(*arguments)

    assert case["doi"]
    # SciPy's CODATA gas constant changed slightly across supported releases.
    assert np.isclose(value, case["expected"], rtol=1.0e-10, atol=1.0e-10)
