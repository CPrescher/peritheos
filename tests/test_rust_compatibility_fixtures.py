"""Lock language-independent compatibility fixtures to the Python 0.5 oracle."""

import json
from pathlib import Path

import numpy as np

from peritheos.eos.rt import (
    BM2,
    BM3,
    BM4,
    Holzapfel,
    ModifiedTait,
    Murnaghan,
    NaturalStrain2,
    NaturalStrain3,
    NaturalStrain4,
    Vinet,
)
from peritheos.eos.thermal import (
    MieGruneisenDebye,
    MieGruneisenEinstein,
    Sokolova2016,
    ThermalModifiedTait,
)
from peritheos.eos.thermal.mie_gruneisen import _debye_function_3

FIXTURE = Path(__file__).parent / "data" / "isothermal_compatibility_cases.json"
THERMAL_FIXTURE = Path(__file__).parent / "data" / "thermal_compatibility_cases.json"


def test_isothermal_compatibility_fixture_matches_python_oracle():
    document = json.loads(FIXTURE.read_text())
    assert document["schema_version"] == 1
    assert document["source_commit"] == "12d033378418cfd6c9ece6050c550fc748ffe02a"
    tolerance = document["relative_tolerance"]
    classes = {
        model.__name__: model
        for model in (
            BM2,
            BM3,
            BM4,
            Holzapfel,
            ModifiedTait,
            Murnaghan,
            NaturalStrain2,
            NaturalStrain3,
            NaturalStrain4,
            Vinet,
        )
    }

    for case in document["cases"]:
        model = classes[case["model"]](**case["parameters"])
        volumes = model.V0 * np.asarray(case["volume_fractions"])
        pressures = np.asarray(case["pressures"])
        bulk_moduli = np.asarray(case["bulk_moduli"])

        assert np.allclose(
            model.pressure(volumes), pressures, rtol=tolerance, atol=tolerance
        )
        assert np.allclose(
            model.bulk_modulus(volumes), bulk_moduli, rtol=tolerance, atol=tolerance
        )
        assert np.allclose(model.volume(pressures), volumes, rtol=1.0e-10)


def _thermal_fixture_models():
    reference = BM3(1.0, 160.0, 4.0)
    holzapfel = Holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0)
    reduced_parameters = (
        298.15,
        684.0,
        0.564,
        1561.0,
        2.436,
        -0.506,
        1.085,
    )
    return {
        "MieGruneisenDebye": MieGruneisenDebye(reference, 300.0, 800.0, 1.5, 1.0, 2.0),
        "MieGruneisenEinstein": MieGruneisenEinstein(
            reference, 300.0, 800.0, 1.5, 1.0, 2.0
        ),
        "ThermalModifiedTait": ThermalModifiedTait(
            ModifiedTait(1.0, 160.0, 4.0, -0.01),
            298.15,
            700.0,
            2.5e-5,
            2.0,
        ),
        "Sokolova2016Reduced": Sokolova2016(
            holzapfel, *reduced_parameters, 0.0, 0.0, 0.0, 0.0
        ),
        "Sokolova2016Complete": Sokolova2016(
            holzapfel,
            *reduced_parameters,
            5.2,
            1.3,
            0.8,
            2.7,
            beta=0.35,
            QBo=480.0,
            d=2.4,
            mb=0.75,
            QB1o=1120.0,
            d1=1.6,
            mb1=0.4,
        ),
    }


def test_thermal_compatibility_fixture_matches_python_oracle():
    document = json.loads(THERMAL_FIXTURE.read_text())
    assert document["schema_version"] == 1
    tolerance = document["relative_tolerance"]
    for case in document["debye_function_3"]:
        assert np.isclose(
            _debye_function_3(case["argument"]),
            case["value"],
            rtol=tolerance,
            atol=tolerance,
        )

    models = _thermal_fixture_models()
    for case in document["cases"]:
        model = models[case["model"]]
        for state in case["states"]:
            volume = state["volume"]
            temperature = state["temperature"]
            for quantity, expected in state["quantities"].items():
                evaluator = getattr(model, quantity)
                value = (
                    evaluator(volume)
                    if quantity == "characteristic_temperature"
                    else evaluator(volume, temperature)
                )
                assert np.isclose(value, expected, rtol=tolerance, atol=tolerance), (
                    case["model"],
                    quantity,
                )
