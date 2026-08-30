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

FIXTURE = Path(__file__).parent / "data" / "isothermal_compatibility_cases.json"


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
