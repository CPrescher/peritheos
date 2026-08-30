"""Direct tests for the private native extension during the additive phase."""

import numpy as np
import pytest

from peritheos import _rust
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


@pytest.mark.parametrize(
    ("python_model", "native_model"),
    [
        (BM2(10.0, 120.0), _rust.RtEos.bm2(10.0, 120.0)),
        (BM3(10.0, 120.0, 4.3), _rust.RtEos.bm3(10.0, 120.0, 4.3)),
        (
            BM4(10.0, 120.0, 4.3, -0.02),
            _rust.RtEos.bm4(10.0, 120.0, 4.3, -0.02),
        ),
        (
            Murnaghan(10.0, 120.0, 4.3),
            _rust.RtEos.murnaghan(10.0, 120.0, 4.3),
        ),
        (
            ModifiedTait(10.0, 120.0, 4.3, -0.02),
            _rust.RtEos.modified_tait(10.0, 120.0, 4.3, -0.02),
        ),
        (
            NaturalStrain2(10.0, 120.0),
            _rust.RtEos.natural_strain2(10.0, 120.0),
        ),
        (
            NaturalStrain3(10.0, 120.0, 4.3),
            _rust.RtEos.natural_strain3(10.0, 120.0, 4.3),
        ),
        (
            NaturalStrain4(10.0, 120.0, 4.3, -0.02),
            _rust.RtEos.natural_strain4(10.0, 120.0, 4.3, -0.02),
        ),
        (Vinet(10.0, 120.0, 4.3), _rust.RtEos.vinet(10.0, 120.0, 4.3)),
        (
            Holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
            _rust.RtEos.holzapfel(0.3414, 441.5, 3.9, 1.0, 6.0),
        ),
    ],
)
def test_native_isothermal_binding_matches_python(python_model, native_model):
    fractions = np.array([[0.55, 0.7], [0.85, 1.0]], dtype=float)
    volumes = python_model.V0 * fractions
    expected_pressures = python_model.pressure(volumes)

    assert native_model.reference_volume == python_model.V0
    assert np.allclose(
        native_model.pressure_array(volumes), expected_pressures, rtol=1.0e-12
    )
    assert np.allclose(
        native_model.bulk_modulus_array(volumes),
        python_model.bulk_modulus(volumes),
        rtol=1.0e-12,
    )
    assert np.allclose(
        native_model.volume_array(expected_pressures), volumes, rtol=1.0e-10
    )
    assert isinstance(native_model.pressure_scalar(float(volumes[0, 0])), float)
    assert native_model.pressure_array(volumes).shape == volumes.shape


def test_native_binding_preserves_error_categories():
    with pytest.raises(ValueError):
        _rust.RtEos.bm2(0.0, 100.0)

    model = _rust.RtEos.bm2(10.0, 100.0)
    with pytest.raises(ValueError):
        model.pressure_scalar(0.0)
    with pytest.raises(ValueError):
        model.volume_scalar(-100.0)
