"""Cross-language tests for the public Peritheos error contract."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from peritheos import (
    EosError,
    EosmatError,
    EosValidationError,
    FitEosNumericalError,
    FitEosValidationError,
    FitError,
    FitValidationError,
    MaterialLookupError,
    PeritheosError,
    _rust,
    get_material_document,
    load_eosmat,
    validate_eosmat_document,
)


def test_error_hierarchy_preserves_builtin_compatibility() -> None:
    eos_error = EosValidationError("bad state")
    fit_error = FitValidationError("bad fit")
    lookup_error = MaterialLookupError("missing")

    assert isinstance(eos_error, (PeritheosError, EosError, ValueError))
    assert isinstance(fit_error, (PeritheosError, FitError, ValueError, RuntimeError))
    assert isinstance(lookup_error, (PeritheosError, KeyError))
    assert isinstance(
        FitEosValidationError("bad EOS state"),
        (FitError, EosError, ValueError, RuntimeError),
    )
    assert isinstance(
        FitEosNumericalError("EOS did not converge"),
        (FitError, EosError, ArithmeticError, RuntimeError),
    )


def test_error_metadata_is_structured_and_context_is_read_only() -> None:
    error = EosValidationError(
        "volume is invalid",
        code="eos.invalid_state",
        operation="pressure",
        field="volume",
        context={"value": -1.0},
    )

    assert error.args == ("volume is invalid",)
    assert error.code == "eos.invalid_state"
    assert error.operation == "pressure"
    assert error.field == "volume"
    assert error.context == {"value": -1.0}
    with pytest.raises(TypeError):
        error.context["value"] = 1.0  # type: ignore[index]


def test_pure_python_domains_raise_specific_errors() -> None:
    with pytest.raises(EosmatError) as captured:
        validate_eosmat_document({})
    assert captured.value.code == "eosmat.invalid_document"

    with pytest.raises(MaterialLookupError) as captured:
        get_material_document("not-a-material")
    assert captured.value.code == "material.not_found"


def test_invalid_eosmat_json_retains_decoder_as_cause(tmp_path: Path) -> None:
    path = tmp_path / "invalid.eosmat"
    path.write_text('{"format_version":', encoding="utf-8")

    with pytest.raises(EosmatError) as captured:
        load_eosmat(path)

    assert captured.value.code == "eosmat.json"
    assert captured.value.operation == "load"
    assert captured.value.context["path"] == str(path)
    assert isinstance(captured.value.__cause__, json.JSONDecodeError)


def test_native_eos_errors_use_the_public_python_hierarchy() -> None:
    with pytest.raises(EosValidationError) as captured:
        _rust.RtEos.bm2(0.0, 100.0)

    error = captured.value
    assert error.code == "eos.invalid_parameter"
    assert error.operation == "eos"
    assert error.field == "V0"
    assert isinstance(error, ValueError)


def test_native_fit_validation_uses_the_public_python_hierarchy() -> None:
    with pytest.raises(FitValidationError) as captured:
        _rust.fit_least_squares(
            lambda parameters: np.array([parameters[0]]),
            np.array([0.0]),
            np.array([-1.0]),
            np.array([1.0]),
            global_parameter_count=1,
        )

    assert captured.value.code == "fit.invalid_input"
    assert isinstance(captured.value, (FitError, ValueError))
