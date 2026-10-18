"""Create native EOSMAT files without database audit boilerplate."""

from __future__ import annotations

import copy
import re
from collections.abc import Mapping, Sequence
from typing import Any

from peritheos.eosmat import (
    _RT_MODELS,
    _THERMAL_MODELS,
    EOSMAT_FORMAT,
    EOSMAT_FORMAT_VERSION,
    _finite_number,
    convert_eosmat,
    validate_eosmat_document,
)
from peritheos.errors import EosmatError
from peritheos.materials import Material


def _component(model: str, models: Mapping[str, str]) -> dict[str, Any]:
    if not isinstance(model, str):
        raise EosmatError("model must be an exact supported model name")
    if model in models:
        return {"type": model, "model": models[model]}
    matches = [key for key, value in models.items() if value == model]
    if len(matches) == 1:
        return {"type": matches[0], "model": model}
    if matches:
        raise EosmatError(
            f"Model {model!r} needs an explicit formulation: {', '.join(matches)}"
        )
    raise EosmatError(f"Unknown EOS model {model!r}; use an exact supported model name")


def create_eosmat(
    identifier: str,
    formula: str,
    model: str,
    parameters: Mapping[str, float],
    *,
    temperature_ref: float,
    name: str | None = None,
    reference: str | Mapping[str, Any] | None = None,
    parameter_errors: Mapping[str, float | None] | None = None,
    fixed_parameters: Sequence[str] = (),
    formula_units_per_cell: float | None = None,
    thermal: Mapping[str, Any] | None = None,
    format_version: int = 4,
) -> dict[str, Any]:
    """Build a native v4, executable but unaudited standalone material document.

    ``model`` accepts an exact EOS type (e.g. ``BM3``) or its canonical model
    identifier. Parameters use GPa, kelvin and conventional-cell Å³. The
    reference temperature is mandatory; no measurement conditions, cell
    contents, uncertainties or validation findings are inferred.

    ``thermal`` accepts the regular component mapping, with either ``type``
    or ``model`` supplied; model-specific options and metadata are preserved.
    Molar-volume thermal models require explicit ``formula_units_per_cell``.
    Hugoniot initial states and advanced record metadata use the full EOSMAT schema.

    The returned dictionary is independent of its inputs. V4 keeps values,
    errors and fixed flags together, and detailed source checks in evidence.
    ``format_version=3`` explicitly exports to historical consumers.
    Loading it for calculations requires explicit unaudited acceptance through
    ``Material.from_eosmat(document, require_primary_validation=False)``.
    """
    if not isinstance(identifier, str) or not re.fullmatch(
        r"[a-z0-9]+(?:_[a-z0-9]+)*", identifier
    ):
        raise EosmatError(
            "identifier must use lower-case words separated by underscores"
        )
    if not isinstance(formula, str) or not formula.strip():
        raise EosmatError("formula must be a nonempty string")
    if name is not None and (not isinstance(name, str) or not name.strip()):
        raise EosmatError("name must be a nonempty string")
    if _finite_number(temperature_ref, "temperature_ref") < 0:
        raise EosmatError("temperature_ref must be nonnegative")
    eos = _component(model, _RT_MODELS)
    if eos["type"] == "LinearUsUpHugoniot":
        raise EosmatError(
            "Hugoniots require the full EOSMAT initial-state and branch metadata"
        )
    if not isinstance(parameters, Mapping):
        raise EosmatError("parameters must be a mapping")
    eos["parameters"] = copy.deepcopy(dict(parameters))
    if parameter_errors is not None and not isinstance(parameter_errors, Mapping):
        raise EosmatError("parameter_errors must be a mapping")
    errors = copy.deepcopy(dict(parameter_errors or {}))
    if set(errors) - set(parameters):
        raise EosmatError("parameter_errors must name supplied EOS parameters")
    for key, value in errors.items():
        if value is not None and _finite_number(value, f"parameter_errors.{key}") < 0:
            raise EosmatError("parameter_errors must be nonnegative when reported")
    if (
        not isinstance(fixed_parameters, Sequence)
        or isinstance(fixed_parameters, str)
        or not all(isinstance(key, str) for key in fixed_parameters)
    ):
        raise EosmatError("fixed_parameters must be a sequence of parameter names")
    fixed = list(fixed_parameters)
    if len(set(fixed)) != len(fixed) or set(fixed) - set(parameters):
        raise EosmatError(
            "fixed_parameters must contain unique supplied parameter names"
        )
    if reference is not None and (
        not isinstance(reference, (str, Mapping)) or not reference
    ):
        raise EosmatError("reference must be a nonempty citation string or mapping")
    record: dict[str, Any] = {
        "identifier": f"{identifier}_{eos['type'].lower()}",
        "label": f"{name or identifier}: {eos['type']}",
        "record_kind": "diagnostic",
        "reference": copy.deepcopy(reference)
        if reference is not None
        else "User-supplied parameters",
        "temperature_ref": temperature_ref,
        "eos": eos,
        "parameter_errors": errors,
        "fixed_parameters": fixed,
        "scientific_validation": {"status": "pending_primary_source_check"},
    }
    if thermal is not None:
        if not isinstance(thermal, Mapping):
            raise EosmatError("thermal must be a component mapping")
        component = copy.deepcopy(dict(thermal))
        selected = component.get("type", component.get("model"))
        if not isinstance(selected, str):
            raise EosmatError("thermal must specify a supported type or model")
        resolved = _component(selected, _THERMAL_MODELS)
        if (
            "type" in component
            and "model" in component
            and component["model"] != resolved["model"]
        ):
            raise EosmatError("thermal.model does not match thermal.type")
        component.update(resolved)
        record["thermal"] = component
    document: dict[str, Any] = {
        "format": EOSMAT_FORMAT,
        "format_version": 3,
        "identifier": identifier,
        "name": name or identifier,
        "formula": formula,
        "units": {
            "pressure": "GPa",
            "temperature": "K",
            "volume": "angstrom^3/conventional_unit_cell",
        },
        "eos_records": [record],
    }
    if formula_units_per_cell is not None:
        document["formula_units_per_cell"] = formula_units_per_cell
    validate_eosmat_document(document)
    # Check that the chosen model can actually be constructed. Keep the authored
    # declaration, rather than replacing it with a model's expanded defaults.
    Material.from_eosmat(document, require_primary_validation=False)
    native = convert_eosmat(document, version=EOSMAT_FORMAT_VERSION)
    record = native["eos_records"][0]
    record.pop("label", None)
    if reference is None:
        record.pop("reference", None)
    record.pop("evidence", None)
    return convert_eosmat(native, version=format_version)
