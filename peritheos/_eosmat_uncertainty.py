"""Uncertainty unit conversion at the EOSMAT version boundary."""

from __future__ import annotations

import math
from collections.abc import Mapping
from typing import Any

from peritheos.errors import EosmatError


def volume_scale(document: Mapping[str, Any], record: Mapping[str, Any]) -> float:
    """Return the declared conversion from coefficient V0 to runtime volume."""
    from peritheos.materials import _MOLAR_VOLUME_THERMAL_MODELS

    volume = record.get("volume")
    explicit = (
        volume.get("public_to_model_scale") if isinstance(volume, Mapping) else None
    )
    if explicit is not None:
        if isinstance(explicit, bool) or not isinstance(explicit, (int, float)):
            raise EosmatError("Uncertainty volume conversion scale must be numeric")
        scale = float(explicit)
    elif (
        record.get("thermal", {}).get("model") in _MOLAR_VOLUME_THERMAL_MODELS
        or record["eos"]["model"] == "maltby_2024_published"
    ):
        z = document.get("formula_units_per_cell")
        if z is None:
            raise EosmatError(
                "Uncertainty conversion requires formula_units_per_cell or "
                "volume.public_to_model_scale for a molar-volume model"
            )
        if (
            isinstance(z, bool)
            or not isinstance(z, (int, float))
            or not math.isfinite(z)
            or z <= 0
        ):
            raise EosmatError("formula_units_per_cell must be positive and finite")
        scale = 0.0602214076 / float(z)
    else:
        scale = 1.0
    if not math.isfinite(scale) or scale <= 0:
        raise EosmatError(
            "Uncertainty volume conversion scale must be positive and finite"
        )
    return scale


def _scaled(value: float, factor: float) -> float:
    if factor == 1:
        return value
    result = value * factor
    if not math.isfinite(result) or value != 0 and result == 0:
        raise EosmatError("Uncertainty unit conversion overflowed or underflowed")
    return result


def mark_parameter_units(document: dict[str, Any]) -> None:
    """Mark a v5 document's expanded v3 construction view explicitly."""
    for record in document["eos_records"]:
        errors = record.get("parameter_errors", {})
        covariance = record.get("parameter_covariance")
        order = (
            covariance.get("parameter_order") if isinstance(covariance, dict) else None
        )
        has_volume_covariance = isinstance(order, list) and any(
            name in ("V0", "rt_eos.V0") for name in order
        )
        if errors.get("V0") is None and not has_volume_covariance:
            continue
        if volume_scale(document, record) == 1:
            continue
        if errors.get("V0") is not None:
            record["parameter_error_volume_basis"] = "public"
        if has_volume_covariance:
            covariance["volume_basis"] = "public"


def convert_units(document: dict[str, Any], *, parameter_units: bool) -> None:
    """Convert a v3-shaped document in place, then remove historical flags."""
    target = "public" if parameter_units else "model"
    for record in document["eos_records"]:
        errors = record.get("parameter_errors", {})
        error_basis = record.pop("parameter_error_volume_basis", "model")
        covariance = record.get("parameter_covariance")
        covariance_basis = (
            covariance.pop("volume_basis", "model")
            if isinstance(covariance, dict)
            else "model"
        )
        convert_error = errors.get("V0") is not None and error_basis != target
        convert_covariance = (
            isinstance(covariance, dict)
            and covariance_basis != target
            and any(
                name in {"V0", "rt_eos.V0"} for name in covariance["parameter_order"]
            )
        )
        if not convert_error and not convert_covariance:
            continue
        scale = volume_scale(document, record)
        factor = 1 / scale if parameter_units else scale
        if convert_error:
            errors["V0"] = _scaled(errors["V0"], factor)
        if convert_covariance:
            order = covariance["parameter_order"]
            factors = [factor if name in {"V0", "rt_eos.V0"} else 1 for name in order]
            covariance["matrix"] = [
                [
                    _scaled(_scaled(value, factors[i]), factors[j])
                    for j, value in enumerate(row)
                ]
                for i, row in enumerate(covariance["matrix"])
            ]
