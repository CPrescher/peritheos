"""Explicit, source-preserving observation views for EOS consumers.

Tables are archives; views describe the sample rows and coordinates to display.
Nothing here fits an EOS, fills a missing temperature, or changes source values.
The packaged registry is a migration layer for existing EOSMAT v3 documents.
New documents can declare the same ``observation_views`` directly on a dataset.
"""

from __future__ import annotations

import copy
import hashlib
import json
from dataclasses import dataclass
from importlib import resources
from typing import Any, Mapping, NoReturn

import numpy as np

from peritheos.datasets import Dataset, _convert_values
from peritheos.errors import DatasetError, EosmatError

EVIDENCE_KINDS = (
    "experimental",
    "digitized",
    "shock",
    "experimental_reduction",
    "simulation",
    "model_output",
    "source_comparison",
)
CELL_UNIT = "angstrom^3/conventional_unit_cell"
_LENGTHS = ("a", "b", "c")
_ANGLES = ("alpha", "beta", "gamma")
_METRIC = _LENGTHS + _ANGLES


def observation_source_fingerprint(dataset: Mapping[str, Any]) -> str:
    """Bind a view to the whole source declaration, including resource checksum."""
    source = {k: v for k, v in dataset.items() if k != "observation_views"}
    return hashlib.sha256(
        json.dumps(
            source,
            sort_keys=True,
            separators=(",", ":"),
            ensure_ascii=False,
            allow_nan=False,
        ).encode()
    ).hexdigest()


def observation_material_fingerprint(document: Mapping[str, Any]) -> str:
    """Bind cell conversions to material identity and the declared cell basis."""
    fields = (
        "identifier",
        "formula",
        "formula_units_per_cell",
        "units",
        "symmetry",
        "phase",
        "lattice",
        "cell_setting",
        "cell_contents",
    )
    return observation_source_fingerprint(
        {k: document[k] for k in fields if k in document}
    )


def observation_registry() -> dict[str, Any]:
    """Return an independent copy of the packaged, source-bound migration views."""
    return json.loads(
        resources.files("peritheos.data")
        .joinpath("observation-views.json")
        .read_text(encoding="utf-8")
    )


def validate_observation_view(
    view: Mapping[str, Any],
    dataset: Mapping[str, Any],
    record_ids: set[str],
) -> None:
    """Validate references as well as structure; no data file is opened here."""

    def fail(reason: str) -> NoReturn:
        raise EosmatError(f"Observation view {view.get('identifier')!r}: {reason}")

    if not isinstance(view, Mapping):
        raise EosmatError("Observation view must be an object")

    def object_keys(value: Any, allowed: set[str], label: str) -> Mapping[str, Any]:
        if not isinstance(value, Mapping):
            fail(f"{label} must be an object")
        if set(value) - allowed:
            fail(f"{label} has unknown properties")
        return value

    object_keys(
        view,
        {
            "identifier",
            "eos_records",
            "evidence_kind",
            "notes",
            "coordinates",
            "row_selection",
            "phase_reviewed",
            "uncertainties",
            "fit_selection",
        },
        "view",
    )
    if not isinstance(view.get("identifier"), str) or not view["identifier"]:
        fail("identifier is required")
    if view.get("evidence_kind") not in EVIDENCE_KINDS:
        fail("an explicit evidence_kind is required")
    if not isinstance(view.get("notes"), str) or not view["notes"]:
        fail("source interpretation notes are required")
    links = view.get("eos_records")
    if (
        not isinstance(links, list)
        or not links
        or not all(isinstance(r, str) for r in links)
        or len(set(links)) != len(links)
    ):
        fail("eos_records must be a nonempty, unique array")
    if any(r not in record_ids for r in links):
        fail("unknown EOS record")
    columns = {c["name"]: c for c in dataset["columns"]}

    def column(name: Any, *, value: bool = True) -> Mapping[str, Any]:
        if not isinstance(name, str) or name not in columns:
            fail(f"unknown column {name!r}")
        c = columns[name]
        if value and c["role"] != "value":
            fail(f"coordinate {name!r} must have role value")
        return c

    def filters(items: Any) -> None:
        if not isinstance(items, list):
            fail("row selection must be an array")
        for item in items:
            object_keys(item, {"column", "values"}, "row selection")
            column(item.get("column"), value=False)
            if not isinstance(item.get("values"), list) or not item["values"]:
                fail("row selection needs explicit values")
            if any(
                not isinstance(v, (str, int, float)) or isinstance(v, bool)
                for v in item["values"]
            ):
                fail("row selection values must be strings or numbers")

    filters(view.get("row_selection", []))
    if not isinstance(view.get("phase_reviewed", False), bool):
        fail("phase_reviewed must be boolean")
    coordinates = view.get("coordinates")
    object_keys(coordinates, {"pressure", "volume", "temperature"}, "coordinates")
    pressure = coordinates.get("pressure", {})
    if not isinstance(pressure, Mapping):
        fail("pressure coordinate must be an object")
    kind = pressure.get("kind")
    if not isinstance(kind, str):
        fail("unsupported pressure coordinate")
    pressure_keys = {
        "column": {"kind", "column"},
        "ruby_scale": {
            "kind",
            "column",
            "source_calibration",
            "target_calibration",
            "reference",
        },
        "constant": {"kind", "value", "unit", "reference"},
        "calibrant_eos": {
            "kind",
            "volume_column",
            "row_selection",
            "pair_by",
            "eos_record",
            "temperature_k",
            "reference",
        },
    }
    object_keys(pressure, pressure_keys.get(kind, {"kind"}), "pressure coordinate")
    if kind in {"column", "ruby_scale"}:
        c = column(pressure.get("column"))
        if c["unit"] not in {"GPa", "Pa", "MPa", "kPa", "kbar", "bar", "Mbar"}:
            fail("pressure needs convertible units")
        if kind == "ruby_scale":
            from peritheos.pressure_calibrations import (
                RubyFluorescenceCalibration,
                get_pressure_calibration,
            )

            for key in ("source_calibration", "target_calibration"):
                if not isinstance(pressure.get(key), str) or not isinstance(
                    get_pressure_calibration(pressure[key]), RubyFluorescenceCalibration
                ):
                    fail(
                        "ruby reduction requires explicit executable ruby calibrations"
                    )
            if not pressure.get("reference"):
                fail("ruby reduction requires a source reference")
    elif kind == "constant":
        _positive_number(pressure.get("value"), allow_zero=True)
        if pressure.get("unit") != "GPa" or not pressure.get("reference"):
            fail("constant pressure needs GPa units and a source reference")
    elif kind == "calibrant_eos":
        if column(pressure.get("volume_column"))["unit"] != CELL_UNIT:
            fail("calibrant volumes must declare conventional-cell units")
        filters(pressure.get("row_selection", []))
        keys = pressure.get("pair_by")
        if not isinstance(keys, list) or not keys or len(set(keys)) != len(keys):
            fail("calibrant pairing needs unique explicit keys")
        for key in keys:
            column(key, value=False)
        if (
            not isinstance(pressure.get("eos_record"), str)
            or not pressure["eos_record"]
        ):
            fail("calibrant EOS record is required")
        _positive_number(pressure.get("temperature_k"))
        if not pressure.get("reference"):
            fail("calibrant reduction needs a source reference")
    else:
        fail("unsupported pressure coordinate")
    volume = coordinates.get("volume", {})
    if not isinstance(volume, Mapping):
        fail("volume coordinate must be an object")
    volume_keys = (
        {"kind", "column", "unit", "factor", "reciprocal", "reference"}
        if volume.get("kind") == "column"
        else {"kind", "columns", "constraints", "geometry", "unit", "reference"}
    )
    object_keys(volume, volume_keys, "volume coordinate")
    if volume.get("unit") != CELL_UNIT:
        fail("volume output must declare the conventional-cell basis")
    if volume.get("kind") == "column":
        c = column(volume.get("column"))
        factor = volume.get("factor", 1.0)
        _positive_number(factor)
        if not isinstance(volume.get("reciprocal", False), bool):
            fail("reciprocal must be boolean")
        if factor != 1 or volume.get("reciprocal") or c["unit"] != CELL_UNIT:
            if not volume.get("reference"):
                fail("a volume conversion needs a source/basis reference")
    elif volume.get("kind") == "lattice":
        supplied = object_keys(volume.get("columns"), set(_METRIC), "lattice columns")
        constraints = object_keys(
            volume.get("constraints"), set(_METRIC), "lattice constraints"
        )
        if set(supplied) - set(_METRIC) or set(constraints) - set(_METRIC):
            fail("unknown cell parameter")
        for key in _METRIC:
            if key in supplied:
                c = column(supplied[key])
                units = (
                    {"angstrom", "Å", "A", "angstroms", "nm", "pm"}
                    if key in _LENGTHS
                    else {"degree", "degrees", "deg", "°", "radian", "rad"}
                )
                if c["unit"] not in units:
                    fail(f"unsupported {key} units")
            if key in constraints:
                constraint = constraints[key]
                if isinstance(constraint, str):
                    # Deterministic, acyclic references to earlier same-dimension metrics.
                    allowed = _LENGTHS if key in _LENGTHS else _ANGLES
                    if constraint not in allowed[: allowed.index(key)]:
                        fail(
                            "cell constraint must reference an earlier compatible metric"
                        )
                else:
                    _positive_number(constraint)
                    if key in _ANGLES and constraint >= 180:
                        fail("cell angle must be less than 180 degrees")
            if key not in supplied and key not in constraints:
                fail(f"missing cell parameter {key}")
        if not volume.get("reference"):
            fail("lattice geometry needs a source/cell-setting reference")
    else:
        fail("unsupported volume coordinate")
    if "temperature" in coordinates:
        temperature = object_keys(
            coordinates["temperature"], {"column"}, "temperature coordinate"
        )
        if column(temperature.get("column"))["unit"] not in {"K", "degC"}:
            fail("temperature must have K or degC units")
    uncertainties = object_keys(
        view.get("uncertainties", {}), {"pressure", "volume"}, "uncertainties"
    )
    for axis, name in uncertainties.items():
        if axis not in {"pressure", "volume"}:
            fail("unknown uncertainty coordinate")
        c = column(name, value=False)
        coordinate = coordinates[axis]
        if (
            c["role"] not in {"standard_deviation", "standard_error", "uncertainty"}
            or coordinate["kind"] not in {"column", "ruby_scale"}
            or c.get("of") != coordinate["column"]
            or c["unit"] != columns[coordinate["column"]]["unit"]
        ):
            fail("uncertainty must describe its selected source coordinate and units")
    fit = object_keys(
        view.get("fit_selection", {}),
        {
            "flag_column",
            "row_selection",
            "pressure_above_gpa",
            "pressure_below_gpa",
            "notes",
        },
        "fit selection",
    )
    if "flag_column" in fit:
        column(fit["flag_column"], value=False)
    filters(fit.get("row_selection", []))
    for key in ("pressure_above_gpa", "pressure_below_gpa"):
        if key in fit:
            _positive_number(fit[key], allow_zero=True)


def _positive_number(value: Any, *, allow_zero: bool = False) -> float:
    if (
        isinstance(value, bool)
        or not isinstance(value, (float, int))
        or not np.isfinite(value)
        or (value < 0 if allow_zero else value <= 0)
    ):
        raise EosmatError("Observation coordinate requires a finite positive number")
    return float(value)


def validate_dataset_observation_metadata(document: Mapping[str, Any]) -> None:
    """Validate optional embedded views without requiring migration of old files."""
    records = {
        r["identifier"] for r in document.get("eos_records", []) if "identifier" in r
    }
    for dataset in document.get("datasets", []):
        views = dataset.get("observation_views", [])
        if not isinstance(views, list):
            raise EosmatError("observation_views must be an array")
        ids: set[str] = set()
        links: set[str] = set()
        for view in views:
            if not isinstance(view, Mapping):
                raise EosmatError("observation view must be an object")
            validate_observation_view(view, dataset, records)
            if view["identifier"] in ids or links.intersection(view["eos_records"]):
                raise EosmatError("Ambiguous observation view for an EOS/dataset pair")
            ids.add(view["identifier"])
            links.update(view["eos_records"])


def _reviewed_view(
    document: Mapping[str, Any],
    dataset: Mapping[str, Any],
    eos_record: str,
    registry: Mapping[str, Any] | None,
) -> dict[str, Any] | None:
    entries = registry if registry is not None else observation_registry()
    matching = [
        e
        for e in entries["datasets"]
        if e["material_id"] == document["identifier"]
        and e["dataset_id"] == dataset["identifier"]
    ]
    candidates = [
        copy.deepcopy(v)
        for v in dataset.get("observation_views", [])
        if eos_record in v["eos_records"]
    ]
    for entry in matching:
        views = [v for v in entry["views"] if eos_record in v["eos_records"]]
        if not views:
            continue
        record = next(
            r for r in document["eos_records"] if r["identifier"] == eos_record
        )
        binding = entry["record_pressure_calibration_sha256"][eos_record]
        if binding != observation_source_fingerprint(
            {"pressure_calibration": record.get("pressure_calibration")}
        ):
            raise DatasetError(
                "The reviewed observation view's pressure calibration changed; review the record binding"
            )
        if entry["declaration_sha256"] != observation_source_fingerprint(
            dataset
        ) or entry["material_identity_sha256"] != observation_material_fingerprint(
            document
        ):
            raise DatasetError(
                "Observation view source binding changed; review the declaration"
            )
        candidates.extend(copy.deepcopy(views))
    if len(candidates) > 1:
        raise DatasetError("Ambiguous observation view for this EOS/dataset pair")
    if candidates:
        validate_observation_view(
            candidates[0], dataset, {r["identifier"] for r in document["eos_records"]}
        )
        return candidates[0]
    return None


def _evidence_kind(dataset: Mapping[str, Any]) -> str | None:
    kind = dataset["kind"].lower()
    if any(
        k in kind
        for k in ("benchmark", "model", "calculated_eos", "parameter_grid", "eos_curve")
    ):
        return "model_output"
    if any(k in kind for k in ("simulation", "theoretical", "ab_initio")):
        return "simulation"
    if "digitized" in kind:
        return "digitized"
    if any(k in kind for k in ("shock", "hugoniot")):
        return "shock"
    # Only exact, typed experimental table kinds get an automatic interpretation.
    if kind in {
        "pressure_volume",
        "pressure_volume_temperature",
        "compression",
        "pvt",
        "pv",
        "experimental",
        "experimental_pv",
        "experimental_pvt",
    }:
        return "experimental"
    return None


def plan_dataset_observations(
    document: Mapping[str, Any],
    eos_record: str,
    dataset_identifier: str,
    *,
    registry: Mapping[str, Any] | None = None,
) -> dict[str, Any]:
    """Return a validated view or a specific unresolved reason, never guess rows.

    Automatic views are limited to unambiguous typed P/V columns in a single
    phase. They are returned as plans, not written over archival declarations.
    """
    records = {r["identifier"]: r for r in document.get("eos_records", [])}
    if eos_record not in records:
        raise DatasetError(f"Unknown EOS record {eos_record!r}")
    dataset = next(
        (
            d
            for d in document.get("datasets", [])
            if d["identifier"] == dataset_identifier
        ),
        None,
    )
    if dataset is None:
        raise DatasetError(f"Unknown dataset {dataset_identifier!r}")
    reviewed = _reviewed_view(document, dataset, eos_record, registry)
    if reviewed:
        return {"status": "ready", "origin": "declared", "view": reviewed}
    if eos_record not in dataset.get("used_by_eos_records", []):
        return {"status": "requires_review", "reason": "missing_record_link"}
    evidence = _evidence_kind(dataset)
    if evidence is None:
        return {"status": "requires_review", "reason": "undeclared_evidence_kind"}
    if any(
        c["quantity"]
        in {"phase", "phase_identifier", "phase_label", "material_identifier"}
        for c in dataset["columns"]
    ):
        return {
            "status": "requires_review",
            "reason": "sample_or_phase_selection_required",
        }
    columns = dataset["columns"]
    p = [
        c
        for c in columns
        if c["role"] == "value"
        and c["quantity"]
        in {
            "pressure",
            "reported_pressure",
            "source_pressure",
            "static_pressure",
            "hugoniot_pressure",
            "shock_pressure",
        }
    ]
    v = [
        c
        for c in columns
        if c["role"] == "value"
        and c["quantity"]
        in {
            "volume",
            "cell_volume",
            "unit_cell_volume",
            "conventional_unit_cell_volume",
            "conventional_cell_volume",
        }
        and c["unit"] == CELL_UNIT
    ]
    if len(p) != 1 or len(v) != 1:
        return {"status": "requires_review", "reason": "coordinate_binding_required"}
    view: dict[str, Any] = {
        "identifier": "typed_source_columns",
        "eos_records": [eos_record],
        "evidence_kind": evidence,
        "notes": "Unambiguous source-declared coordinates; no recalibration or fit.",
        "coordinates": {
            "pressure": {"kind": "column", "column": p[0]["name"]},
            "volume": {"kind": "column", "column": v[0]["name"], "unit": CELL_UNIT},
        },
    }
    temperatures = [
        c
        for c in columns
        if c["role"] == "value"
        and c["quantity"] == "temperature"
        and c["unit"] in {"K", "degC"}
    ]
    if len(temperatures) == 1:
        view["coordinates"]["temperature"] = {"column": temperatures[0]["name"]}
    elif len(temperatures) > 1:
        return {"status": "requires_review", "reason": "temperature_binding_required"}
    for axis, c in (("pressure", p[0]), ("volume", v[0])):
        errors = [
            x
            for x in columns
            if x.get("of") == c["name"]
            and x["role"] in {"standard_deviation", "standard_error", "uncertainty"}
        ]
        if len(errors) == 1:
            view.setdefault("uncertainties", {})[axis] = errors[0]["name"]
        elif len(errors) > 1:
            return {
                "status": "requires_review",
                "reason": "uncertainty_binding_required",
            }
    validate_observation_view(view, dataset, set(records))
    return {"status": "ready", "origin": "typed_columns", "view": view}


@dataclass(frozen=True)
class DatasetObservations:
    """Resolved display coordinates with original source rows and omissions.

    Uncertainties retain their declared role. Calibration and lattice errors
    are not propagated without a specified covariance model. Missing T is NaN.
    """

    dataset: Dataset
    view: Mapping[str, Any]
    row_indices: Any
    pressure_gpa: Any
    volume_a3: Any
    temperature_k: Any
    pressure_uncertainty: Any
    volume_uncertainty: Any
    pressure_uncertainty_role: str | None
    volume_uncertainty_role: str | None
    selected_for_fit: Any
    calibrant_row_indices: Any
    omissions: tuple[Mapping[str, Any], ...]
    provenance: Mapping[str, Any]


def resolve_dataset_observations(
    document: Mapping[str, Any],
    eos_record: str,
    dataset_identifier: str,
    *,
    resource_root: Any | None = None,
    registry: Mapping[str, Any] | None = None,
    calibration_records: Mapping[str, Any] | None = None,
) -> DatasetObservations:
    """Resolve a plan using the existing complete EOSRecord pressure evaluator.

    Pairing requires exact declared keys. Duplicate or missing calibrants cause
    explicit per-row omissions. Source temperatures never inherit EOS Tr.
    ``calibration_records`` permits an immutable catalog snapshot or test fixture.
    """
    plan = plan_dataset_observations(
        document, eos_record, dataset_identifier, registry=registry
    )
    if plan["status"] != "ready":
        raise DatasetError(f"Unresolved observation view: {plan['reason']}")
    metadata = next(
        d for d in document["datasets"] if d["identifier"] == dataset_identifier
    )
    dataset = Dataset.from_mapping(metadata, resource_root=resource_root)
    view = plan["view"]
    coordinates = view["coordinates"]
    p, v = coordinates["pressure"], coordinates["volume"]

    def numeric(column: str, unit: str | None = None) -> np.ndarray:
        # Source tokens such as '--' stay in Dataset; only the derived coordinate
        # treats them as missing. One placeholder must not hide the whole table.
        values = np.full(len(dataset), np.nan)
        for i, value in enumerate(dataset[column]):
            try:
                values[i] = float(value)
            except (TypeError, ValueError):
                pass
        meta = dataset.get_column(column)
        if unit is not None and unit != meta.unit:
            return _convert_values(
                values, meta.unit, unit, uncertainty=meta.role != "value"
            )
        return values

    def matches(index: int, filters: list[dict[str, Any]]) -> bool:
        return all(dataset[f["column"]][index] in f["values"] for f in filters)

    selected = [
        i for i in range(len(dataset)) if matches(i, view.get("row_selection", []))
    ]
    omissions: list[dict[str, Any]] = []
    pressure = np.full(len(dataset), np.nan)
    paired = np.full(len(dataset), -1, dtype=np.int64)
    calibration_provenance: dict[str, Any] = {}
    if p["kind"] in {"column", "ruby_scale"}:
        pressure = numeric(p["column"], "GPa")
        if p["kind"] == "ruby_scale":
            from peritheos.pressure_calibrations import recalculate_ruby_pressure

            finite = np.isfinite(pressure) & (pressure >= 0)
            reduced = np.full(len(dataset), np.nan)
            reduced[finite] = recalculate_ruby_pressure(
                pressure[finite], p["source_calibration"], p["target_calibration"]
            )
            pressure = reduced
            calibration_provenance = {
                "source_calibration": p["source_calibration"],
                "target_calibration": p["target_calibration"],
                "pressure_uncertainty": "not_propagated",
            }
    elif p["kind"] == "constant":
        pressure.fill(p["value"])
    else:
        from peritheos.catalog import get_eos_record

        calibration = (
            get_eos_record(p["eos_record"])
            if calibration_records is None
            else calibration_records[p["eos_record"]]
        )
        if calibration.identifier != p["eos_record"]:
            raise DatasetError("Calibrant record identifier does not match the view")
        # Source calibration identity must agree with the target record metadata.
        record = next(
            r for r in document["eos_records"] if r["identifier"] == eos_record
        )
        methods = record.get("pressure_calibration", {}).get("methods", [])
        if not any(
            m.get("reference_eos_record") == calibration.identifier for m in methods
        ):
            raise DatasetError(
                "Observation calibrant differs from EOS pressure calibration"
            )
        index: dict[tuple[Any, ...], list[int]] = {}

        def pair_key(i: int) -> tuple[Any, ...] | None:
            key = tuple(dataset[k][i] for k in p["pair_by"])
            if any(
                x is None
                or str(x).strip() == ""
                or (isinstance(x, (float, np.floating)) and not np.isfinite(x))
                for x in key
            ):
                return None
            return key

        for i in range(len(dataset)):
            if matches(i, p.get("row_selection", [])):
                calibrant_key = pair_key(i)
                if calibrant_key is not None:
                    index.setdefault(calibrant_key, []).append(i)
        volumes = numeric(p["volume_column"])
        cache: dict[int, float] = {}
        for i in selected:
            sample_key = pair_key(i)
            candidates = [] if sample_key is None else index.get(sample_key, [])
            if len(candidates) != 1:
                omissions.append(
                    {
                        "row_index": i,
                        "reason": "missing_calibrant"
                        if not candidates
                        else "ambiguous_calibrant",
                    }
                )
                continue
            j = candidates[0]
            if not np.isfinite(volumes[j]) or volumes[j] <= 0:
                omissions.append({"row_index": i, "reason": "invalid_calibrant_volume"})
                continue
            if j not in cache:
                cache[j] = float(calibration.pressure(volumes[j], p["temperature_k"]))
            pressure[i], paired[i] = cache[j], j
        calibration_provenance = {
            "eos_record": calibration.identifier,
            "temperature_k": p["temperature_k"],
            "record": copy.deepcopy(dict(calibration.eosmat_metadata)),
            "pressure_uncertainty": "not_propagated",
        }
    if v["kind"] == "column":
        original = numeric(v["column"])
        with np.errstate(divide="ignore", invalid="ignore"):
            volume = v.get("factor", 1.0) * (
                1 / original if v.get("reciprocal") else original
            )
    else:
        metrics: dict[str, np.ndarray] = {}
        length_factors = {
            "angstrom": 1,
            "Å": 1,
            "A": 1,
            "angstroms": 1,
            "nm": 10,
            "pm": 0.01,
        }
        for key in _METRIC:
            if key in v.get("columns", {}):
                name = v["columns"][key]
                metric_column = dataset.get_column(name)
                factor = (
                    length_factors[metric_column.unit]
                    if key in _LENGTHS
                    else (180 / np.pi if metric_column.unit in {"rad", "radian"} else 1)
                )
                metrics[key] = numeric(name) * factor
            else:
                value = v["constraints"][key]
                metrics[key] = (
                    metrics[value]
                    if isinstance(value, str)
                    else np.full(len(dataset), value)
                )
        a, b, c, alpha, beta, gamma = (metrics[k] for k in _METRIC)
        s = (alpha + beta + gamma) / 2
        terms = [s, s - alpha, s - beta, s - gamma]
        valid = np.logical_and.reduce([(t > 0) & (t < 180) for t in terms])
        valid &= np.logical_and.reduce(
            [
                np.isfinite(metrics[k])
                & (metrics[k] > 0)
                & ((metrics[k] < 180) if k in _ANGLES else True)
                for k in _METRIC
            ]
        )
        for key, constraint in v.get("constraints", {}).items():
            expected = (
                metrics[constraint] if isinstance(constraint, str) else constraint
            )
            valid &= np.isclose(metrics[key], expected, rtol=1e-8, atol=1e-8)
        with np.errstate(invalid="ignore"):
            volume = (
                2
                * a
                * b
                * c
                * np.sqrt(np.prod([np.sin(t * np.pi / 180) for t in terms], axis=0))
            )
        volume[~valid] = np.nan
    temperature = (
        numeric(coordinates["temperature"]["column"], "K").copy()
        if "temperature" in coordinates
        else np.full(len(dataset), np.nan)
    )
    temperature[~np.isfinite(temperature) | (temperature < 0)] = np.nan
    # Missing T never erases otherwise valid pressure/volume measurements.
    rejected = {o["row_index"] for o in omissions}
    keep: list[int] = []
    for i in selected:
        if i in rejected:
            continue
        if not np.isfinite(pressure[i]) or not np.isfinite(volume[i]) or volume[i] <= 0:
            omissions.append({"row_index": i, "reason": "invalid_pressure_or_volume"})
        else:
            keep.append(i)
    indices = np.asarray(keep, dtype=np.int64)

    def errors(axis: str) -> tuple[Any, str | None]:
        name = view.get("uncertainties", {}).get(axis)
        if name is None:
            return np.full(len(indices), np.nan), None
        if axis == "pressure" and p["kind"] == "ruby_scale":
            # A source uncertainty is not a target-scale uncertainty.
            return np.full(len(indices), np.nan), None
        column = dataset.get_column(name)
        values = numeric(name, "GPa" if axis == "pressure" else None).copy()
        if axis == "volume":
            values *= (
                np.abs(volume / original)
                if v.get("reciprocal")
                else v.get("factor", 1.0)
            )
        values[~np.isfinite(values) | (values < 0)] = np.nan
        return values[indices], column.role

    pressure_error, pressure_role = errors("pressure")
    volume_error, volume_role = errors("volume")
    fit = view.get("fit_selection")
    fit_mask = np.zeros(len(indices), dtype=bool)
    if fit is not None:
        for j, i in enumerate(indices):
            fit_mask[j] = (
                matches(i, fit.get("row_selection", []))
                and (
                    "pressure_above_gpa" not in fit
                    or pressure[i] > fit["pressure_above_gpa"]
                )
                and (
                    "pressure_below_gpa" not in fit
                    or pressure[i] < fit["pressure_below_gpa"]
                )
                and (
                    "flag_column" not in fit
                    or str(dataset[fit["flag_column"]][i])
                    in {"1", "1.0", "True", "true", "yes"}
                )
            )
    arrays = [
        indices,
        pressure[indices],
        volume[indices],
        temperature[indices],
        pressure_error,
        volume_error,
        fit_mask,
        paired[indices],
    ]
    for array in arrays:
        array.setflags(write=False)
    return DatasetObservations(
        dataset=dataset,
        view=copy.deepcopy(view),
        row_indices=arrays[0],
        pressure_gpa=arrays[1],
        volume_a3=arrays[2],
        temperature_k=arrays[3],
        pressure_uncertainty=arrays[4],
        volume_uncertainty=arrays[5],
        pressure_uncertainty_role=pressure_role,
        volume_uncertainty_role=volume_role,
        selected_for_fit=None if fit is None else arrays[6],
        calibrant_row_indices=arrays[7],
        omissions=tuple(omissions),
        provenance={
            "declaration_sha256": observation_source_fingerprint(metadata),
            "material_identity_sha256": observation_material_fingerprint(document),
            "resource_sha256": dataset.resource.sha256 if dataset.resource else None,
            "volume_unit": CELL_UNIT,
            "pressure_unit": "GPa",
            "temperature_policy": "source_column_or_unknown",
            "calibration": calibration_provenance,
        },
    )


__all__ = [
    "DatasetObservations",
    "observation_registry",
    "plan_dataset_observations",
    "resolve_dataset_observations",
]
