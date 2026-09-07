"""Typed access to observation datasets stored in ``.eosmat`` documents.

The loader intentionally uses only the standard library and NumPy.  Dataset
columns remain generic; domain-specific views such as pressure--volume data are
small conveniences built on top of the column metadata.
"""

from __future__ import annotations

import copy
import csv
import hashlib
import io
from collections.abc import Mapping, Sequence
from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path, PurePosixPath
from types import MappingProxyType
from typing import Any

import numpy as np
from numpy.typing import NDArray

from peritheos.errors import DatasetError, DatasetLookupError

_UNCERTAINTY_ROLES = (
    "standard_deviation",
    "standard_error",
    "uncertainty",
)
_PRESSURE_QUANTITIES = (
    "pressure",
    "hugoniot_pressure",
    "shock_pressure",
    "reported_pressure",
    "calculated_pressure",
    "source_pressure",
    "static_pressure",
)
_VOLUME_QUANTITIES = (
    "conventional_unit_cell_volume",
    "unit_cell_volume",
    "cell_volume",
    "conventional_cell_volume",
    "volume",
    "formula_unit_volume",
    "atomic_volume",
    "molar_volume",
    "specific_volume",
    "relative_volume",
    "volume_ratio",
)


def _read_only(values: Sequence[Any], *, textual: bool = False) -> NDArray[Any]:
    if textual:
        array = np.asarray(values, dtype=object)
    else:
        try:
            array = np.asarray(values, dtype=float)
        except (TypeError, ValueError):
            array = np.asarray(values, dtype=object)
    array.setflags(write=False)
    return array


def _plain_mapping(value: Mapping[str, Any]) -> Mapping[str, Any]:
    return MappingProxyType(copy.deepcopy(dict(value)))


@dataclass(frozen=True)
class DatasetColumn:
    """Schema metadata for one dataset column."""

    name: str
    quantity: str
    unit: str
    role: str
    of: str | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", _plain_mapping(self.metadata))


@dataclass(frozen=True)
class DatasetResource:
    """Provenance for a packaged dataset resource."""

    path: str
    sha256: str
    media_type: str
    metadata: Mapping[str, Any] = field(default_factory=dict, compare=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "metadata", _plain_mapping(self.metadata))


@dataclass(frozen=True)
class PressureVolumeData:
    """A pressure--volume view with converted values and uncertainties.

    ``pressure_uncertainty`` and ``volume_uncertainty`` retain the declared
    uncertainty kind in their companion ``*_uncertainty_role`` attributes.
    The ``*_sigma`` aliases are populated only for columns explicitly declared
    as standard deviations.
    """

    dataset: Dataset
    pressure: NDArray[np.float64] = field(repr=False, compare=False)
    volume: NDArray[np.float64] = field(repr=False, compare=False)
    pressure_unit: str
    volume_unit: str
    pressure_uncertainty: NDArray[np.float64] | None = field(
        default=None, repr=False, compare=False
    )
    volume_uncertainty: NDArray[np.float64] | None = field(
        default=None, repr=False, compare=False
    )
    pressure_uncertainty_role: str | None = None
    volume_uncertainty_role: str | None = None

    @property
    def pressure_sigma(self) -> NDArray[np.float64] | None:
        """Pressure standard deviations, or ``None`` when not declared."""
        if self.pressure_uncertainty_role == "standard_deviation":
            return self.pressure_uncertainty
        return None

    @property
    def volume_sigma(self) -> NDArray[np.float64] | None:
        """Volume standard deviations, or ``None`` when not declared."""
        if self.volume_uncertainty_role == "standard_deviation":
            return self.volume_uncertainty
        return None


@dataclass(frozen=True)
class Dataset:
    """An immutable, loaded observation table and its source metadata."""

    identifier: str
    kind: str
    source_location: str
    reference: str | Mapping[str, Any]
    columns: tuple[DatasetColumn, ...]
    data: Mapping[str, NDArray[Any]] = field(repr=False, compare=False)
    description: str | None = None
    used_by_eos_records: tuple[str, ...] = ()
    notes: str | tuple[str, ...] | None = None
    resource: DatasetResource | None = None
    checksum_verified: bool | None = None
    source_column_names: tuple[str, ...] | None = None
    metadata: Mapping[str, Any] = field(default_factory=dict, compare=False, repr=False)

    def __post_init__(self) -> None:
        object.__setattr__(self, "reference", copy.deepcopy(self.reference))
        object.__setattr__(self, "data", MappingProxyType(dict(self.data)))
        object.__setattr__(self, "metadata", _plain_mapping(self.metadata))

    @classmethod
    def from_mapping(
        cls,
        metadata: Mapping[str, Any],
        *,
        resource_root: Any | None = None,
    ) -> Dataset:
        """Load a dataset mapping, including checksum verification for resources.

        ``resource_root`` defaults to the packaged :mod:`peritheos.data`
        directory.  A path or ``importlib.resources`` traversable may be passed
        when loading a compatible external document.
        """
        try:
            raw_columns = metadata["columns"]
            columns = tuple(
                DatasetColumn(
                    name=str(raw["name"]),
                    quantity=str(raw["quantity"]),
                    unit=str(raw["unit"]),
                    role=str(raw["role"]),
                    of=None if raw.get("of") is None else str(raw["of"]),
                    metadata={
                        key: copy.deepcopy(value)
                        for key, value in raw.items()
                        if key not in {"name", "quantity", "unit", "role", "of"}
                    },
                )
                for raw in raw_columns
            )
        except (KeyError, TypeError) as error:
            raise DatasetError(f"Invalid dataset column metadata: {error}") from error

        names = tuple(column.name for column in columns)
        if len(names) != len(set(names)):
            raise DatasetError("Dataset column names must be unique")

        has_rows = "rows" in metadata
        has_resource = "resource" in metadata
        if has_rows == has_resource:
            raise DatasetError("Dataset must contain exactly one of rows or resource")

        resource_info: DatasetResource | None = None
        checksum_verified: bool | None = None
        if has_rows:
            rows = metadata["rows"]
            checksum_verified = None
            source_column_names = None
        else:
            raw_resource = metadata["resource"]
            try:
                resource_info = DatasetResource(
                    path=str(raw_resource["path"]),
                    sha256=str(raw_resource["sha256"]),
                    media_type=str(raw_resource["media_type"]),
                    metadata={
                        key: copy.deepcopy(value)
                        for key, value in raw_resource.items()
                        if key not in {"path", "sha256", "media_type"}
                    },
                )
            except (KeyError, TypeError) as error:
                raise DatasetError(
                    f"Invalid dataset resource metadata: {error}"
                ) from error
            rows, source_column_names = _load_resource_rows(
                resource_info, columns, resource_root
            )
            checksum_verified = True
        if not isinstance(rows, Sequence) or isinstance(rows, (str, bytes)):
            raise DatasetError("Dataset rows must be a sequence")
        by_column: dict[str, list[Any]] = {name: [] for name in names}
        for index, row in enumerate(rows):
            if len(row) != len(columns):
                raise DatasetError(
                    f"Dataset row {index} has {len(row)} values; expected {len(columns)}"
                )
            for column, value in zip(columns, row):
                by_column[column.name].append(value)

        arrays = {
            column.name: _read_only(
                by_column[column.name],
                textual=column.unit in {"text", "identifier"},
            )
            for column in columns
        }
        known = {
            "identifier",
            "kind",
            "source_location",
            "reference",
            "columns",
            "rows",
            "resource",
            "description",
            "used_by_eos_records",
            "notes",
        }
        try:
            return cls(
                identifier=str(metadata["identifier"]),
                kind=str(metadata["kind"]),
                source_location=str(metadata["source_location"]),
                reference=copy.deepcopy(metadata["reference"]),
                columns=columns,
                data=arrays,
                description=(
                    None
                    if metadata.get("description") is None
                    else str(metadata["description"])
                ),
                used_by_eos_records=tuple(metadata.get("used_by_eos_records", ())),
                notes=copy.deepcopy(metadata.get("notes")),
                resource=resource_info,
                checksum_verified=checksum_verified,
                source_column_names=source_column_names,
                metadata={
                    key: copy.deepcopy(value)
                    for key, value in metadata.items()
                    if key not in known
                },
            )
        except KeyError as error:
            raise DatasetError(f"Missing dataset metadata: {error}") from error

    def __len__(self) -> int:
        """Return the number of observations."""
        return len(next(iter(self.data.values()), ()))

    def __getitem__(self, name: str) -> NDArray[Any]:
        """Return one column by name."""
        try:
            return self.data[name]
        except KeyError as error:
            raise DatasetLookupError(
                f"Unknown column {name!r} in dataset {self.identifier!r}; "
                f"available: {list(self.data)}",
                operation="lookup_dataset_column",
                field="name",
                context={"dataset": self.identifier, "column": name},
            ) from error

    def get_column(self, name: str) -> DatasetColumn:
        """Return typed metadata for one column by name."""
        for column in self.columns:
            if column.name == name:
                return column
        self[name]
        raise AssertionError("unreachable")

    def find_columns(
        self, *, quantity: str | None = None, role: str | None = None
    ) -> tuple[DatasetColumn, ...]:
        """Find columns by semantic quantity and/or role."""
        return tuple(
            column
            for column in self.columns
            if (quantity is None or column.quantity == quantity)
            and (role is None or column.role == role)
        )

    def values(self, name: str, *, unit: str | None = None) -> NDArray[Any]:
        """Return a column, optionally converted to a compatible unit."""
        column = self.get_column(name)
        values = self[name]
        if unit is None or unit == column.unit:
            return values
        if values.dtype.kind not in "iuf":
            raise DatasetError(
                f"Column {name!r} is not numeric and cannot be converted"
            )
        return _convert_values(
            values, column.unit, unit, uncertainty=column.role != "value"
        )

    def as_pressure_volume(
        self,
        *,
        pressure_unit: str | None = None,
        volume_unit: str | None = None,
        pressure_column: str | None = None,
        volume_column: str | None = None,
    ) -> PressureVolumeData:
        """Return the table's primary pressure and volume columns.

        Explicit column names resolve datasets containing multiple sample or
        calibrant pressure/volume series.  Associated uncertainty columns are
        selected through their schema ``of`` relationship.
        """
        pressure_meta = self._select_value_column(
            pressure_column, _PRESSURE_QUANTITIES, "pressure"
        )
        volume_meta = self._select_value_column(
            volume_column, _VOLUME_QUANTITIES, "volume"
        )
        requested_pressure_unit = pressure_unit or pressure_meta.unit
        requested_volume_unit = volume_unit or volume_meta.unit
        pressure = self.values(pressure_meta.name, unit=requested_pressure_unit)
        volume = self.values(volume_meta.name, unit=requested_volume_unit)
        pressure_error = self._uncertainty_for(pressure_meta)
        volume_error = self._uncertainty_for(volume_meta)
        return PressureVolumeData(
            dataset=self,
            pressure=pressure,
            volume=volume,
            pressure_unit=requested_pressure_unit,
            volume_unit=requested_volume_unit,
            pressure_uncertainty=(
                None
                if pressure_error is None
                else self.values(pressure_error.name, unit=requested_pressure_unit)
            ),
            volume_uncertainty=(
                None
                if volume_error is None
                else self.values(volume_error.name, unit=requested_volume_unit)
            ),
            pressure_uncertainty_role=(
                None if pressure_error is None else pressure_error.role
            ),
            volume_uncertainty_role=(
                None if volume_error is None else volume_error.role
            ),
        )

    def _select_value_column(
        self,
        explicit_name: str | None,
        preferred_quantities: tuple[str, ...],
        label: str,
    ) -> DatasetColumn:
        if explicit_name is not None:
            column = self.get_column(explicit_name)
            if column.role != "value":
                raise DatasetError(f"{label}_column must select a value column")
            return column
        for quantity in preferred_quantities:
            candidates = self.find_columns(quantity=quantity, role="value")
            if len(candidates) == 1:
                return candidates[0]
            if len(candidates) > 1:
                names = [column.name for column in candidates]
                raise DatasetError(
                    f"Dataset {self.identifier!r} has multiple {label} columns "
                    f"{names}; pass {label}_column explicitly"
                )
        raise DatasetError(
            f"Dataset {self.identifier!r} has no recognized {label} value column"
        )

    def _uncertainty_for(self, value: DatasetColumn) -> DatasetColumn | None:
        for role in _UNCERTAINTY_ROLES:
            candidates = tuple(
                column
                for column in self.columns
                if column.of == value.name and column.role == role
            )
            if len(candidates) > 1:
                raise DatasetError(
                    f"Dataset {self.identifier!r} has multiple {role} columns "
                    f"for {value.name!r}"
                )
            if candidates:
                return candidates[0]
        return None


def _load_resource_rows(
    resource: DatasetResource,
    columns: tuple[DatasetColumn, ...],
    resource_root: Any | None,
) -> tuple[list[list[Any]], tuple[str, ...]]:
    path = PurePosixPath(resource.path)
    if path.is_absolute() or ".." in path.parts:
        raise DatasetError(f"Dataset resource path {resource.path!r} is not local")
    root = resources.files("peritheos.data") if resource_root is None else resource_root
    resource_path: Any
    if isinstance(root, (str, Path)):
        resource_path = Path(root).joinpath(*path.parts)
    else:
        resource_path = root.joinpath(*path.parts)
    try:
        payload = resource_path.read_bytes()
    except (OSError, FileNotFoundError) as error:
        raise DatasetError(
            f"Could not read dataset resource {resource.path!r}: {error}"
        ) from error
    actual = hashlib.sha256(payload).hexdigest()
    if actual != resource.sha256:
        raise DatasetError(
            f"SHA-256 mismatch for dataset resource {resource.path!r}: "
            f"expected {resource.sha256}, got {actual}",
            code="dataset.checksum_mismatch",
            operation="load_dataset",
            field="resource.sha256",
            context={
                "path": resource.path,
                "expected": resource.sha256,
                "actual": actual,
            },
        )
    if resource.media_type != "text/csv":
        raise DatasetError(
            f"Unsupported dataset media type {resource.media_type!r}",
            code="dataset.unsupported_media_type",
            operation="load_dataset",
            field="resource.media_type",
        )
    reader = csv.reader(io.StringIO(payload.decode("utf-8-sig")))
    try:
        header = tuple(next(reader))
    except StopIteration as error:
        raise DatasetError(f"Dataset resource {resource.path!r} is empty") from error
    if len(header) != len(columns):
        raise DatasetError(
            f"Dataset resource {resource.path!r} header has {len(header)} columns; "
            f"expected {len(columns)} from its metadata"
        )
    rows: list[list[Any]] = []
    for row_index, row in enumerate(reader, start=2):
        if len(row) != len(header):
            raise DatasetError(
                f"Dataset resource {resource.path!r} row {row_index} has "
                f"{len(row)} values; expected {len(header)}"
            )
        rows.append(
            [
                _parse_csv_value(value, textual=column.unit in {"text", "identifier"})
                for column, value in zip(columns, row)
            ]
        )
    return rows, header


def _parse_csv_value(value: str, *, textual: bool) -> Any:
    stripped = value.strip()
    if not stripped:
        return None if textual else np.nan
    if textual:
        return value
    try:
        return float(stripped)
    except ValueError:
        return value


_LINEAR_UNITS = {
    "Pa": ("pressure", 1.0),
    "kPa": ("pressure", 1.0e3),
    "MPa": ("pressure", 1.0e6),
    "GPa": ("pressure", 1.0e9),
    "bar": ("pressure", 1.0e5),
    "kbar": ("pressure", 1.0e8),
    "Mbar": ("pressure", 1.0e11),
    "m": ("length", 1.0),
    "mm": ("length", 1.0e-3),
    "nm": ("length", 1.0e-9),
    "pm": ("length", 1.0e-12),
    "angstrom": ("length", 1.0e-10),
    "m/s": ("velocity", 1.0),
    "km/s": ("velocity", 1.0e3),
    "kg/m^3": ("density", 1.0),
    "Mg/m^3": ("density", 1.0e3),
    "g/cm^3": ("density", 1.0e3),
    "1": ("dimensionless", 1.0),
    "dimensionless": ("dimensionless", 1.0),
}
_VOLUME_BASES = {
    "m^3": 1.0,
    "cm^3": 1.0e-6,
    "nm^3": 1.0e-27,
    "angstrom^3": 1.0e-30,
}


def _unit_definition(unit: str) -> tuple[str, str, float]:
    if unit in _LINEAR_UNITS:
        dimension, scale = _LINEAR_UNITS[unit]
        return dimension, "", scale
    for base, scale in _VOLUME_BASES.items():
        if unit == base:
            return "volume", "", scale
        prefix = base + "/"
        if unit.startswith(prefix):
            return "volume", unit[len(prefix) :], scale
    raise DatasetError(f"Unit {unit!r} is not supported for conversion")


def _convert_values(
    values: NDArray[Any], source: str, target: str, *, uncertainty: bool
) -> NDArray[np.float64]:
    if source in {"K", "degC"} and target in {"K", "degC"}:
        result = np.asarray(values, dtype=float).copy()
        if not uncertainty:
            if source == "degC" and target == "K":
                result += 273.15
            elif source == "K" and target == "degC":
                result -= 273.15
    else:
        source_dimension, source_basis, source_scale = _unit_definition(source)
        target_dimension, target_basis, target_scale = _unit_definition(target)
        if (source_dimension, source_basis) != (target_dimension, target_basis):
            raise DatasetError(
                f"Cannot convert {source!r} to incompatible unit {target!r}"
            )
        result = np.asarray(values, dtype=float) * (source_scale / target_scale)
    result.setflags(write=False)
    return result


__all__ = [
    "Dataset",
    "DatasetColumn",
    "DatasetResource",
    "PressureVolumeData",
]
