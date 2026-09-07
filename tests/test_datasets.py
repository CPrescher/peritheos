import hashlib

import numpy as np
import pytest

from peritheos import (
    Dataset,
    DatasetError,
    DatasetLookupError,
    PressureVolumeData,
    get_material,
)


def test_packaged_coesite_dataset_loads_and_converts_pressure():
    material = get_material("coesite")

    dataset = material.get_dataset("coesite_levien_1981_table7_pv")
    pressure_volume = dataset.as_pressure_volume(pressure_unit="GPa")

    assert isinstance(dataset, Dataset)
    assert isinstance(pressure_volume, PressureVolumeData)
    assert dataset.checksum_verified is True
    assert dataset.resource is not None
    assert dataset.resource.path == "datasets/coesite-levien-1981-table7-pv.csv"
    assert dataset.source_column_names == tuple(
        column.name for column in dataset.columns
    )
    assert dataset.source_location == "Table 7, page 330"
    assert dataset.used_by_eos_records == ("coesite_levien_1981_bm3_1",)
    assert len(dataset) == 9
    np.testing.assert_allclose(pressure_volume.pressure[:3], [0.0001, 2.18, 2.24])
    np.testing.assert_allclose(pressure_volume.volume[:3], [546.46, 535.47, 535.1])
    np.testing.assert_allclose(pressure_volume.pressure_sigma[:3], [np.nan, 0.05, 0.05])
    np.testing.assert_allclose(pressure_volume.volume_sigma[:3], [0.05, 0.07, 0.2])
    assert pressure_volume.pressure_unit == "GPa"
    assert pressure_volume.volume_unit == "angstrom^3"


def test_raw_material_datasets_remain_available():
    material = get_material("coesite")

    assert isinstance(material.datasets[0], dict) is False
    assert material.datasets[0]["identifier"] == "coesite_levien_1981_table7_pv"
    assert material.datasets[0]["resource"]["sha256"] == (
        "ba8c7fcdd377eeb28f24a3ba0586884120159b6ecbf9eb6d87caa2c59e75f1cc"
    )


def test_embedded_rows_use_the_same_typed_interface():
    dataset = get_material("b4c").get_dataset("b4c_somayazulu_2023_table4_pvt")

    pressure_volume = dataset.as_pressure_volume()

    assert dataset.resource is None
    assert dataset.checksum_verified is None
    assert len(dataset) == 51
    np.testing.assert_allclose(pressure_volume.pressure[:2], [10.1, 11.9])
    np.testing.assert_allclose(pressure_volume.volume[:2], [314.0, 314.7])
    assert dataset.get_column("temperature_k").quantity == "temperature"
    np.testing.assert_allclose(dataset["temperature_k"][:2], [300.0, 1027.0])
    with pytest.raises(ValueError):
        dataset["temperature_k"][0] = 1.0


def test_find_columns_and_lookup_errors_are_typed():
    material = get_material("coesite")
    dataset = material.get_dataset("coesite_levien_1981_table7_pv")

    assert [column.name for column in dataset.find_columns(quantity="pressure")] == [
        "pressure_kbar",
        "pressure_sigma_kbar",
    ]
    with pytest.raises(DatasetLookupError, match="Unknown column"):
        dataset["not_a_column"]
    with pytest.raises(DatasetLookupError, match="Unknown dataset"):
        material.get_dataset("not_a_dataset")


def test_resource_checksum_is_verified_before_csv_parsing(tmp_path):
    payload = b"pressure_gpa,volume_a3\n1,10\n"
    (tmp_path / "sample.csv").write_bytes(payload)
    metadata = {
        "identifier": "sample",
        "kind": "pressure_volume",
        "source_location": "test",
        "reference": "test",
        "columns": [
            {
                "name": "pressure_gpa",
                "quantity": "pressure",
                "unit": "GPa",
                "role": "value",
            },
            {
                "name": "volume_a3",
                "quantity": "volume",
                "unit": "angstrom^3",
                "role": "value",
            },
        ],
        "resource": {
            "path": "sample.csv",
            "sha256": hashlib.sha256(payload).hexdigest(),
            "media_type": "text/csv",
        },
    }

    assert len(Dataset.from_mapping(metadata, resource_root=tmp_path)) == 1
    metadata["resource"]["sha256"] = "0" * 64
    with pytest.raises(DatasetError) as error:
        Dataset.from_mapping(metadata, resource_root=tmp_path)
    assert error.value.code == "dataset.checksum_mismatch"


def test_pressure_volume_column_override_and_incompatible_units():
    dataset = Dataset.from_mapping(
        {
            "identifier": "multiple",
            "kind": "pressure_volume",
            "source_location": "test",
            "reference": "test",
            "columns": [
                {
                    "name": "sample_pressure",
                    "quantity": "pressure",
                    "unit": "GPa",
                    "role": "value",
                },
                {
                    "name": "calibrant_pressure",
                    "quantity": "pressure",
                    "unit": "GPa",
                    "role": "value",
                },
                {
                    "name": "volume",
                    "quantity": "volume",
                    "unit": "angstrom^3",
                    "role": "value",
                },
            ],
            "rows": [[1.0, 1.1, 10.0]],
        }
    )

    with pytest.raises(DatasetError, match="multiple pressure columns"):
        dataset.as_pressure_volume()
    view = dataset.as_pressure_volume(pressure_column="sample_pressure")
    np.testing.assert_allclose(view.pressure, [1.0])
    with pytest.raises(DatasetError, match="incompatible"):
        dataset.values("sample_pressure", unit="angstrom")


def test_general_columns_and_unit_conversion_edges():
    dataset = Dataset.from_mapping(
        {
            "identifier": "general",
            "kind": "mixed",
            "source_location": "test",
            "reference": "test",
            "columns": [
                {
                    "name": "temperature_c",
                    "quantity": "temperature",
                    "unit": "degC",
                    "role": "value",
                },
                {
                    "name": "temperature_sigma_c",
                    "quantity": "temperature",
                    "unit": "degC",
                    "role": "standard_deviation",
                    "of": "temperature_c",
                },
                {
                    "name": "label",
                    "quantity": "observation_identifier",
                    "unit": "text",
                    "role": "flag",
                },
                {
                    "name": "pressure",
                    "quantity": "pressure",
                    "unit": "GPa",
                    "role": "value",
                },
                {
                    "name": "volume",
                    "quantity": "volume",
                    "unit": "angstrom^3",
                    "role": "value",
                },
            ],
            "rows": [[26.85, 2.0, "001", 1.0, 1000.0]],
        }
    )

    np.testing.assert_allclose(dataset.values("temperature_c", unit="K"), [300.0])
    np.testing.assert_allclose(dataset.values("temperature_sigma_c", unit="K"), [2.0])
    np.testing.assert_allclose(dataset.values("volume", unit="nm^3"), [1.0])
    assert dataset["label"].tolist() == ["001"]
    with pytest.raises(DatasetError, match="not numeric"):
        dataset.values("label", unit="identifier")
    with pytest.raises(DatasetError, match="not supported"):
        dataset.values("pressure", unit="psi")

    pressure_volume = dataset.as_pressure_volume()
    assert pressure_volume.pressure_sigma is None
    assert pressure_volume.volume_sigma is None


@pytest.mark.parametrize(
    ("change", "message"),
    [
        ({"columns": [{}]}, "Invalid dataset column metadata"),
        ({"rows": "not rows"}, "rows must be a sequence"),
        ({"rows": [[1.0, 2.0]]}, "expected 1"),
        ({"resource": {}}, "exactly one"),
    ],
)
def test_invalid_embedded_dataset_metadata(change, message):
    metadata = {
        "identifier": "invalid",
        "kind": "test",
        "source_location": "test",
        "reference": "test",
        "columns": [{"name": "x", "quantity": "x", "unit": "1", "role": "value"}],
        "rows": [[1.0]],
    }
    metadata.update(change)

    with pytest.raises(DatasetError, match=message):
        Dataset.from_mapping(metadata)


def test_duplicate_and_non_value_columns_are_rejected():
    metadata = {
        "identifier": "invalid",
        "kind": "test",
        "source_location": "test",
        "reference": "test",
        "columns": [
            {"name": "x", "quantity": "pressure", "unit": "GPa", "role": "value"},
            {"name": "x", "quantity": "volume", "unit": "angstrom^3", "role": "value"},
        ],
        "rows": [[1.0, 2.0]],
    }
    with pytest.raises(DatasetError, match="must be unique"):
        Dataset.from_mapping(metadata)

    metadata["columns"][1]["name"] = "volume_error"
    metadata["columns"][1]["role"] = "standard_deviation"
    dataset = Dataset.from_mapping(metadata)
    with pytest.raises(DatasetError, match="must select a value column"):
        dataset.as_pressure_volume(volume_column="volume_error")
    with pytest.raises(DatasetError, match="no recognized volume"):
        dataset.as_pressure_volume()


def test_resource_structure_failures_are_reported(tmp_path):
    column = {"name": "x", "quantity": "x", "unit": "1", "role": "value"}

    def resource_metadata(path, payload, *, media_type="text/csv"):
        return {
            "identifier": "resource",
            "kind": "test",
            "source_location": "test",
            "reference": "test",
            "columns": [column],
            "resource": {
                "path": path,
                "sha256": hashlib.sha256(payload).hexdigest(),
                "media_type": media_type,
            },
        }

    with pytest.raises(DatasetError, match="not local"):
        Dataset.from_mapping(
            resource_metadata("../sample.csv", b""), resource_root=tmp_path
        )
    with pytest.raises(DatasetError, match="Could not read"):
        Dataset.from_mapping(
            resource_metadata("missing.csv", b""), resource_root=tmp_path
        )

    payload = b"x\n1\n"
    (tmp_path / "sample.csv").write_bytes(payload)
    with pytest.raises(DatasetError, match="Unsupported dataset media type"):
        Dataset.from_mapping(
            resource_metadata("sample.csv", payload, media_type="application/json"),
            resource_root=tmp_path,
        )

    for filename, invalid_payload, message in (
        ("empty.csv", b"", "is empty"),
        ("header.csv", b"x,y\n", "header has 2 columns"),
        ("row.csv", b"x\n1,2\n", "row 2 has 2 values"),
    ):
        (tmp_path / filename).write_bytes(invalid_payload)
        with pytest.raises(DatasetError, match=message):
            Dataset.from_mapping(
                resource_metadata(filename, invalid_payload), resource_root=tmp_path
            )
