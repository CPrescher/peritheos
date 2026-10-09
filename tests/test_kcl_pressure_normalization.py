"""KCl observation pressure contracts and conditional refit capability."""

import json
from pathlib import Path

import jsonschema
import numpy as np
import pytest

from peritheos import (
    DatasetError,
    Material,
    eosmat_schema,
    get_eos_record_document,
    get_material_document,
    recalculate_pressure_calibration_path,
    resolve_dataset_pressure,
)
from scripts.audit_kcl_pressure_normalization import CASES, OUTPUT, SOURCE, audit


@pytest.mark.parametrize("material,record_id,dataset_id", CASES)
def test_dewaele_all_rows_have_executable_pressure_contract(
    material, record_id, dataset_id
):
    document = get_material_document(material)
    jsonschema.validate(document, eosmat_schema())
    selected = Material.from_eosmat(
        document, record_identifiers=[record_id]
    ).to_eosmat()
    assert (
        selected["datasets"][0]["pressure_reductions"]
        == next(d for d in document["datasets"] if d["identifier"] == dataset_id)[
            "pressure_reductions"
        ]
    )
    resolved = resolve_dataset_pressure(selected, record_id, dataset_id)
    expected_count = 123 if material == "kcl" else 30
    assert len(resolved.pressure_gpa) == expected_count
    np.testing.assert_array_equal(resolved.row_indices, np.arange(expected_count))
    np.testing.assert_array_equal(resolved.pressure_gpa, resolved.reported_pressure)
    assert resolved.dataset.checksum_verified
    assert resolved.source_calibration["identifier"] == SOURCE
    assert resolved.uncertainty_treatment == "not_propagated"
    assert (
        get_eos_record_document(record_id)["pressure_calibration"]["recalculation"][
            "status"
        ]
        == "ready"
    )

    # Independent algebra: solve 1884*x*(1+5.5*x)=P, then apply Mao.
    pressure = resolved.pressure_gpa
    shift = 2 * pressure / 1884 / (1 + np.sqrt(1 + 22 * pressure / 1884))
    expected = 1904 / 7.665 * ((1 + shift) ** 7.665 - 1)
    graph = recalculate_pressure_calibration_path(
        pressure, record_id, "ruby_mao_1986", 298.0
    )
    np.testing.assert_allclose(graph.target_pressure_gpa, expected, atol=1e-11)


def test_diagnostic_refits_are_reproducible_and_keep_source_constraints():
    report = audit()
    assert report == json.loads(OUTPUT.read_text(encoding="utf-8"))
    for record in report["records"]:
        published = get_eos_record_document(record["record_identifier"])
        assert record["original_published_parameters"] == published["eos"]["parameters"]
        for reduction in record["reductions"].values():
            fit = reduction["conditional_unweighted_vinet_refit"]
            for name, value in fit["fixed_parameters"].items():
                assert fit["parameters_conventional_cell"][name] == value
            assert fit["uncertainty_treatment"] == "not_propagated"
        assert (
            record["reductions"][SOURCE]["conditional_unweighted_vinet_refit"]
            != record["reductions"]["ruby_mao_1986"][
                "conditional_unweighted_vinet_refit"
            ]
        )
    assert (
        len(
            report["records"][0]["reductions"]["ruby_mao_1986"][
                "pressure_extrapolated_row_indices"
            ]
        )
        == 23
    )


@pytest.mark.parametrize(
    "record_id,dataset_id",
    [
        ("kcl_walker_2002_bm3_2", "kcl_walker_2002_table2_pvt"),
    ],
)
def test_unverified_exact_scales_do_not_authorize_a_reduction(record_id, dataset_id):
    document = get_material_document("kcl")
    with pytest.raises(DatasetError, match="Unresolved"):
        resolve_dataset_pressure(document, record_id, dataset_id)
    calibration = get_eos_record_document(record_id)["pressure_calibration"]
    assert calibration["methods"][0]["reference"]["year"] == 1986
    assert calibration["recalculation"]["status"] == "reference_eos_not_bundled"


def test_manifest_recalculation_counts_match_current_records():
    from collections import Counter

    root = Path(__file__).resolve().parents[1] / "peritheos/data/materials"
    counts = Counter(
        record["pressure_calibration"]["recalculation"]["status"]
        for path in root.glob("*.eosmat")
        for record in json.loads(path.read_text(encoding="utf-8"))["eos_records"]
    )
    assert json.loads((root / "manifest.json").read_text(encoding="utf-8"))[
        "pressure_calibration"
    ]["recalculation_counts"] == dict(counts)
