"""New outputs remain independent of held evidence and scientifically scoped."""

import csv
import hashlib
import io
import json
from collections import Counter
from pathlib import Path

import numpy as np
import pytest
from scipy.constants import R

from peritheos import Dataset, get_material_document
from peritheos_build import ROOT, digest, stage_datasets
from scripts import generate_campbell_2009_pressure_integral as integration
from scripts import generate_sakai_2025_lineage as lineage

PROOF = json.loads((ROOT / "docs/data/independent-lineage-buffers.json").read_bytes())


def rows(payload):
    return list(csv.DictReader(io.StringIO(payload.decode())))


def test_lineage_generation_reads_only_bound_primary_resources(monkeypatch):
    config = json.loads(lineage.SOURCE.read_bytes())
    allowed = {lineage.SOURCE} | {
        lineage.DATA / entry["file"] for entry in config["sources"].values()
    }
    read = Path.read_bytes

    def guarded(path):
        assert path in allowed, f"Undeclared/held input read: {path}"
        return read(path)

    expected = lineage.OUTPUT.read_bytes()
    monkeypatch.setattr(Path, "read_bytes", guarded)
    assert lineage.generate() == expected


def test_lineage_distinguishes_gauges_from_co_measurement_and_missing_evidence():
    graph = rows(lineage.OUTPUT.read_bytes())
    assert len(graph) == 37
    assert all(int(r["observations"]) > 0 for r in graph)
    fei = [r for r in graph if r["source_doi"] == "10.1002/2016GL069456"]
    assert {r["calibrant_or_partner"]: int(r["observations"]) for r in fei} == {
        "neon": 49,
        "mgo": 18,
        "nacl_b2": 7,
    }
    assert all(r["relation"] == "reported_pressure_gauge" for r in fei)
    assert not any(
        r["source_doi"] == "10.1016/j.pepi.2013.12.010"
        and r["calibrant_or_partner"] == "mgo"
        for r in graph
    )
    present_study = [
        r for r in graph if r["source_doi"] == "10.1038/s43246-025-00792-5"
    ]
    assert len(present_study) == 7
    assert sum(int(r["observations"]) for r in present_study) == 393
    assert all(r["relation"] == "paired_volumes" for r in present_study)
    assert all(
        r["relation"] == "co_measured"
        for r in graph
        if r["source_doi"] == "10.1029/2012JB009292"
    )
    assert not any("ruby2020" in r["calibrant_or_partner"] for r in graph)
    assert len({tuple(r.values()) for r in graph}) == 37


def test_lineage_input_change_fails_closed(tmp_path):
    config = json.loads(lineage.SOURCE.read_bytes())
    config["sources"]["sakai_s1"]["sha256"] = "0" * 64
    path = tmp_path / "modified.json"
    path.write_text(json.dumps(config))
    with pytest.raises(ValueError, match="Changed lineage input"):
        lineage.generate(source=path)


def test_pressure_integrator_units_and_quadrature_have_analytic_reference(monkeypatch):
    monkeypatch.setattr(integration, "volume", lambda p, t, c: c[0] + c[1] * p)
    pressures = np.array([0, 1, 10, 100])
    actual = integration.integrate(pressures, 1500, [1, 0.02], [3, -0.01])
    expected = (2 * pressures - 0.015 * pressures**2) * 0.8686 * 1000 / (R * 1500)
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)
    for p, t in [([1, 2], 1000), ([0, 0], 1000), ([0, np.nan], 1000), ([0, 1], 0)]:
        with pytest.raises(ValueError):
            integration.integrate(p, t, [1, 0.02], [3, -0.01])


@pytest.fixture(scope="module")
def independent_nno():
    return integration.generate()


def test_nno_generation_never_reads_original_output_or_polynomial(monkeypatch):
    expected = integration.OUTPUT.read_bytes()
    read = Path.read_bytes

    def guarded(path):
        assert path == integration.SOURCE, f"Non-parameter input read: {path}"
        return read(path)

    monkeypatch.setattr(Path, "read_bytes", guarded)
    np.testing.assert_allclose(
        [
            float(r["delta_log10_fugacity_from_zero_gpa"])
            for r in rows(integration.generate())
        ],
        [float(r["delta_log10_fugacity_from_zero_gpa"]) for r in rows(expected)],
        rtol=0,
        atol=1e-9,
    )


def test_nno_pressure_only_scope_convergence_and_source_mismatch(independent_nno):
    generated = rows(independent_nno)
    assert len(generated) == 1717
    assert set(generated[0]) == {
        "temperature_k",
        "pressure_gpa",
        "delta_log10_fugacity_from_zero_gpa",
    }
    assert Counter(r["temperature_k"] for r in generated) == {
        str(t): 101 for t in range(1000, 2601, 100)
    }
    assert all(
        float(r["delta_log10_fugacity_from_zero_gpa"]) == 0
        for r in generated
        if r["pressure_gpa"] == "0"
    )
    configuration = json.loads(integration.SOURCE.read_bytes())
    for temperature in [1000, 2600]:
        args = (
            [0, 1, 10, 60, 100],
            temperature,
            configuration["metal_parameters"],
            configuration["oxide_parameters"],
        )
        np.testing.assert_allclose(
            integration.integrate(*args, order=8),
            integration.integrate(*args, order=16),
            rtol=0,
            atol=2e-9,
        )
    comparison = integration.compare_original(independent_nno)
    assert comparison["maximum_absolute_difference"] == pytest.approx(
        0.01425912761639836, abs=1e-8
    )
    assert comparison["rmse_delta_log10_fugacity"] == pytest.approx(
        0.005282584020164474, abs=1e-8
    )
    assert PROOF["campbell"]["new_absolute_tables"] == 0
    assert PROOF["campbell"]["original_table_hold_retained"] is True


def test_only_three_new_declarations_change_the_catalog(
    independent_generation_baseline,
):
    assert independent_generation_baseline["proof"] == PROOF
    for key, fingerprint in PROOF["added_dataset_declarations"].items():
        material, identifier = key.split("/")
        declaration = next(
            d
            for d in get_material_document(material)["datasets"]
            if d["identifier"] == identifier
        )
        assert digest(declaration) == fingerprint
        assert declaration["used_by_eos_records"] == []
        assert len(Dataset.from_mapping(declaration)) == (
            37 if material == "gold" else 1717
        )
    for path, fingerprint in PROOF["source_bindings"].items():
        assert hashlib.sha256((ROOT / path).read_bytes()).hexdigest() == fingerprint


def test_installed_graph_and_pressure_output_do_not_admit_held_originals(tmp_path):
    package = stage_datasets(tmp_path / "datasets") / "peritheos_datasets"
    review = json.loads((package / "redistribution-review.json").read_bytes())
    manifest = json.loads((package / "manifest.json").read_bytes())
    assert manifest["dataset_count"] == 451
    assert review["independent_lineage_buffers"] == PROOF
    assert len([e for e in review["entries"] if e["decision"] == "hold"]) == 5
    for name in [
        "sakai-2025-table-s8-calibration-lineage.csv",
        "campbell-2009-tables-s4-s5-buffer-results.csv",
    ]:
        assert not (package / "data/datasets" / name).exists()
    for source in [lineage.OUTPUT, integration.OUTPUT]:
        assert (
            package / "data/datasets" / source.name
        ).read_bytes() == source.read_bytes()
    assert (
        package / "notices/independent-lineage-buffers.LICENSE.md"
    ).read_bytes() == (
        ROOT / "peritheos/data/datasets/independent-lineage-buffers.LICENSE.md"
    ).read_bytes()
