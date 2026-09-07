import csv
import json
import subprocess
import sys
from pathlib import Path

from peritheos import recalculate_eos_pressure_scale
from peritheos.eos import rt, thermal

ROOT = Path(__file__).resolve().parents[1]


def test_validation_narrative_matches_primary_refit_ledger():
    summary = json.loads(
        (ROOT / "docs/data/primary-eos-refits.json").read_text(encoding="utf-8")
    )["summary"]
    validation = " ".join(
        (ROOT / "docs/validation.md").read_text(encoding="utf-8").split()
    )

    assert f"{summary['total']} records" in validation
    assert f"{summary['parity']} uncertainty-parity matches" in validation
    assert f"{summary['similar']} additional numerically similar results" in validation
    assert f"[{summary['parity_not_achieved']} direct refits]" in validation
    assert f"while {summary['not_refittable']} records cannot be" in validation


def test_material_readme_inventory_matches_bundled_datasets():
    data_root = ROOT / "peritheos/data"
    datasets: dict[str, int] = {}
    material_links = 0
    record_identifiers: set[str] = set()

    for path in sorted((data_root / "materials").glob("*.eosmat")):
        document = json.loads(path.read_text(encoding="utf-8"))
        for dataset in document.get("datasets", []):
            material_links += 1
            record_identifiers.update(dataset.get("used_by_eos_records", []))
            identifier = dataset["identifier"]
            if identifier in datasets:
                continue
            if "rows" in dataset:
                row_count = len(dataset["rows"])
            else:
                resource = data_root / dataset["resource"]["path"]
                with resource.open(encoding="utf-8", newline="") as handle:
                    row_count = sum(1 for _ in csv.reader(handle)) - 1
            datasets[identifier] = row_count

    readme = " ".join(
        (data_root / "materials/README.md").read_text(encoding="utf-8").split()
    )
    expected = (
        f"contains {len(datasets)} distinct primary datasets with "
        f"{sum(datasets.values()):,} observation rows, represented by "
        f"{material_links} material-document links to "
        f"{len(record_identifiers)} EOS records"
    )
    assert expected in readme


def test_source_audit_index_is_complete_and_current():
    subprocess.run(
        [sys.executable, "scripts/generate_source_audit_index.py", "--check"],
        cwd=ROOT,
        check=True,
    )


def test_reference_dois_are_clickable():
    references = (ROOT / "docs/references.md").read_text(encoding="utf-8")
    for line in references.splitlines():
        if "doi:" in line.casefold():
            assert "https://doi.org/" in line.casefold()


def test_api_reference_lists_all_public_equation_classes():
    api_reference = (ROOT / "docs/api.md").read_text(encoding="utf-8")
    for name in (*rt.__all__, *thermal.__all__):
        assert f"    {name}," in api_reference or f"`{name}`" in api_reference


def test_strict_pressure_normalization_example_executes():
    normalized = recalculate_eos_pressure_scale(
        source_eos_record="tungsten_dewaele_2004_vinet_2",
        sample_volume=[31.0, 29.0, 27.0],
        target_standard_eos_record="gold_sokolova_2013_holzapfel_4",
        sample_temperature_k=298.0,
        standard_temperature_k=298.0,
        check_validity=True,
    )

    assert normalized.source_node == "tungsten_dewaele_2004_vinet_2"
    assert normalized.target_node == "gold_sokolova_2013_holzapfel_4"
    assert len(normalized.target_pressure_gpa) == 3
