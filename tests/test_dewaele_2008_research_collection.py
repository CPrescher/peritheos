"""Research calculations remain preserved and cannot leak into catalog discovery."""

import hashlib
import json
from importlib.resources import files
from pathlib import Path

import pytest

from peritheos import (
    Dataset,
    Material,
    get_eos_record,
    get_material_document,
    list_eos_record_documents,
    list_eos_records,
    load_eosmat,
    search_eos_records,
)
from scripts.reproduce_dewaele_2008_metals import Z, table4_specs, table6_specs

ROOT = Path(__file__).resolve().parents[1]
COLLECTION = ROOT / "peritheos/data/research/dewaele-2008-paw"


def test_paw_records_are_absent_from_all_catalog_surfaces():
    archived_ids = {rid for _, rid, _, _, _ in table6_specs()}
    assert len(archived_ids) == 54
    assert not archived_ids.intersection(list_eos_record_documents())
    assert not archived_ids.intersection(
        r.identifier for r in list_eos_records(include_reference_variants=True)
    )
    assert not archived_ids.intersection(
        r.identifier for r in search_eos_records(include_reference_variants=True)
    )
    for material in Z:
        document = get_material_document(material)
        assert not any(r["identifier"] in archived_ids for r in document["eos_records"])
        assert not any(
            d["identifier"] == f"{material}_dewaele_2008_figure4"
            for d in document.get("datasets", [])
        )
    for rid in archived_ids:
        with pytest.raises(KeyError):
            get_eos_record(rid)
    experimental = {rid for _, rid, _, _, _ in table4_specs()}
    assert len(experimental) == 13
    assert experimental.issubset(list_eos_record_documents())
    assert all(get_eos_record(rid) for rid in experimental)


def test_native_research_files_preserve_source_coefficients_and_evidence():
    documents = {m: load_eosmat(COLLECTION / f"{m}.eosmat") for m in Z}
    assert sum(len(d["eos_records"]) for d in documents.values()) == 54
    for material, rid, q, functional, state in table6_specs():
        r = next(
            r for r in documents[material]["eos_records"] if r["identifier"] == rid
        )
        expected = [q[0] * Z[material], q[1], q[2]]
        assert [
            r["parameters"][k]["value"] for k in ("V0", "K0", "K0_prime")
        ] == expected
        assert all(p["error"] is None for p in r["parameters"].values())
        assert r["temperature_ref"] == (0 if state == "static" else 298)
        assert r["reference"]["doi"] == "10.1103/PhysRevB.78.104102"
        evidence = r["evidence"]
        assert evidence["parameter_provenance"]
        assert evidence["scientific_validation"]["primary_source_check"][
            "source_pdf_sha256"
        ]
        reproduction = evidence["reproduction"]
        assert (
            reproduction["source_equation_check"]["independent_paw_refit_available"]
            is False
        )
        assert reproduction["assessment"]["scope"] == "source_equation"
        assert reproduction["previous_catalog_report"]
        assert reproduction["assessment"]["report"] == "reproduction.json"
        assert (
            "docs/research/dewaele-2008-paw.md#record-"
            in reproduction["assessment"]["report_path"]
        )


def test_collection_is_packaged_self_contained_and_explicitly_loadable():
    collection = files("peritheos.data").joinpath("research/dewaele-2008-paw")
    manifest = json.loads(collection.joinpath("manifest.json").read_text())
    assert manifest["catalog_scope"] == "research_only"
    assert manifest["eos_records"] == 54
    assert len(manifest["material_files"]) == 9
    for name, checksum in manifest["files"].items():
        assert (
            hashlib.sha256(collection.joinpath(name).read_bytes()).hexdigest()
            == checksum
        )
    for name in manifest["material_files"]:
        document = load_eosmat(COLLECTION / name)
        ids = {r["identifier"] for r in document["eos_records"]}
        for metadata in document["datasets"]:
            dataset = Dataset.from_mapping(metadata, resource_root=collection)
            assert dataset.checksum_verified
            assert len(dataset) > 0
            assert set(dataset.used_by_eos_records) == ids
    silver = Material.from_eosmat(load_eosmat(COLLECTION / "silver.eosmat"))
    assert len(silver.eos_records) == 6
    assert (
        silver.eos_records[0].pressure(0.8 * silver.eos_records[0].reference_volume) > 0
    )
    # Explicit loading has no side effect on database registration.
    assert not any("_dewaele_2008_paw_" in r.identifier for r in list_eos_records())


def test_archived_audits_and_all_research_report_anchors_are_retained():
    ids = {rid for _, rid, _, _, _ in table6_specs()}
    reconstruction = json.loads((COLLECTION / "reproduction.json").read_text())
    source = json.loads((COLLECTION / "source-audit.json").read_text())
    assert {r["record_identifier"] for r in reconstruction["records"]} == ids
    assert {r["record"] for r in source["records"]} == ids
    page = (ROOT / "docs/research/dewaele-2008-paw.md").read_text()
    for rid in ids:
        assert f'<a id="record-{rid}"></a>' in page
    assert sum(r["status"] == "parity" for r in reconstruction["records"]) == 50
    assert sum(r["status"] == "deferred" for r in source["records"]) == 2
    assert sum(r["status"] == "not_reproduced" for r in source["records"]) == 2
