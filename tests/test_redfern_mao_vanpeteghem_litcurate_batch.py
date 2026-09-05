import csv
import hashlib
import importlib.util
from pathlib import Path

import pytest

from peritheos import Material, get_material_document

ROOT = Path(__file__).parents[1]


def _load_script(name):
    path = ROOT / "scripts" / name
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _record(material_identifier, record_identifier):
    document = get_material_document(material_identifier)
    matches = [
        record
        for record in document["eos_records"]
        if record["identifier"] == record_identifier
    ]
    assert len(matches) == 1
    return document, matches[0]


@pytest.mark.parametrize(
    ("record_identifier", "k0", "kp", "fixed"),
    [
        ("mg0991fe0008mn0001co3_redfern_1993_bm2_1", 142.0, 4.0, True),
        ("mg0991fe0008mn0001co3_redfern_1993_bm3_2", 151.0, 2.5, False),
    ],
)
def test_redfern_natural_magnesite_records_execute(record_identifier, k0, kp, fixed):
    document, record = _record("mg0991fe0008mn0001co3_magnesite", record_identifier)
    assert document["formula"] == "Mg0.991Fe0.008Mn0.001CO3"
    assert record["reference"]["doi"] == "10.1029/93GL02507"
    assert record["eos"]["parameters"] == {"V0": 279.4, "K0": k0, "K0_prime": kp}
    assert ("K0_prime" in record["fixed_parameters"]) is fixed
    assert record["volume_basis"]["molar_mass_g_mol"] == pytest.approx(84.595953)
    eos = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).get_eos_record(record_identifier)
    assert eos.pressure(279.4) == pytest.approx(0.0, abs=1e-12)
    assert eos.volume(eos.pressure(0.9 * 279.4)) == pytest.approx(
        0.9 * 279.4, rel=1e-10
    )


def test_redfern_reproduction_has_two_distinct_source_curves():
    result = _load_script("reproduce_redfern_1993_magnesite.py").reproduce()
    assert result["curves_gpa"]["fixed_kp"]["0.9"] == pytest.approx(18.47440617497)
    assert result["curves_gpa"]["free_kp"]["0.9"] == pytest.approx(18.03711940794)


@pytest.mark.parametrize(
    ("material_identifier", "record_identifier", "v0", "mass"),
    [
        ("bridgmanite", "bridgmanite_mao_1991_bm2_1", 162.49, 100.387),
        (
            "mg09fe01sio3_bridgmanite",
            "mg09fe01sio3_bridgmanite_mao_1991_bm2_2",
            162.79,
            103.541,
        ),
        (
            "mg080fe020sio3_bridgmanite",
            "mg080fe020sio3_bridgmanite_mao_1991_bm2_1",
            163.53,
            106.695,
        ),
    ],
)
def test_mao_1991_composition_records_execute(
    material_identifier, record_identifier, v0, mass
):
    document, record = _record(material_identifier, record_identifier)
    assert record["reference"]["doi"] == "10.1029/91JB00176"
    assert record["eos"]["parameters"] == {"V0": v0, "K0": 261.0, "K0_prime": 4.0}
    assert record["fixed_parameters"] == ["K0_prime"]
    assert record["volume_basis"]["molar_mass_g_mol"] == pytest.approx(mass)
    eos = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).get_eos_record(record_identifier)
    assert eos.pressure(v0) == pytest.approx(0.0, abs=1e-12)
    assert eos.pressure(0.9 * v0) == pytest.approx(33.9564789554)


def test_mao_reproduction_preserves_shared_modulus_and_three_volumes():
    result = _load_script("reproduce_mao_1991_bridgmanites.py").reproduce()
    assert result["shared_K0_gpa"] == 261.0
    assert set(result["curves_gpa"]) == {"MgSiO3", "Mg0.9Fe0.1SiO3", "Mg0.8Fe0.2SiO3"}


@pytest.mark.parametrize(
    ("record_identifier", "k0", "kp", "fixed"),
    [
        ("delta_alooh_vanpeteghem_2002_bm2_1", 252.0, 4.0, True),
        ("delta_alooh_vanpeteghem_2002_bm3_2", 228.0, 7.0, False),
    ],
)
def test_vanpeteghem_delta_alooh_records_execute(record_identifier, k0, kp, fixed):
    document, record = _record("delta_alooh", record_identifier)
    assert record["reference"]["doi"] == "10.1029/2001GL014224"
    assert record["eos"]["parameters"] == {"V0": 56.54, "K0": k0, "K0_prime": kp}
    assert ("K0_prime" in record["fixed_parameters"]) is fixed
    eos = Material.from_eosmat(
        document, record_identifiers=[record_identifier]
    ).get_eos_record(record_identifier)
    assert eos.pressure(56.54) == pytest.approx(0.0, abs=1e-12)
    assert eos.pressure(52.03) > 20.0


def test_vanpeteghem_table1_transcription_and_curve_diagnostics():
    resource = (
        ROOT / "peritheos/data/datasets/delta-alooh-vanpeteghem-2002-table1-pv.csv"
    )
    assert (
        hashlib.sha256(resource.read_bytes()).hexdigest()
        == "252a908bc05855868d8c933315d7d0afcd3485641f65b3597b139641dfaff826"
    )
    with resource.open(newline="", encoding="utf-8") as stream:
        rows = list(csv.DictReader(stream))
    assert len(rows) == 9
    assert rows[0]["pressure_gpa"] == "22.51"
    assert rows[-1]["volume_a3"] == "56.54"

    result = _load_script("reproduce_vanpeteghem_2002_delta_alooh.py").reproduce()
    assert result["observations"] == 9
    assert result["results"]["fixed_kp"]["pressure_rmse_gpa"] == pytest.approx(
        1.4332160943
    )
    assert result["results"]["free_kp"]["pressure_rmse_gpa"] == pytest.approx(
        1.3626960446
    )


def test_four_audit_documents_cover_all_twenty_two_litcurate_rows():
    expected = {
        "redfern-1993-magnesite.md": 4,
        "mao-1991-bridgmanites.md": 8,
        "vanpeteghem-2002-delta-alooh.md": 5,
        "xiao-2013-srsio3.md": 5,
    }
    for name, count in expected.items():
        text = (ROOT / "docs/literature-reproductions" / name).read_text(
            encoding="utf-8"
        )
        assert text.count("`litcurate_") == count
