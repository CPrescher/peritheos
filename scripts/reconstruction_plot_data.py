"""Replay stored curves on recoverable audit rows, without running any optimizer.

These repository-only report adapters reuse the audit's unit/row readers. They
refuse mismatched row counts or recorded pressure RMS rather than silently
plotting a different experiment as the reconstruction's input.
"""

from __future__ import annotations

import copy
import csv
import json
from dataclasses import dataclass, replace
from pathlib import Path
from typing import Any

import numpy as np

from peritheos import Material, convert_eosmat
from scripts.reconstruction_diagnostics import (
    DiagnosticComparison,
    diagnostic_comparison,
)
from scripts.validate_primary_eos_refits import (
    SHEN_SELECTIONS,
    Series,
    _masked_series,
    _series,
    _shen_series,
)

ROOT = Path(__file__).resolve().parents[1]


@dataclass
class CurveComparison:
    source: Any
    refit: Any
    series: Series
    source_pressure: np.ndarray
    refit_pressure: np.ndarray
    source_rmse: float
    refit_rmse: float
    volume_unit: str
    independent_refit: bool = False
    note: str = ""


def executable(document, record, parameters, side):
    """Build a private plotting copy; preserve all unrefitted fixed coefficients."""
    candidate = copy.deepcopy(document)
    target = next(
        r for r in candidate["eos_records"] if r["identifier"] == record["identifier"]
    )
    for row in parameters:
        name = row["parameter"]
        if name.startswith("eos."):
            name = "rt_eos." + name.removeprefix("eos.")
        elif name.startswith("thermal."):
            name = name.removeprefix("thermal.")
        value = row.get(side)
        if value is None:
            raise ValueError(f"{side} coefficient {name} is missing")
        if (
            name == "alpha0"
            and record["identifier"]
            == "superhydrous_phase_b_lt_litasov_2007_thermal_pressure_a89"
        ):
            continue  # Convert the audit's alpha0*K0 after all elastic updates.
        if name.startswith("rt_eos."):
            block, key = target["eos"]["parameters"], name.removeprefix("rt_eos.")
        elif name in target["eos"]["parameters"]:
            block, key = target["eos"]["parameters"], name
        elif name in target.get("thermal", {}).get("parameters", {}):
            block, key = target["thermal"]["parameters"], name
        else:
            raise ValueError(f"Audit coefficient {name} has no direct model mapping")
        if key not in block:
            raise ValueError(f"Audit coefficient {name} has no direct model mapping")
        block[key] = value
    if (
        record["identifier"]
        == "superhydrous_phase_b_lt_litasov_2007_thermal_pressure_a89"
    ):
        alpha0 = next(row[side] for row in parameters if row["parameter"] == "alpha0")
        target["thermal"]["parameters"]["alpha_KT_ref"] = (
            alpha0 * target["eos"]["parameters"]["K0"]
        )
    material = Material.from_eosmat(
        candidate,
        require_primary_validation=False,
        record_identifiers=[record["identifier"]],
    )
    return next(r for r in material.eos_records if r.identifier == record["identifier"])


def observations(document, record, fit):
    datasets = {d["identifier"]: d for d in document.get("datasets", [])}
    identifiers = fit.get("dataset_identifiers", [])
    rid = record["identifier"]

    def direct(
        pressure, volume, temperature=None, pressure_sigma=None, volume_sigma=None
    ):
        return Series(
            "+".join(identifiers),
            np.asarray(pressure, dtype=float),
            np.asarray(volume, dtype=float),
            None if temperature is None else np.asarray(temperature, dtype=float),
            pressure_sigma,
            volume_sigma,
            None,
            "audit pressure (GPa)",
            "audit volume in record units",
            "temperature (K)" if temperature is not None else None,
            fit.get("selection", "Named source-audit row reader"),
        )

    def table_rows(filename):
        with (ROOT / "peritheos/data/datasets" / filename).open(
            newline="", encoding="utf-8"
        ) as stream:
            return list(csv.DictReader(stream))

    def column(rows, name):
        return np.array([float(row[name]) for row in rows])

    molar_to_cell = float(document.get("formula_units_per_cell") or 1) / 0.602214076
    if rid == "argon_fcc_finger_1981_murnaghan2_debye":
        from peritheos.eos.finger1981 import CELL_TO_MOLAR

        rows = table_rows("argon-fcc-finger-1981-table1.csv")
        return direct(
            column(rows, "pressure_gpa"),
            column(rows, "molar_volume_cm3_mol") / CELL_TO_MOLAR,
            np.full(len(rows), 293.0),
        )
    if "_dewaele_2008_table4_" in rid:
        from scripts.reproduce_dewaele_2008_metals import experiment_rows

        rows = experiment_rows(document["identifier"])
        pressure_name = (
            "pressure_revised_gpa"
            if document["identifier"] == "fe"
            else "ruby_pressure_dewaele_gpa"
        )
        return direct(
            column(rows, pressure_name),
            column(rows, "atomic_volume_a3") * document["formula_units_per_cell"],
        )
    if fit.get("fit_kind") == "finger_1979_fixed_measured_v0":
        rows = table_rows("ferrous-ringwoodite-finger-1979-table2.csv")[1:]
        return direct(
            column(rows, "pressure_kbar") / 10,
            column(rows, "volume_a3"),
            pressure_sigma=column(rows, "pressure_esd_kbar") / 10,
            volume_sigma=column(rows, "volume_esd_a3"),
        )
    if fit.get("fit_kind") == "partial_298k_rydberg_stacey_diagnostic":
        phase = (
            "bridgmanite"
            if rid.startswith("bridgmanite_")
            else "post_perovskite"
            if rid.startswith("mgsio3_")
            else "akimotoite"
        )
        phases = (
            {"akimotoite_wang", "akimotoite_ruby"} if phase == "akimotoite" else {phase}
        )
        rows = [
            r
            for r in table_rows("dorogokupets-2015-298k-literature-reconstruction.csv")
            if r["phase"] in phases
        ]
        return direct(column(rows, "pressure_gpa"), column(rows, "volume_a3"))
    if rid in {
        "akimotoite_wang_2004_t0133_bm3",
        "akimotoite_wang_2004_t0150_bm3",
        "mgsio3_post_perovskite_komabayashi_2008_thermal_bm3",
    }:
        from scripts.reproduce_wang_zhou_komabayashi_mgsio3 import (
            observations as source_rows,
        )

        data = np.asarray(source_rows()[rid], dtype=float)
        return direct(
            data[:, 2], data[:, 1], data[:, 3] if record.get("thermal") else None
        )
    if rid == "ca_perovskite_noguchi_2013_bm2_mgd_1":
        from scripts.reproduce_noguchi_2013_casio3 import (
            SOURCE_TABLE_PATH,
            load_source_table,
        )

        rows = [r for r in load_source_table(SOURCE_TABLE_PATH) if r["fit_included"]]
        return direct(
            column(rows, "pressure_fei_gpa"),
            column(rows, "volume_a3"),
            column(rows, "temperature_k"),
        )
    if rid == "diamond_datchi_2007_vinet_1":
        rows = [
            r
            for r in table_rows("diamond-occelli-2003-figure2-digitized.csv")
            if r["best_supported_table_ii_selection"].lower() == "true"
        ]
        return direct(
            column(rows, "pressure_h2005_gpa"),
            column(rows, "atomic_volume_a3")
            * float(document["formula_units_per_cell"]),
        )
    if rid == "nickel_oxide_noguchi_1999_bm3_1":
        from scripts.reproduce_noguchi_1999_nio import reduce_table1

        reduced = reduce_table1()
        return replace(
            direct(
                reduced.isothermal_pressure_gpa,
                reduced.volume_ratio * record["eos"]["parameters"]["V0"],
            ),
            pressure_column="Hugoniot pressures thermally reduced to 300 K; HEL precursors excluded",
        )
    if "_ye_2017_" in rid:
        from scripts.reproduce_zhu_2025_pressure_standards import (
            ANCHORS,
            _load,
            _vinet_pressure,
        )

        material = (
            "Au"
            if rid.startswith("gold_")
            else "MgO"
            if rid.startswith("mgo_")
            else "Pt"
        )
        data = _load(
            "ye_2017_data_s2_pt_mgo_300k"
            if material == "Pt"
            else "ye_2017_data_s1_au_mgo_300k"
        )
        standard = "Au" if material == "MgO" else "MgO"
        anchor = dict(ANCHORS[standard])
        if material == "Pt":
            # Stored preceding MgO-from-Au stage, not a new optimization.
            from scripts.reconstruction_diagnostics import ledger

            anchor["K0_prime"] = next(
                p["refit"]
                for p in ledger()["mgo_ye_2017_vinet_300k"]["parameters"]
                if p["parameter"] == "K0_prime"
            )
        return replace(
            direct(
                _vinet_pressure(data[standard.lower() + "_volume_a3"], anchor),
                data[material.lower() + "_volume_a3"],
            ),
            pressure_column="Pressure from paired calibrant volumes, following the stored Ye stage",
        )
    if fit.get("fit_kind") == "static_vinet_with_explicit_ruby_scale_conversion":
        from scripts.audit_dewaele_2019_static_dac import (
            SOURCES,
            _load_rows,
            _target_pressure,
        )

        spec = next((s for s in SOURCES if rid.startswith(s.record_prefix + "_")), None)
        if spec is not None:
            pressure, volume, _ = _load_rows(spec)
            scale = "mao" if "_mao_" in rid else "dor"
            return direct(
                _target_pressure(pressure, spec.source_scale, scale),
                volume * float(document["formula_units_per_cell"]),
            )
    if rid.startswith("brucite_fukui_2003_") and not record.get("thermal"):
        from scripts.brucite_reproduce import observations as brucite_rows

        rows = [
            r
            for r in brucite_rows("fukui_2003", "table4")
            if float(r["temperature_k"]) == record["temperature_ref"]
        ]
        return direct(
            column(rows, "pressure_gpa"),
            column(rows, "cell_volume_a3"),
            pressure_sigma=column(rows, "pressure_error_gpa"),
            volume_sigma=column(rows, "cell_volume_error_a3"),
        )
    if rid in {
        "platinum_fei_2007_vinet_mgd",
        "neon_fcc_fei_2007_bm3_1",
        "neon_fcc_fei_2007_vinet_2",
    }:
        if rid.startswith("platinum"):
            from scripts.reproduce_fei_2007_platinum import observations as fei_rows
        else:
            from scripts.reproduce_fei_2007_neon import observations as fei_rows
        volume, temperature, pressure, _, _ = fei_rows()
        return direct(pressure, volume, temperature)
    if rid in {
        "iron_miozzi_2020_bm3",
        "iron_miozzi_2020_vinet",
        "iron_miozzi_2020_bm3_mgd",
    }:
        from scripts.reproduce_miozzi_2020_iron import observations as miozzi_rows

        data = miozzi_rows()
        mask = (
            np.ones(len(data["pressure_gpa"]), dtype=bool)
            if record.get("thermal")
            else data["temperature_k"] <= 300
        )
        if not record.get("thermal"):
            # The named effective-error diagnostic excludes the one 300 K row
            # without a positive effective coordinate error.
            mask &= (np.nan_to_num(data["pressure_error_gpa"]) > 0) | (
                np.nan_to_num(data["volume_error_a3"]) > 0
            )
        return direct(
            data["pressure_gpa"][mask],
            data["volume_a3"][mask],
            data["temperature_k"][mask] if record.get("thermal") else None,
        )
    if rid in {"gold_hirose_2008_bm3_fit1", "gold_hirose_2008_bm3_fit2"}:
        rows = [
            r
            for r in table_rows("gold-hirose-2008-table1.csv")
            if float(r["temperature_k"]) > 300
        ]
        if rid.endswith("fit2"):
            rows = table_rows("gold-fei-2004-table1.csv") + rows
        return direct(
            column(rows, "pressure_gpa"),
            column(rows, "volume_a3"),
            column(rows, "temperature_k"),
        )
    if rid == "argon_fcc_dewaele_2021_vinet_mgd":
        rows = [
            r
            for r in table_rows("argon-fcc-dewaele-2021-supplement.csv")
            if r["selected_for_room_temperature_fit"] == "1"
        ]
        return direct(
            column(rows, "pressure_gpa"),
            column(rows, "argon_a_angstrom") ** 3,
            np.full(len(rows), 296.0),
        )
    if rid.startswith("kcl_b2_tateno_2019_") and (
        "_holmes_" in rid or "_sokolova_" in rid
    ):
        from scripts.reproduce_kcl_variants import vinet

        rows = table_rows("kcl-tateno-2019-official-table-s1.csv")
        temperature = column(rows, "temperature_k")
        pressure = column(rows, "pressure_gpa")
        if "_holmes_" in rid:
            pressure = vinet(
                column(rows, "platinum_unit_cell_volume_a3"),
                60.4000884,
                798.31 / 3,
                1 + 7.2119 / 1.5,
            ) + 0.0069426 * (temperature - 300)
        return direct(pressure, column(rows, "kcl_unit_cell_volume_a3"), temperature)
    if rid.startswith("klb1_"):
        from scripts.audit_ricolleau_2009_klb1_eos import _read_numeric

        data = _read_numeric()
        phase = (
            "mg_perovskite"
            if "mg_perovskite" in rid
            else "ca_perovskite"
            if "ca_perovskite" in rid
            else "ferropericlase"
        )
        mask = np.isfinite(data[phase + "_volume_a3"])
        if "low_spin_300k" in rid:
            mask &= (data["temperature_k"] == 300) & (data["pressure_gpa"] >= 70)
        elif phase == "ferropericlase":
            mask &= data["temperature_k"] > 300
        return direct(
            data["pressure_gpa"][mask],
            data[phase + "_volume_a3"][mask],
            data["temperature_k"][mask] if record.get("thermal") else None,
            pressure_sigma=data["pressure_sigma_gpa"][mask],
            volume_sigma=data[phase + "_volume_sigma_a3"][mask],
        )
    if "_campbell_2009_" in rid:
        from scripts.reproduce_campbell_2009_buffers import current_data

        material = rid.split("_campbell_2009_", 1)[0]
        volume, temperature, pressure, sv, _, sp = current_data(material)
        return direct(
            pressure,
            volume * molar_to_cell,
            temperature,
            pressure_sigma=sp,
            volume_sigma=sv * molar_to_cell,
        )
    if "_farla_2025_" in rid:
        from scripts.reproduce_rubidium_halides_farla_2025 import read_rows

        compound = rid.split("_", 1)[0]
        rows = [
            r
            for r in read_rows(compound, "observations")
            if r["volume_a3"]
            and (r["source_table"] == "S1" or float(r["temperature_k"]) > 300)
        ]
        return direct(
            column(rows, "pressure_gpa"),
            column(rows, "volume_a3"),
            column(rows, "temperature_k"),
        )
    if "_chidester_2018_" in rid:
        material = "thorianite" if rid.startswith("thorianite") else "tho2_cotunnite"
        rows = [
            r
            for r in table_rows(f"chidester_2018_{material}_pvt.csv")
            if r["included_in_fit"] == "1"
        ]
        return direct(
            column(rows, "pressure_gpa"),
            column(rows, "molar_volume_cm3_mol") * molar_to_cell,
            column(rows, "temperature_k"),
            pressure_sigma=column(rows, "pressure_error_gpa"),
            volume_sigma=column(rows, "molar_volume_error_cm3_mol") * molar_to_cell,
        )
    if rid in {
        "forsterite_katsura_2022_author_script_bm3_mgd",
        "forsterite_katsura_2022_bm3_mgd",
    }:
        from scripts.forsterite_reproduce_katsura_2022 import FIT_INPUT

        with FIT_INPUT.open(newline="", encoding="utf-8") as stream:
            rows = list(csv.DictReader(stream))
        return direct(
            column(rows, "pressure_gpa"),
            column(rows, "volume_ratio") * record["eos"]["parameters"]["V0"],
            column(rows, "temperature_k"),
        )

    from scripts.reproduce_iron_source_papers import belongs, observation_set

    if belongs(rid) and "brown" not in rid:
        rows, pressure_column, reason = observation_set(record)
        if reason:
            raise ValueError(reason)
        return direct(
            [r[pressure_column] for r in rows],
            [r["volume_a3_conventional_cell"] for r in rows],
            [r["temperature_k"] for r in rows] if record.get("thermal") else None,
        )
    if "_nisr_2017_" in rid:
        from scripts.reproduce_nisr_2017_silica import PUBLISHED, load_rows

        rows = [r for r in load_rows() if r["table"] == PUBLISHED[rid][0]]
        return direct(
            [r["pressure_gpa"] for r in rows],
            [r["volume_a3"] for r in rows],
            volume_sigma=np.array([float(r["volume_sigma_a3"]) for r in rows]),
        )
    if rid.startswith("superhydrous_phase_b_lt_litasov_2007_"):
        from scripts.reproduce_litasov_2007_superhydrous_phase_b import (
            observations as litasov_rows,
        )

        _, (volume, temperature, pressure) = litasov_rows(
            rid.removeprefix("superhydrous_phase_b_lt_litasov_2007_")
        )
        return direct(pressure, volume, temperature if record.get("thermal") else None)
    if rid == "gold_fei_2007_vinet_2":
        from scripts.reproduce_fei_2007_gold import observations as fei_gold_rows

        _, data = fei_gold_rows()
        cold = table_rows("gold-dewaele-2004-table1-compression.csv")
        return direct(
            np.r_[column(cold, "ruby_pressure_revised_gpa"), data["p"]],
            np.r_[4 * column(cold, "atomic_volume_a3"), data["v"]],
            np.r_[np.full(len(cold), 300.0), data["t"]],
        )
    if rid == "gold_sakai_2025_rydberg_stacey_1":
        report = json.loads(
            (ROOT / "docs/data/sakai-2025-gold-reconstruction.json").read_text(
                encoding="utf-8"
            )
        )
        rows = [r for r in report["pressure_reductions"] if r["included_in_diagnostic"]]
        return direct(
            [r["calculated_pressure_gpa"] for r in rows], [r["volume_a3"] for r in rows]
        )
    if rid == "rhenium_sakai_2018_yokoo_pt_vinet":
        from scripts.reproduce_sakai_2018_rhenium import load_points

        rows = [r for r in load_points() if r["series"] == "rp01_yokoo_pt"]
        return direct([r["pressure_gpa"] for r in rows], [r["volume_a3"] for r in rows])
    if rid in SHEN_SELECTIONS:
        return _shen_series(record, datasets[identifiers[0]])
    if rid.startswith("nacl_b2_sakai_2011_"):
        from scripts.reproduce_sakai_2011_nacl_b2 import observations as sakai_rows

        scale = rid.removeprefix("nacl_b2_sakai_2011_").split("_", 1)[0]
        pressure, volume, sigma, _ = sakai_rows(scale)
        return Series(
            "+".join(identifiers),
            pressure,
            volume,
            None,
            None,
            sigma,
            None,
            "source/reconstructed Pt pressure",
            "volume_a3",
            None,
            "Sakai Table 1 plus selected Sata rows, as in the named audit",
        )
    if rid in {"nacl_b2_fei_2007", "nacl_b2_fei_2007_staged_thermal_refit"}:
        from scripts.reproduce_fei_2007_nacl_b2 import SOURCE, subset

        source = json.loads(SOURCE.read_text(encoding="utf-8"))
        cold, hot = subset(source, 3, 300), subset(source, 4, 1000)
        volume, temperature, pressure = [
            np.concatenate(values) for values in zip(cold, hot)
        ]
        return Series(
            "+".join(identifiers),
            pressure,
            volume,
            temperature,
            None,
            None,
            None,
            "pressure_gpa",
            "volume_a3",
            "temperature_k",
            "12 Figure 3 cold and 12 Figure 4 hot digitized markers",
        )
    if not identifiers or any(identifier not in datasets for identifier in identifiers):
        raise ValueError(
            "Selected audit datasets are not registered as replayable P–V rows"
        )
    if fit.get("fit_kind") == "shared_bm2_fe10_pressure_residuals":
        from scripts.validate_primary_eos_refits import _load_rows

        rows = [
            r
            for r in _load_rows(datasets[identifiers[0]])
            if str(r["fit_included"]) == "1"
        ]
        # These records explicitly inherit Fe10's relative compression, not its
        # absolute cell size. Keep that normalization visible in the plot caption.
        series = direct(
            column(rows, "pressure_gpa"),
            column(rows, "volume_a3") * record["eos"]["parameters"]["V0"] / 162.79,
        )
        return replace(
            series, volume_column="Fe10 V/V0 mapped to this record's source V0"
        )
    parts = [
        _series(document, record, datasets[identifier]) for identifier in identifiers
    ]
    if rid == "feo_b8_2_fischer_2011_bm3_1":
        from scripts.validate_primary_eos_refits import _load_rows

        ozawa = datasets["feo_ozawa_2010_table2_pvt"]
        parts = [
            _series(document, record, datasets["feo_fischer_2011_table_s1_pvt"]),
            _masked_series(
                _series(document, record, ozawa),
                np.array([r["phase"] == "B8" for r in _load_rows(ozawa)]),
                "all eight printed Ozawa B8 rows",
            ),
        ]
    if record.get("thermal") and len(parts) > 1:
        parts = [
            replace(
                part,
                temperature=np.full(part.pressure.shape, record["temperature_ref"]),
            )
            if part.temperature is None
            else part
            for part in parts
        ]
    if len(parts) == 1 and (
        rid.startswith("b4c_somayazulu_2023_") or rid == "mgo_dewaele_2000_bm3_mgd_5"
    ):
        part = parts[0]
        parts[0] = _masked_series(part, part.temperature > 300, "41 heated rows")
    if rid == "ice_vi_bezacier_2014_bm2_1":
        part = parts[0]
        parts[0] = _masked_series(
            part,
            (np.abs(part.temperature - 300) <= 1.5)
            | (np.abs(part.temperature - 340) <= 1),
            "source endpoint isotherms",
        )
    if rid == "e_feooh_gleason_2008_bm2_1":
        part = parts[0]
        correction = record["fit_protocol"]["temperature_handling"]
        parts[0] = replace(
            part,
            volume=part.volume
            * np.exp(
                -correction["alpha0_per_k"]
                * (part.temperature - correction["reference_temperature_k"])
            ),
            temperature=None,
            volume_column=part.volume_column + " reduced to 300 K",
        )
    if not record.get("thermal"):
        for index, part in enumerate(parts):
            if part.temperature is not None:
                delta = np.abs(part.temperature - record.get("temperature_ref", 300))
                parts[index] = replace(
                    _masked_series(
                        part,
                        delta <= max(5.0, float(np.min(delta)) + 1e-9),
                        part.selection + "; reference-T rows",
                    ),
                    temperature=None,
                    temperature_sigma=None,
                )
    if len(parts) == 1:
        return parts[0]

    def joined(field):
        values = [getattr(part, field) for part in parts]
        return (
            None if any(value is None for value in values) else np.concatenate(values)
        )

    return Series(
        "+".join(identifiers),
        joined("pressure"),
        joined("volume"),
        joined("temperature"),
        joined("pressure_sigma"),
        joined("volume_sigma"),
        joined("temperature_sigma"),
        parts[0].pressure_column,
        parts[0].volume_column,
        parts[0].temperature_column,
        "; ".join(part.selection for part in parts),
    )


def curve_comparison(
    document, record
) -> tuple[CurveComparison | DiagnosticComparison | None, str]:
    """Return comparable source/refit curves or the specific missing evidence."""
    document = convert_eosmat(document, version=3)
    record = next(
        r for r in document["eos_records"] if r["identifier"] == record["identifier"]
    )
    assessment = record["reproduction"]["assessment"]
    parameters, fit = assessment.get("parameters", []), assessment.get("fit", {})
    if not parameters:
        return (
            None,
            "No reconstructed coefficient set is recorded; a source/refit residual comparison is unavailable.",
        )
    try:
        diagnostic = diagnostic_comparison(document, record)
        if diagnostic is not None:
            diagnostic.independent_refit = assessment["scope"] == "independent_refit"
            return diagnostic, ""
        source = executable(document, record, parameters, "published")
        refit = executable(document, record, parameters, "refit")
        series = observations(document, record, fit)
        if len(series.pressure) != fit.get("observations"):
            raise ValueError(
                "Recovered row count differs from the recorded audit selection"
            )
        if not len(series.pressure):
            raise ValueError("No selected P–V observations are available")
        source_pressure = np.asarray(
            source.pressure(series.volume, series.temperature), dtype=float
        )
        refit_pressure = np.asarray(
            refit.pressure(series.volume, series.temperature), dtype=float
        )
        if not np.all(np.isfinite([source_pressure, refit_pressure])):
            raise ValueError(
                "The stored models are not finite on these observation coordinates"
            )
        source_rmse = float(np.sqrt(np.mean((source_pressure - series.pressure) ** 2)))
        refit_rmse = float(np.sqrt(np.mean((refit_pressure - series.pressure) ** 2)))
        # EIV objectives may report residuals at latent coordinates. In that case
        # source RMS still checks the row replay; the plotted refit RMS is explicitly
        # at observed coordinates and is not relabelled as the optimizer statistic.
        source_check = source_rmse
        note = ""
        if record["identifier"] == "ca_perovskite_noguchi_2013_bm2_mgd_1":
            from scripts.reproduce_noguchi_2013_casio3 import (
                SOURCE_TABLE_PATH,
                THERMAL_SOURCE_ORDERS,
                load_source_table,
            )

            rows = [
                r for r in load_source_table(SOURCE_TABLE_PATH) if r["fit_included"]
            ]
            hot = np.array([r["source_order"] in THERMAL_SOURCE_ORDERS for r in rows])
            source_check = float(
                np.sqrt(np.mean((source_pressure[hot] - series.pressure[hot]) ** 2))
            )
            note = "This plot includes all nine accepted 700 K static-stage rows and all 42 thermal-stage rows. Its RMS summarizes all 51 rows; the separately recorded thermal-stage RMS covers 42 rows only."
        if record["identifier"] == "gold_fei_2007_vinet_2":
            hot = np.arange(len(series.pressure)) >= 37
            source_check = float(
                np.sqrt(np.mean((source_pressure[hot] - series.pressure[hot]) ** 2))
            )
        if record["identifier"] == "fes_vi_morard_2026_bm3_mgd":
            with (
                ROOT
                / "peritheos/data/datasets/fes-morard-2026-author-without-additional-cold.csv"
            ).open(newline="", encoding="utf-8") as stream:
                rows = list(csv.DictReader(stream))
            hot = np.array([r["source_group"] == "morard_thermal" for r in rows])
            source_check = float(
                np.sqrt(np.mean((source_pressure[hot] - series.pressure[hot]) ** 2))
            )
            note = "All 21 literature cold rows and 146 thermal rows are shown. The plotted RMS covers all 167 inputs; the audit's separately recorded pressure RMS describes only the 146 thermal rows."
        for key, actual in (("published_rmse_gpa", source_check),):
            expected = fit.get(key)
            if expected is not None and not np.isclose(
                actual, expected, rtol=2e-4, atol=2e-6
            ):
                raise ValueError(
                    "Recovered source-curve RMS differs from the named audit; its row/calibration replay needs a dedicated adapter"
                )
        return CurveComparison(
            source,
            refit,
            series,
            source_pressure,
            refit_pressure,
            source_rmse,
            refit_rmse,
            source.volume_unit,
            assessment["scope"] == "independent_refit",
            note,
        ), ""
    except (
        ValueError,
        KeyError,
        TypeError,
        AttributeError,
        NotImplementedError,
    ) as error:
        return None, str(error)
