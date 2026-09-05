import csv
import hashlib
from importlib import resources

import numpy as np
import pytest
from scipy.optimize import least_squares

from peritheos import Material, get_material_document
from peritheos.eos.rt import BM3
from peritheos.eos.thermal import ThermalReferenceStateEOS

RECORD_ID = "bridgmanite_funamori_1996_bm3_4"
WANG_RECORD_ID = "bridgmanite_wang_1994_bm3_5"
DATASET_EXPECTATIONS = {
    "bridgmanite_wang_1994_table3_pvt": (
        "48af346ab6224affa1e207d49ee86166475982944da27e815abf7ba2a2a08ff1",
        79,
        79,
    ),
    "bridgmanite_utsumi_1995_table1_pvt": (
        "fb879dd8240eb3d152bcfdc03b64d345f9ca8d948b2bfb125202f9b4863eed3d",
        71,
        34,
    ),
    "bridgmanite_funamori_1996_table4_pvt": (
        "c2798ae8a4250a16548db90bbbd7fe2786a5d7ea212643a096d95459a2e437d0",
        83,
        15,
    ),
    "bridgmanite_funamori_1996_table5_25gpa": (
        "9f96da5317872f2aa54515f9d5aecd93cfdf7e3a602afc55f891f768988812d9",
        18,
        None,
    ),
}


def _record_and_datasets():
    document = get_material_document("bridgmanite")
    record = next(
        item for item in document["eos_records"] if item["identifier"] == RECORD_ID
    )
    datasets = {
        item["identifier"]: item
        for item in document["datasets"]
        if item["identifier"] in DATASET_EXPECTATIONS
    }
    return document, record, datasets


def _rows(dataset):
    path = resources.files("peritheos.data").joinpath(dataset["resource"]["path"])
    payload = path.read_bytes()
    return payload, list(csv.DictReader(payload.decode("utf-8").splitlines()))


def test_funamori_1996_record_preserves_the_published_constrained_model():
    document, record, datasets = _record_and_datasets()
    executable = Material.from_eosmat(
        document, record_identifiers=[RECORD_ID]
    ).get_eos_record(RECORD_ID)

    assert record["reference"]["doi"].lower() == "10.1029/95jb03732"
    assert record["eos"]["parameters"] == {
        "V0": 162.32,
        "K0": 261.0,
        "K0_prime": 4.0,
    }
    thermal = record["thermal"]
    assert thermal["thermal_expansion_law"] == (
        "linear_temperature_inverse_square"
    )
    assert thermal["reference_volume_law"] == "integrated_expansivity"
    assert thermal["parameters"] == {
        "Tr": 300.0,
        "alpha0": 1.982e-5,
        "alpha1": 0.818e-8,
        "alpha_inverse_square": 0.474,
        "dK_dT": -0.0280,
    }
    alpha_300 = (
        thermal["parameters"]["alpha0"]
        + 300.0 * thermal["parameters"]["alpha1"]
        - thermal["parameters"]["alpha_inverse_square"] / 300.0**2
    )
    # a0 is printed to only four significant digits in Table 6.
    assert alpha_300 == pytest.approx(1.7e-5, abs=1.0e-8)
    assert set(record["fit_datasets"]) == set(DATASET_EXPECTATIONS) - {
        "bridgmanite_funamori_1996_table5_25gpa"
    }
    assert set(datasets) == set(DATASET_EXPECTATIONS)
    assert executable.eos.thermal_expansion_law == (
        "linear_temperature_inverse_square"
    )
    assert record["pressure_calibration"]["status"] == "partially_resolved"
    assert record["pressure_calibration"]["recalculation"]["status"] == (
        "reference_eos_not_bundled"
    )


def test_funamori_1996_primary_tables_are_complete_selected_and_checksummed():
    _, _, datasets = _record_and_datasets()
    for identifier, (checksum, row_count, fit_count) in DATASET_EXPECTATIONS.items():
        dataset = datasets[identifier]
        payload, rows = _rows(dataset)
        assert hashlib.sha256(payload).hexdigest() == checksum
        assert len(rows) == row_count
        assert list(rows[0]) == [column["name"] for column in dataset["columns"]]
        expected_records = [RECORD_ID]
        if identifier == "bridgmanite_wang_1994_table3_pvt":
            expected_records.append(WANG_RECORD_ID)
        assert dataset["used_by_eos_records"] == expected_records
        if fit_count is not None:
            assert sum(row["fit_included"] == "1" for row in rows) == fit_count

    _, utsumi_rows = _rows(datasets["bridgmanite_utsumi_1995_table1_pvt"])
    assert all(
        row["fit_included"] == "0"
        for row in utsumi_rows
        if row["cycle"].startswith("first_") and float(row["temperature_k"]) < 650
    )
    _, funamori_rows = _rows(datasets["bridgmanite_funamori_1996_table4_pvt"])
    selected = [row for row in funamori_rows if row["fit_included"] == "1"]
    assert all(row["table"] == "4a" and row["run"] == "1" for row in selected)
    assert all(row["pressure_nd1_gpa"] for row in selected)
    assert [row["source_row"] for row in selected[:2]] == ["17", "18"]
    assert selected[-1]["source_row"] == "32"


def test_funamori_1996_independent_equal_pressure_residual_refit():
    _, record, datasets = _record_and_datasets()
    observations = []
    for identifier, pressure_column in (
        ("bridgmanite_wang_1994_table3_pvt", "pressure_gpa"),
        ("bridgmanite_utsumi_1995_table1_pvt", "pressure_gpa"),
        ("bridgmanite_funamori_1996_table4_pvt", "pressure_nd1_gpa"),
    ):
        _, rows = _rows(datasets[identifier])
        observations.extend(
            (
                float(row[pressure_column]),
                float(row["temperature_k"]),
                float(row["unit_cell_volume_a3"])
                / float(row["normalization_volume_a3"]),
            )
            for row in rows
            if row["fit_included"] == "1"
        )

    pressure = np.array([row[0] for row in observations])
    temperature = np.array([row[1] for row in observations])
    volume_ratio = np.array([row[2] for row in observations])

    def residual(parameters):
        alpha1, alpha_inverse_square, dK_dT = parameters
        alpha0 = (
            1.7e-5
            - 300.0 * alpha1
            + alpha_inverse_square / 300.0**2
        )
        eos = ThermalReferenceStateEOS(
            BM3(1.0, 261.0, 4.0),
            300.0,
            alpha0,
            dK_dT,
            alpha1,
            thermal_expansion_law="linear_temperature_inverse_square",
            alpha_inverse_square=alpha_inverse_square,
        )
        return np.asarray(eos.pressure(volume_ratio, temperature)) - pressure

    fit = least_squares(
        residual,
        [0.818e-8, 0.474, -0.0280],
        x_scale="jac",
        xtol=1.0e-13,
        ftol=1.0e-13,
        gtol=1.0e-13,
    )
    published = np.array([0.818e-8, 0.474, -0.0280])
    published_sigma = np.array([0.257e-8, 0.118, 0.0032])
    diagnostics = record["scientific_validation"]["numerical_reproduction"]

    assert fit.success
    assert len(observations) == 128
    assert fit.x == pytest.approx(
        [8.53841399e-9, 0.465685536, -0.0283729204], abs=3.0e-7
    )
    assert np.all(np.abs(fit.x - published) < published_sigma)
    residual_std = np.std(fit.fun, ddof=3)
    assert residual_std == pytest.approx(
        diagnostics["independent_refit_residual_standard_deviation_gpa"],
        abs=2.0e-9,
    )
    assert abs(residual_std - diagnostics["published_residual_standard_deviation_gpa"]) < 0.02


def test_funamori_1996_combined_eos_matches_table5_within_cumulative_errors():
    document, _, datasets = _record_and_datasets()
    executable = Material.from_eosmat(
        document, record_identifiers=[RECORD_ID]
    ).get_eos_record(RECORD_ID)
    _, rows = _rows(datasets["bridgmanite_funamori_1996_table5_25gpa"])

    standardized_residuals = []
    for row in rows:
        calculated = executable.volume(25.0, float(row["temperature_k"])) / 162.32
        observed = float(row["volume_ratio"])
        error = float(row["volume_ratio_cumulative_error"])
        standardized_residuals.append(abs(calculated - observed) / error)

    # Table 5 is an independently derived AA1-pressure-scale interpolation,
    # not an input to the combined ND1 fit. The source compares them as a
    # consistency check, so agreement is tested against its cumulative errors.
    assert max(standardized_residuals) < 0.33


def test_wang_1994_record_preserves_the_published_bm3_thermal_eos():
    document = get_material_document("bridgmanite")
    record = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == WANG_RECORD_ID
    )
    executable = Material.from_eosmat(
        document, record_identifiers=[WANG_RECORD_ID]
    ).get_eos_record(WANG_RECORD_ID)

    assert record["reference"]["doi"] == "10.1016/0031-9201(94)90109-0"
    assert record["eos"]["parameters"] == {
        "V0": 162.2953463,
        "K0": 261.0,
        "K0_prime": 4.0,
    }
    assert record["fixed_parameters"] == ["V0", "K0", "K0_prime"]
    assert record["thermal"]["parameters"] == {
        "Tr": 300.0,
        "alpha0": 1.647e-5,
        "alpha1": 0.821e-8,
        "dK_dT": -0.0209,
    }
    assert record["thermal"]["parameter_errors"] == {
        "Tr": None,
        "alpha0": 0.180e-5,
        "alpha1": 0.337e-8,
        "dK_dT": 0.0102,
    }
    assert record["thermal"]["thermal_expansion_law"] == "linear_temperature"
    assert record["thermal"]["reference_volume_law"] == "integrated_expansivity"
    assert record["fit_datasets"] == ["bridgmanite_wang_1994_table3_pvt"]
    assert executable.eos.thermal_expansion_law == "linear_temperature"

    # Wang normalized each experimental block. The executable absolute anchor
    # is transparently derived from the most precise 295 K Table 3 measurement.
    integral_295_to_300 = 1.647e-5 * 5 + 0.5 * 0.821e-8 * (300**2 - 295**2)
    assert 162.280 * np.exp(integral_295_to_300) == pytest.approx(
        record["eos"]["parameters"]["V0"], abs=5.0e-8
    )


def test_wang_1994_independent_printed_volume_uncertainty_refit():
    document = get_material_document("bridgmanite")
    record = next(
        item
        for item in document["eos_records"]
        if item["identifier"] == WANG_RECORD_ID
    )
    dataset = next(
        item
        for item in document["datasets"]
        if item["identifier"] == "bridgmanite_wang_1994_table3_pvt"
    )
    _, rows = _rows(dataset)
    pressure = np.array([float(row["pressure_gpa"]) for row in rows])
    temperature = np.array([float(row["temperature_k"]) for row in rows])
    normalization = np.array(
        [float(row["normalization_volume_a3"]) for row in rows]
    )
    volume_ratio = np.array(
        [float(row["unit_cell_volume_a3"]) for row in rows]
    ) / normalization
    volume_ratio_uncertainty = np.array(
        [float(row["unit_cell_volume_printed_uncertainty_a3"]) for row in rows]
    ) / normalization

    def eos(parameters):
        return ThermalReferenceStateEOS(
            BM3(1.0, 261.0, 4.0),
            300.0,
            parameters[0],
            parameters[2],
            parameters[1],
            thermal_expansion_law="linear_temperature",
        )

    def residual(parameters):
        calculated = np.asarray(eos(parameters).volume(pressure, temperature))
        return (calculated - volume_ratio) / volume_ratio_uncertainty

    published = np.array([1.647e-5, 0.821e-8, -0.0209])
    published_sigma = np.array([0.180e-5, 0.337e-8, 0.0102])
    fit = least_squares(
        residual,
        published,
        x_scale="jac",
        xtol=1.0e-13,
        ftol=1.0e-13,
        gtol=1.0e-13,
        max_nfev=2000,
    )
    diagnostics = record["scientific_validation"]
    refit = diagnostics["independent_refit"]
    calculated_published = np.asarray(
        eos(published).volume(pressure, temperature)
    )
    published_volume_rmse = np.sqrt(
        np.mean((calculated_published - volume_ratio) ** 2)
    )

    assert fit.success
    assert len(rows) == 79
    np.testing.assert_allclose(
        fit.x,
        [
            refit["parameters"]["alpha0"],
            refit["parameters"]["alpha1"],
            refit["parameters"]["dK_dT"],
        ],
        rtol=5.0e-7,
        atol=1.0e-12,
    )
    assert np.max(np.abs((fit.x - published) / published_sigma)) < 0.44
    assert published_volume_rmse == pytest.approx(
        diagnostics["numerical_reproduction"][
            "published_curve_normalized_volume_rmse"
        ],
        abs=1.0e-9,
    )
    assert refit["objective_limitation"].startswith(
        "Wang et al. show propagated volume, pressure, and systematic error bars"
    )
