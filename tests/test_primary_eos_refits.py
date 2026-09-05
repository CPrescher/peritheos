import json
from collections import Counter
from pathlib import Path

import pytest

from peritheos import list_eos_record_documents

ROOT = Path(__file__).resolve().parents[1]
LEDGER_PATH = ROOT / "docs" / "data" / "primary-eos-refits.json"
MARKDOWN_PATH = ROOT / "docs" / "primary-eos-refits.md"


def load_ledger():
    return json.loads(LEDGER_PATH.read_text(encoding="utf-8"))


def test_primary_refit_ledger_covers_every_bundled_record_once():
    ledger = load_ledger()
    identifiers = [item["record_identifier"] for item in ledger["records"]]

    assert ledger["format"] == "peritheos.primary-eos-refit-validation"
    assert ledger["format_version"] == 1
    assert len(identifiers) == len(set(identifiers)) == 253
    assert set(identifiers) == set(list_eos_record_documents())


def test_primary_refit_summary_and_results_are_internally_consistent():
    ledger = load_ledger()
    statuses = Counter(item["status"] for item in ledger["records"])

    assert ledger["summary"] == {"total": 253, **dict(sorted(statuses.items()))}
    assert statuses == {
        "parity": 131,
        "similar": 52,
        "parity_not_achieved": 18,
        "not_refittable": 52,
    }
    assert all(
        item.get("reason")
        for item in ledger["records"]
        if item["status"] == "not_refittable"
    )


def test_primary_refit_regression_examples_and_documentation_coverage():
    ledger = load_ledger()
    by_identifier = {item["record_identifier"]: item for item in ledger["records"]}
    markdown = MARKDOWN_PATH.read_text(encoding="utf-8")

    luo = by_identifier["mgo_b1_luo_2023_vinet_thermal_5"]
    assert luo["status"] == "not_refittable"
    assert luo["dataset_identifiers"] == ["mgo_luo_2023_table1_shock"]
    assert "derived EOS output" in luo["reason"]
    reynard_ruby = by_identifier["akimotoite_reynard_1996_bm3_ruby_2"]
    reynard_ice = by_identifier["akimotoite_reynard_1996_bm3_ice_vii_3"]
    assert reynard_ruby["status"] == reynard_ice["status"] == "parity"
    assert reynard_ruby["observations"] == 16
    assert reynard_ice["observations"] == 14
    assert reynard_ruby["columns"]["pressure"] == "ruby_pressure_gpa"
    assert reynard_ice["columns"]["pressure"] == "ice_vii_pressure_gpa"
    assert reynard_ruby["fixed_parameters"] == reynard_ice["fixed_parameters"] == ["K0"]
    assert [item["refit"] for item in reynard_ruby["parameters"]] == pytest.approx(
        [262.421420404, 6.915511476]
    )
    assert [item["refit"] for item in reynard_ice["parameters"]] == pytest.approx(
        [262.580091892, 5.055863026]
    )
    assert (
        "same 16-row compression/decompression experiment"
        in reynard_ruby["source_notes"]
    )
    assert "distinct fit" in reynard_ice["source_notes"]

    assert by_identifier["aragonite_martinez_1996_bm2_2"]["status"] == "parity"
    assert by_identifier["kcl_campbell_1991_bm2_1"]["status"] == "parity"
    phase_egg = by_identifier["phase_egg_mookherjee_2019_bm3_lp_1"]
    assert phase_egg["status"] == "parity"
    assert phase_egg["observations"] == 11
    assert phase_egg["published_rmse_gpa"] == pytest.approx(0.03011320577)
    assert phase_egg["rmse_gpa"] == pytest.approx(0.02339254240)
    assert [item["refit"] for item in phase_egg["parameters"]] == pytest.approx(
        [210.2176007194, 164.8339602796, 7.0813721266]
    )
    assert (
        "cannot reconstruct the source energy-fit protocol"
        in phase_egg["qualification"]
    )
    walker = by_identifier["kcl_walker_2002_bm3_2"]
    assert walker["status"] == "similar"
    assert walker["free_parameters"] == [
        "rt_eos.K0",
        "rt_eos.K0_prime",
        "alpha_KT",
    ]
    assert [item["refit"] for item in walker["parameters"]] == pytest.approx(
        [23.7753978405, 4.4158737620, 0.002766245688]
    )
    assert walker["stages"][0]["observations"] == 8
    assert walker["stages"][1]["observations"] == 39
    assert "not a simultaneous four-parameter fit" in walker["qualification"]
    simultaneous = walker["simultaneous_fit_diagnostic"]
    simultaneous_parameters = simultaneous["peritheos_four_parameter_refit"]
    assert simultaneous_parameters["rt_eos.V0"] == pytest.approx(55.39232501)
    assert simultaneous_parameters["rt_eos.K0"] == pytest.approx(14.188725516)
    assert simultaneous_parameters["rt_eos.K0_prime"] == pytest.approx(7.156559906)
    assert simultaneous_parameters["alpha_KT"] == pytest.approx(0.002696920764)
    tateno = by_identifier["kcl_b2_tateno_2019_vinet_4"]
    assert tateno["status"] == "parity"
    assert tateno["dataset_identifiers"] == ["kcl_tateno_2019_table_s1_pvt"]
    assert [item["refit"] for item in tateno["parameters"]] == pytest.approx(
        [18.344616371, 5.600955151, 2.295185193, 0.824899722]
    )
    assert "Corrected final-publication reproduction" in tateno["qualification"]
    chidester = by_identifier["kcl_b2_chidester_2021_bm3_5"]
    assert chidester["status"] == "parity"
    assert chidester["observations"] == 278
    assert chidester["dataset_identifiers"] == [
        "kcl_dewaele_2012_table1_compression",
        "kcl_chidester_2021_supplemental_pvt",
    ]
    assert [item["refit"] for item in chidester["parameters"]] == pytest.approx(
        [53.203552571, 23.972141965, 4.557980658, 2.917137053, 0.965242582]
    )
    assert chidester["high_temperature_refit_rmse_gpa"] == pytest.approx(1.581924568)
    assert chidester["source_reported_pressure_rmse_gpa"] == 1.6
    assert "155 new high-temperature rows" in chidester["qualification"]
    assert by_identifier["b4c_somayazulu_2023_bm3_1"]["status"] == (
        "parity_not_achieved"
    )
    assert by_identifier["b4c_somayazulu_2023_berman_refit"]["status"] == ("parity")
    assert by_identifier["gold_shen_2026_vinet_3"]["status"] == "not_refittable"
    assert by_identifier["ca_perovskite_sun_2016_bm3_3"]["status"] == ("not_refittable")
    assert by_identifier["ca_perovskite_tetragonal_sun_2022_bm3_1"]["status"] == (
        "not_refittable"
    )
    fu_bm2 = by_identifier["mg088fe010al014si090o3_bridgmanite_fu_2024_bm2_1"]
    fu_bm3 = by_identifier["mg088fe010al014si090o3_bridgmanite_fu_2024_bm3_2"]
    assert fu_bm2["status"] == "similar"
    assert fu_bm3["status"] == "parity"
    nio_hugoniot = by_identifier["nickel_oxide_noguchi_1999_linear_hugoniot_2"]
    assert nio_hugoniot["status"] == "parity"
    assert nio_hugoniot["observations"] == 8
    assert nio_hugoniot["rmse_shock_velocity_km_s"] == pytest.approx(0.0556367146096)
    assert by_identifier["mgo_b1_duffy_ahrens_1995_hugoniot_5"]["status"] == (
        "not_refittable"
    )
    neon_bm3 = by_identifier["neon_fcc_fei_2007_bm3_1"]
    neon_vinet = by_identifier["neon_fcc_fei_2007_vinet_2"]
    assert neon_bm3["status"] == "parity"
    assert neon_vinet["status"] == "parity"
    assert neon_bm3["observations"] == 34
    assert neon_vinet["observations"] == 34
    assert neon_bm3["dataset_identifiers"] == [
        "neon_fei_2007_figure5_digitized",
        "neon_hemley_1989_table1_fei_recalculated",
    ]
    assert neon_bm3["observed_pressure_range_gpa"] == pytest.approx(
        [10.047704716, 115.715945233]
    )
    assert [item["refit"] for item in neon_bm3["parameters"]] == pytest.approx(
        [1.4775172078, 7.8503378205]
    )
    assert [item["refit"] for item in neon_vinet["parameters"]] == pytest.approx(
        [1.1439324121, 8.2580743609]
    )
    assert neon_bm3["free_parameters"] == ["K0", "K0_prime"]
    assert neon_vinet["free_parameters"] == ["rt_eos.K0", "rt_eos.K0_prime"]
    assert "Conditional partial reproduction" in neon_bm3["qualification"]
    assert "Conditional partial reference-isotherm" in neon_vinet["qualification"]
    assert (
        "Finger's low-pressure rows remain unavailable" in (neon_bm3["qualification"])
    )
    hemley = by_identifier["neon_fcc_hemley_1989_bm3_refit"]
    assert hemley["status"] == "parity"
    assert hemley["observations"] == 21
    assert hemley["free_parameters"] == ["K0", "K0_prime"]
    stishovite = by_identifier["sio2_stv_andr_wang_2012_vinet_mgd_2"]
    assert stishovite["status"] == "parity"
    assert stishovite["observations"] == 56
    assert stishovite["objective"] == "errors_in_variables"
    assert stishovite["fixed_parameters"] == ["V0", "Tr", "a", "n"]
    assert [item["refit"] for item in stishovite["parameters"]] == pytest.approx(
        [294.6773736503, 4.8789956491, 1134.709214856, 1.6555196866, 3.0469258011]
    )
    assert all(item["within_combined_2sigma"] for item in stishovite["parameters"])
    mo2c = by_identifier["molybenum_carbide_mo2c_haines_2001_bm3_1"]
    assert mo2c["status"] == "similar"
    assert mo2c["observations"] == 16
    assert mo2c["fixed_parameters"] == ["V0"]
    assert mo2c["free_parameters"] == ["K0", "K0_prime"]
    assert [item["refit"] for item in mo2c["parameters"]] == pytest.approx(
        [325.874185450, 4.909198583]
    )
    assert all(item["within_combined_2sigma"] for item in mo2c["parameters"])
    assert "Corrected source-scope reproduction" in mo2c["qualification"]
    mo2c_refit = by_identifier["molybenum_carbide_mo2c_haines_2001_bm3_refit"]
    assert mo2c_refit["status"] == "parity"
    assert mo2c_refit["observations"] == 16
    assert mo2c_refit["fixed_parameters"] == ["V0"]
    assert mo2c_refit["free_parameters"] == ["K0", "K0_prime"]
    assert [item["refit"] for item in mo2c_refit["parameters"]] == pytest.approx(
        [325.8744132414, 4.9091866068]
    )
    assert "Explicit Peritheos refit record" in mo2c_refit["qualification"]
    palladium = by_identifier["palladium_baty_2024_bm3_1"]
    assert palladium["status"] == "parity_not_achieved"
    assert palladium["observations"] == 78
    assert "dedicated palladium reproduction" in palladium["qualification"]
    palladium_diagnostic = palladium["fit_protocol_diagnostic"]
    assert (
        palladium_diagnostic["source_table_rows_verified_against_official_latex"] == 78
    )
    assert palladium_diagnostic["source_equation"] == (
        "standard third-order Birch-Murnaghan, Equation (1)"
    )
    assert palladium_diagnostic["unweighted_volume_all_rows"][
        "parameters"
    ] == pytest.approx(
        {"V0": 59.470126905, "K0": 178.823784858, "K0_prime": 5.009936202}
    )
    assert palladium_diagnostic["fixed_published_V0_all_rows"][
        "parameters"
    ] == pytest.approx({"V0": 58.88, "K0": 207.476522488, "K0_prime": 4.296362911})
    assert palladium_diagnostic["fixed_published_V0_pressure_at_least_40_gpa"][
        "parameters"
    ] == pytest.approx({"V0": 58.88, "K0": 193.664031670, "K0_prime": 5.024310596})
    assert palladium_diagnostic["upper_bound_errors_in_variables_all_rows"][
        "reduced_chi_square"
    ] == pytest.approx(8.107727533)
    frost_cross_check = palladium_diagnostic["frost_2023_vinet_cross_check"]
    assert frost_cross_check["parameters"] == pytest.approx(
        {"V0": 58.678, "K0": 189.3, "K0_prime": 5.473}
    )
    assert frost_cross_check["pressure_residual_rmse_gpa"] == pytest.approx(2.039354740)
    assert frost_cross_check["observation_minus_curve_volume_atomic_a3"] == (
        pytest.approx(
            {
                "rows_positive": 77,
                "mean": 0.076869666,
                "rmse": 0.093553190,
                "maximum": 0.200033281,
            }
        )
    )
    assert palladium_diagnostic["table_monotonicity"] == pytest.approx(
        {
            "adjacent_volume_increases_after_sorting_by_pressure": 12,
            "largest_increase_a3_conventional_cell": 0.1678738,
        }
    )
    fedotenko_bm3 = by_identifier["palladium_fedotenko_2020_bm3_1"]
    assert fedotenko_bm3["status"] == "parity"
    assert fedotenko_bm3["observations"] == 15
    assert [item["refit"] for item in fedotenko_bm3["parameters"]] == (
        pytest.approx([58.853666470, 158.668115182, 9.783509441])
    )
    fedotenko_bm2 = by_identifier["palladium_fedotenko_2020_bm2_2"]
    assert fedotenko_bm2["status"] == "parity"
    assert [item["refit"] for item in fedotenko_bm2["parameters"]] == (
        pytest.approx([202.667325846])
    )
    assert by_identifier["palladium_baty_2024_bm3_dft_2"]["status"] == (
        "not_refittable"
    )
    assert by_identifier["palladium_guigue_2020_vinet_1"]["status"] == (
        "not_refittable"
    )
    frost_vinet = by_identifier["palladium_frost_2023_vinet_1"]
    assert frost_vinet["status"] == "similar"
    assert frost_vinet["observations"] == 93
    assert frost_vinet["dataset_identifiers"] == [
        "palladium_frost_2023_tables_i_iii_compression"
    ]
    assert [item["refit"] for item in frost_vinet["parameters"]] == pytest.approx(
        [58.766517848, 184.190604957, 5.603711668]
    )
    frost_bm3 = by_identifier["palladium_frost_2023_bm3_2"]
    assert frost_bm3["status"] == "similar"
    assert frost_bm3["observations"] == 93
    assert frost_bm3["dataset_identifiers"] == [
        "palladium_frost_2023_tables_i_iii_compression"
    ]
    assert [item["refit"] for item in frost_bm3["parameters"]] == pytest.approx(
        [58.684580564, 191.476720898, 5.147381676]
    )
    rbcl = by_identifier["rbcl_b2_campbell_1994_bm3_1"]
    assert rbcl["status"] == "parity"
    assert rbcl["observations"] == 24
    assert rbcl["dataset_identifiers"] == ["rbcl_campbell_1994_table1_compression"]
    assert [item["refit"] for item in rbcl["parameters"]] == pytest.approx(
        [17.8808436600, 5.2381519592]
    )
    cscl = by_identifier["cscl_campbell_1994_bm3_1"]
    assert cscl["status"] == "parity"
    assert cscl["observations"] == 22
    assert cscl["dataset_identifiers"] == [
        "cscl_campbell_1994_table1_compression",
        "cscl_yagi_1978_table1_compression",
    ]
    assert [item["refit"] for item in cscl["parameters"]] == pytest.approx(
        [17.6967334851, 5.2026008383]
    )
    assert "Complete source-data reproduction" in cscl["qualification"]
    dewaele_mgo = by_identifier["mgo_dewaele_2000_bm3_mgd_5"]
    assert dewaele_mgo["status"] == "parity"
    assert dewaele_mgo["observations"] == 41
    assert dewaele_mgo["selection"] == "41 heated Table 2 rows"
    assert dewaele_mgo["free_parameters"] == ["q"]
    assert dewaele_mgo["fixed_parameters"] == [
        "V0",
        "K0",
        "K0_prime",
        "Tr",
        "theta0",
        "gamma0",
        "n",
    ]
    assert dewaele_mgo["parameters"][0]["refit"] == pytest.approx(0.847089125)
    assert (
        "Conditional current-study thermal reproduction"
        in (dewaele_mgo["qualification"])
    )
    phase_egg = by_identifier["phase_egg_schulze_2018_bm3_1"]
    assert phase_egg["status"] == "parity"
    assert phase_egg["observations"] == 15
    assert phase_egg["selection"] == "used_in_published_fit=1"
    assert phase_egg["dataset_identifiers"] == [
        "phase_egg_schulze_2018_table1_compression"
    ]
    assert [item["refit"] for item in phase_egg["parameters"]] == pytest.approx(
        [214.0810595, 152.8459201, 8.597514758]
    )
    assert phase_egg["reduced_chi_square"] == pytest.approx(5.486063469)
    rh2o3 = by_identifier["alumina_rh2o3_ii_shi_2022_bm3_mgd_1"]
    assert rh2o3["status"] == "similar"
    assert rh2o3["observations"] == 75
    assert rh2o3["dataset_identifiers"] == ["alumina_rh2o3_ii_shi_2022_table_s2_pvt"]
    assert rh2o3["fixed_parameters"] == ["K0_prime", "Tr", "q", "n"]
    assert rh2o3["free_parameters"] == [
        "rt_eos.V0",
        "rt_eos.K0",
        "theta0",
        "gamma0",
    ]
    assert [item["refit"] for item in rh2o3["parameters"]] == pytest.approx(
        [167.1940050282, 239.4153140110, 766.2583682877, 1.5502121628]
    )
    assert rh2o3["published_rmse_gpa"] == pytest.approx(1.1433652954)
    assert rh2o3["rmse_gpa"] == pytest.approx(0.8659771871)
    assert all(item["within_combined_2sigma"] for item in rh2o3["parameters"])
    coo = by_identifier["coo_clendenen_1966_murnaghan_1"]
    assert coo["status"] == "parity_not_achieved"
    tradeoff = coo["coefficient_tradeoff_diagnostic"]
    assert tradeoff["K0_K0_prime_correlation"] == pytest.approx(-0.971026845)
    assert tradeoff["conditional_fits"]["K0_fixed"]["refit_value"] == (
        pytest.approx(3.839109480)
    )
    assert tradeoff["conditional_fits"]["K0_prime_fixed"]["refit_value"] == (
        pytest.approx(189.0898192)
    )
    assert tradeoff["rounding_intervals"]["rows_outside_published_curve"] == 4
    assert tradeoff["rounding_intervals"]["minimum_pressure_gaps_kbar"] == (
        pytest.approx([1.975447794, 4.436419298, 0.633380963, 3.635154586])
    )
    assert tradeoff["published_pair_delta_chi_square"] == pytest.approx(3.561705438)
    assert "properly balanced" in coo["source_notes"]
    goethite = by_identifier["goethite_gleason_2008_bm3_1"]
    assert goethite["status"] == "parity_not_achieved"
    assert goethite["fit_kind"] == "joint_pvt"
    assert goethite["observations"] == 65
    assert goethite["free_parameters"] == [
        "rt_eos.K0",
        "rt_eos.K0_prime",
    ]
    assert [item["refit"] for item in goethite["parameters"]] == pytest.approx(
        [183.337839186, 0.0]
    )
    anomaly = goethite["source_table_anomaly_diagnostic"]
    assert anomaly["depository_row"] == 32
    assert anomaly["published_model_volume_a3"] == pytest.approx(133.105228337)
    assert anomaly["published_model_pressure_residual_gpa"] == pytest.approx(
        16.293430720
    )
    assert anomaly["published_curve_rmse_gpa"] == pytest.approx(
        {"all_65_rows": 2.540000975, "excluding_row_32": 1.550596848}
    )
    clean_goethite = anomaly["refit_sensitivity"][
        "errors_in_variables_excluding_row_32"
    ]
    assert clean_goethite["rt_eos.K0"] == pytest.approx(177.642304461)
    assert clean_goethite["rt_eos.K0_prime"] == pytest.approx(1.541879108)
    assert "not recommended for quantitative" in goethite["source_notes"]
    ice_vi = by_identifier["ice_vi_bezacier_2014_bm2_1"]
    assert ice_vi["status"] == "parity"
    assert ice_vi["observations"] == 30
    assert ice_vi["selection"] == (
        "23 rows at 298.7-300.7 K and seven rows at 340.0-340.7 K"
    )
    assert [item["refit"] for item in ice_vi["parameters"]] == pytest.approx(
        [235.353768434, 14.015825054, 0.000148152628]
    )
    assert "Figure 2 displays" in ice_vi["qualification"]
    all_ice_vi_rows = ice_vi["all_table_rows_diagnostic"]
    assert all_ice_vi_rows["observations"] == 45
    assert all_ice_vi_rows["parameters"]["alpha0"] == pytest.approx(0.000052865775)
    assert all_ice_vi_rows["parameters"]["rt_eos.K0"] == pytest.approx(15.562571244)
    explained = [
        item
        for item in ledger["records"]
        if item["status"] in {"similar", "parity_not_achieved", "refit_failed"}
    ]
    assert markdown.count("### `") == len(explained) == 70
    assert all(identifier in markdown for identifier in by_identifier)
    failed = [
        item
        for item in ledger["records"]
        if item["status"] in {"parity_not_achieved", "refit_failed"}
    ]
    for item in failed:
        anchor = f"investigation-{item['record_identifier']}"
        assert markdown.count(f'<a id="{anchor}"></a>') == 1
        assert markdown.count(f"](#{anchor})") >= 2
