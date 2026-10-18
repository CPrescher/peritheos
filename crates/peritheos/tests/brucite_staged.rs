use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn brucite_catalog_preserves_only_the_three_retained_published_isotherms() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/brucite.eosmat");
    let material = load_eosmat(path).unwrap();
    assert_eq!(material.document["atom_sites"][2]["z"], 0.4108);
    assert_eq!(material.document["formula_units_per_cell"], 1);
    let static_record = material.record("brucite_fukui_2003_300k_bm3").unwrap();
    assert_eq!(
        static_record.document["validity"]["pressure_gpa"],
        serde_json::json!([0, 22.56])
    );
    assert_eq!(
        static_record.document["evidence"]["scientific_validation"]["status"],
        "primary_source_validated"
    );
    assert!(material.record("brucite_fukui_2003_873k_bm3").is_none());
    assert!(material
        .record("brucite_fukui_2003_300k_anchored_thermal")
        .is_none());
    assert!(material
        .record("brucite_peritheos_fukui_2003_thermal_refit")
        .is_none());
    for t in [300, 473, 673] {
        let record = material
            .record(&format!("brucite_fukui_2003_{t}k_bm3"))
            .unwrap();
        let volume = record.volume(10.0, f64::from(t)).unwrap();
        assert!((record.pressure(volume, f64::from(t)).unwrap() - 10.0).abs() < 1e-9);
    }
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}

#[test]
fn brucite_archived_diagnostics_preserve_equations_without_recommendation() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/brucite-diagnostic-archive.eosmat");
    let material = load_eosmat(path).unwrap();
    let thermal = material
        .record("brucite_fukui_2003_300k_anchored_thermal")
        .unwrap();
    assert_eq!(
        thermal.document["scientific_validation"]["status"],
        "not_reproduced"
    );
    assert_eq!(
        thermal.document["validity"]["temperature_k"],
        serde_json::json!([300, 873])
    );
    let t: f64 = 400.0;
    let v0 = 40.917 * (6.69e-5 * (t - 300.0)).exp();
    let eta = (v0 / 40.0).powf(1.0 / 3.0);
    let pressure = 1.5
        * (41.8 - 0.0319 * (t - 300.0))
        * (eta.powi(7) - eta.powi(5))
        * (1.0 + 0.75 * (6.6 + 0.00571 * (t - 300.0) - 4.0) * (eta.powi(2) - 1.0));
    assert!((thermal.pressure(40.0, t).unwrap() - pressure).abs() < 1e-12);
    let refit = material
        .record("brucite_peritheos_fukui_2003_873k_bm3_refit")
        .unwrap();
    assert_eq!(
        refit.document["scientific_validation"]["status"],
        "deferred"
    );
    // Every original 873 K Table 4 row; new fit is not a correction to the source curve.
    let rows = [
        (36.081, 8.75),
        (35.276, 10.46),
        (34.749, 12.44),
        (34.280, 14.0),
        (33.890, 15.48),
        (37.695, 5.21),
        (36.593, 7.52),
        (35.848, 9.49),
    ];
    let sse: f64 = rows
        .iter()
        .map(|(v, p)| (refit.pressure(*v, 873.0).unwrap() - p).powi(2))
        .sum();
    assert!((sse / 8.0).sqrt() < 0.16);
    let continuous = material
        .record("brucite_peritheos_fukui_2003_thermal_refit")
        .unwrap();
    for (p, t) in [(5.95, 300.0), (10.0, 573.0), (15.48, 873.0)] {
        let v = continuous.volume(p, t).unwrap();
        assert!((continuous.pressure(v, t).unwrap() - p).abs() < 1e-9);
    }
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}

#[test]
fn brucite_staged_isotherms_execute_and_preserve_pending_assessments() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/brucite-staged.eosmat");
    let material = load_eosmat(path).unwrap();
    assert_eq!(material.document["formula_units_per_cell"], 1);
    assert_eq!(material.document["space_group_number"], 164);
    for temperature in [300, 473, 673, 873] {
        let identifier = format!("brucite_fukui_2003_{temperature}k_bm3");
        let record = material.record(&identifier).unwrap();
        assert_eq!(
            record.document["scientific_validation"]["status"],
            "pending_primary_source_check"
        );
        let v0 = record.document["eos"]["parameters"]["V0"].as_f64().unwrap();
        let t = f64::from(temperature);
        assert!(record.pressure(v0, t).unwrap().abs() < 1e-12);
        for pressure in [1.0, 6.0, 15.0, 22.0] {
            let volume = record.volume(pressure, t).unwrap();
            assert!((record.pressure(volume, t).unwrap() - pressure).abs() < 1e-9);
        }
    }
    let hot = material.record("brucite_fukui_2003_473k_bm3").unwrap();
    // Fukui Table 4, p. 514, row 19; pressure error is 0.03 GPa.
    assert!((hot.pressure(36.883, 473.0).unwrap() - 5.95).abs() < 0.03);
    let round_trip = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(round_trip.document, material.document);
}

#[test]
fn brucite_curved_expansion_preserves_moduli_and_maps_reference_alpha() {
    let path =
        Path::new(env!("CARGO_MANIFEST_DIR")).join("tests/data/brucite-diagnostic-archive.eosmat");
    let material = load_eosmat(path).unwrap();
    for (t, k0, kp) in [
        (300, 41.8, 6.6),
        (473, 32.6, 8.3),
        (673, 28.8, 9.0),
        (873, 23.9, 10.0),
    ] {
        let identifier = format!("brucite_peritheos_fukui_2003_curved_expansion_{t}k_bm3_refit");
        let r = material.record(&identifier).unwrap();
        assert_eq!(r.document["eos"]["parameters"]["K0"], k0);
        assert_eq!(r.document["eos"]["parameters"]["K0_prime"], kp);
        assert_eq!(r.document["record_kind"], "refit");
        assert_eq!(r.document["scientific_validation"]["status"], "deferred");
        let volume = r.volume(10.0, f64::from(t)).unwrap();
        assert!((r.pressure(volume, f64::from(t)).unwrap() - 10.0).abs() < 1e-9);
    }
    let thermal = material
        .record("brucite_peritheos_fukui_2003_curved_expansion_thermal_refit")
        .unwrap();
    assert_eq!(
        thermal.document["thermal"]["thermal_expansion_law"],
        "linear_reference_temperature"
    );
    let delta: f64 = 573.0;
    let v0 = 40.926_216_959_111_46
        * (6.066_529_28e-5 * delta + 0.5 * 4.389_573_89e-8 * delta * delta).exp();
    let eta = (v0 / 35.5).powf(1.0 / 3.0);
    let expected = 1.5
        * (41.8 - 0.0319 * delta)
        * (eta.powi(7) - eta.powi(5))
        * (1.0 + 0.75 * (6.6 + 0.00571 * delta - 4.0) * (eta.powi(2) - 1.0));
    assert!((thermal.pressure(35.5, 873.0).unwrap() - expected).abs() < 1e-8);
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
