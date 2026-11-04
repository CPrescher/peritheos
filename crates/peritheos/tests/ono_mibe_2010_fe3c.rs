use peritheos::{load_eosmat, load_eosmat_str};
use std::path::Path;

#[test]
fn original_ono_mibe_experiment_is_separate_and_finite_in_scope() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/cementite.eosmat");
    if !path.is_file() {
        return;
    }
    let material = load_eosmat(path).unwrap();
    let r = material.record("cementite_ono_mibe_2010_bm3").unwrap();
    assert_eq!(r.reference_temperature, 300.0);
    assert!((r.reference_volume() - 155.4).abs() < 1e-12);
    assert_eq!(
        material.default_record().unwrap().identifier,
        "cementite_scott_2001_bm3_1"
    );
    assert!((r.volume(30.0, 300.0).unwrap() - 137.7551).abs() < 0.41);
    for p in [0.0, 10.0, 20.0, 34.3] {
        let v = r.volume(p, 300.0).unwrap();
        assert!((r.pressure(v, 300.0).unwrap() - p).abs() < 1e-8);
    }
    let raw = material.document["eos_records"]
        .as_array()
        .unwrap()
        .iter()
        .find(|v| v["identifier"] == r.identifier)
        .unwrap();
    assert_eq!(raw["parameters"]["K0"]["value"], 167);
    assert_eq!(raw["parameters"]["K0"]["error"], serde_json::Value::Null);
    assert_eq!(
        raw["validity"]["pressure_gpa"],
        serde_json::json!([0, 34.3])
    );
    assert_eq!(
        raw["evidence"]["reproduction"]["assessment"]["scope"],
        "source_equation"
    );
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
