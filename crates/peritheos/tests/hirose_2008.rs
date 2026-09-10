use peritheos::isothermal::BM3;
use peritheos::thermal::{ReferenceVolumeLaw, ThermalExpansionLaw, ThermalReferenceState};
use peritheos::ThermalEos;

#[test]
fn hirose_preferred_fit_reproduces_table1_and_reference_state() {
    let model = ThermalReferenceState::new(
        BM3::new(67.85, 167.0, 5.58).unwrap(),
        300.0,
        3.824e-5,
        0.0,
        1.499e-8,
        ThermalExpansionLaw::LinearTemperature,
        ReferenceVolumeLaw::IntegratedExpansivity,
    )
    .unwrap()
    .with_temperature_laws(Some([1.03e-6, 3.95e-10, 1.61e-13]), 3.61e-4)
    .unwrap();
    assert!((model.pressure(3.7152_f64.powi(3), 2070.0).unwrap() - 119.7).abs() < 0.05);
    assert!(model.thermal_pressure(50.0, 300.0).unwrap().abs() < 1e-12);
    assert!((model.bulk_modulus(67.85, 300.0, 1e-6).unwrap() - 167.0).abs() < 1e-10);
    let invalid = model
        .with_temperature_laws(Some([-0.01, 0.0, 0.0]), 0.0)
        .unwrap();
    assert!(invalid.pressure(50.0, 400.0).is_err());
}

#[test]
fn hirose_gold_records_load_and_roundtrip_through_native_eosmat() {
    use peritheos::eosmat::load_eosmat_str;
    let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/gold.eosmat");
    let mut document: serde_json::Value =
        serde_json::from_str(&std::fs::read_to_string(path).unwrap()).unwrap();
    document["eos_records"]
        .as_array_mut()
        .unwrap()
        .retain(|record| {
            record["identifier"]
                .as_str()
                .unwrap()
                .starts_with("gold_hirose_2008_")
        });
    document["datasets"]
        .as_array_mut()
        .unwrap()
        .retain(|dataset| dataset["identifier"].as_str().unwrap() == "gold_hirose_2008_table1");
    assert_eq!(document["eos_records"].as_array().unwrap().len(), 3);
    let material = load_eosmat_str(&document.to_string()).unwrap();
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(material.document, restored.document);
    let record = material
        .eos_records
        .iter()
        .find(|r| r.identifier == "gold_hirose_2008_bm3_fit2")
        .unwrap();
    assert!((record.pressure(3.7152_f64.powi(3), 2070.0).unwrap() - 119.7).abs() < 0.05);
}
