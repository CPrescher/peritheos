use std::path::Path;

use peritheos::{load_eosmat, load_eosmat_str};

#[test]
fn bundled_phase_e_loads_executes_and_preserves_primary_data_in_rust() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/phase_e_fe0072.eosmat");
    if !path.is_file() {
        return; // Bundled data are not shipped with the standalone crate.
    }
    let material = load_eosmat(path).unwrap();
    let record = material
        .record("phase_e_fe0072_crichton_ross_2000_bm3")
        .unwrap();
    assert!((record.pressure(100.0, 298.0).unwrap() - 6.196_492_213_716_72).abs() < 1e-11);
    assert!(record.pressure(105.627, 298.0).unwrap().abs() < 1e-12);
    for pressure in [0.0, 0.296, 2.973, 6.671] {
        let volume = record.volume(pressure, 298.0).unwrap();
        assert!((record.pressure(volume, 298.0).unwrap() - pressure).abs() < 1e-9);
    }
    assert_eq!(material.document["formula"], "Mg1.96Fe0.072Si1.04H3.7O6");
    assert!(material.document["datasets"][0]["rows"][0][1].is_null());
    assert_eq!(
        material.document["datasets"][0]["rows"]
            .as_array()
            .unwrap()
            .len(),
        11
    );
    assert_eq!(material.document["peaks"].as_array().unwrap().len(), 8);
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
