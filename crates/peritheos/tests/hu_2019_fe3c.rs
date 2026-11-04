use peritheos::{load_eosmat, load_eosmat_str};
use std::path::Path;

#[test]
fn hu_original_hugoniot_preserves_porous_precursor_and_full_experimental_domain() {
    let path = Path::new(env!("CARGO_MANIFEST_DIR"))
        .join("../../peritheos/data/materials/cementite_nonmagnetic.eosmat");
    if !path.is_file() {
        return;
    }
    let material = load_eosmat(path).unwrap();
    assert_eq!(
        material.default_record().unwrap().identifier,
        "cementite_nonmagnetic_takahashi_2019_kbr_bm3_mgd"
    );
    let r = material.default_hugoniot_record().unwrap();
    assert_eq!(
        r.record.identifier,
        "cementite_nonmagnetic_hu_2019_linear_hugoniot"
    );
    let state = r.state_from_particle_velocity(2.164).unwrap();
    assert!((state.pressure - 132.221_941_633_6).abs() < 1e-8);
    assert!((state.density - 10.259_206_044_1).abs() < 1e-8);
    for up in [1.585, 2.164, 2.571, 2.955, 3.138, 3.286] {
        assert!(r.state_from_particle_velocity(up).is_ok());
    }
    let v = r.record.volume(245.0, 300.0).unwrap();
    assert!((r.record.pressure(v, 300.0).unwrap() - 245.0).abs() < 1e-8);
    assert!(r.state_from_particle_velocity(3.4).is_err());
    let restored = load_eosmat_str(&material.to_json().unwrap()).unwrap();
    assert_eq!(restored.document, material.document);
}
