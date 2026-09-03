use peritheos::{
    diamond_raman_calibration, recalculate_diamond_raman_pressure, recalculate_ruby_pressure,
    ruby_calibration, DiamondRamanCalibrationModel, RubyCalibrationModel,
    DIAMOND_RAMAN_CALIBRATIONS, RUBY_CALIBRATIONS,
};

#[test]
fn bundled_ruby_scales_have_stable_identifiers_and_models() {
    assert_eq!(RUBY_CALIBRATIONS.len(), 6);
    assert_eq!(RUBY_CALIBRATIONS[0].identifier, "ruby_mao_1978");
    assert!(matches!(
        ruby_calibration("ruby_dorogokupets_oganov_2007")
            .expect("bundled scale")
            .model,
        RubyCalibrationModel::QuadraticShift { .. }
    ));
    assert!(ruby_calibration("unknown").is_none());
    assert!(matches!(
        ruby_calibration("ruby_shen_2020")
            .expect("bundled scale")
            .model,
        RubyCalibrationModel::QuadraticShift { .. }
    ));
}

#[test]
fn bundled_diamond_raman_scales_are_executable_and_reversible() {
    assert_eq!(DIAMOND_RAMAN_CALIBRATIONS.len(), 2);
    assert!(matches!(
        diamond_raman_calibration("diamond_raman_akahama_2006")
            .expect("bundled scale")
            .model,
        DiamondRamanCalibrationModel::AkahamaQuadratic { .. }
    ));
    assert!(diamond_raman_calibration("unknown").is_none());

    for scale in DIAMOND_RAMAN_CALIBRATIONS {
        for pressure in [0.0, 50.0, 150.0, 300.0] {
            let ratio = scale.wavenumber_ratio(pressure).expect("valid ratio");
            let recovered = scale.pressure_from_ratio(ratio).expect("valid pressure");
            assert!(
                (recovered - pressure).abs() < 1.0e-9,
                "{}",
                scale.identifier
            );
        }
    }

    let source = diamond_raman_calibration("diamond_raman_akahama_2006").unwrap();
    let target = diamond_raman_calibration("diamond_raman_eremets_2023").unwrap();
    let converted = recalculate_diamond_raman_pressure(150.0, source, target).unwrap();
    let recovered = recalculate_diamond_raman_pressure(converted, target, source).unwrap();
    assert!((recovered - 150.0).abs() < 1.0e-9);
}

#[test]
fn every_ruby_scale_round_trips_pressure_and_wavelength_ratio() {
    for scale in RUBY_CALIBRATIONS {
        for pressure in [0.0, 1.0, 50.0, 150.0, 300.0] {
            let ratio = scale.wavelength_ratio(pressure).expect("valid ratio");
            let recovered = scale.pressure_from_ratio(ratio).expect("valid pressure");
            assert!(
                (recovered - pressure).abs() < 1.0e-9,
                "{}",
                scale.identifier
            );
        }
    }
}

#[test]
fn ruby_pressure_can_be_recalculated_between_published_scales() {
    let source = ruby_calibration("ruby_mao_1986").unwrap();
    let target = ruby_calibration("ruby_dorogokupets_oganov_2007").unwrap();
    let converted = recalculate_ruby_pressure(100.0, source, target).unwrap();
    assert!(converted > 100.0);

    let recovered = recalculate_ruby_pressure(converted, target, source).unwrap();
    assert!((recovered - 100.0).abs() < 1.0e-9);
}

#[test]
fn ruby_scale_rejects_unphysical_inputs() {
    let scale = ruby_calibration("ruby_mao_1986").unwrap();
    assert!(scale.pressure_from_ratio(0.99).is_err());
    assert!(scale.wavelength_ratio(-1.0).is_err());
    assert!(scale.pressure_from_wavelength(0.0).is_err());

    let wavelength = scale.wavelength_from_pressure(50.0).unwrap();
    assert!((scale.pressure_from_wavelength(wavelength).unwrap() - 50.0).abs() < 1.0e-9);

    let holzapfel = ruby_calibration("ruby_holzapfel_2005").unwrap();
    assert!(holzapfel.wavelength_ratio(2_000.0).is_err());
}
