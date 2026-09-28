use peritheos::{isothermal::OddInversePower, load_eosmat, IsothermalEos};

#[test]
fn published_argon_series_matches_molar_equation_and_inverts() {
    let factor = 4.0e24 / 6.022_140_76e23;
    let material = load_eosmat(concat!(
        env!("CARGO_MANIFEST_DIR"),
        "/tests/data/argon_fcc.eosmat"
    ))
    .unwrap();
    let table = [
        (
            "4p2",
            4.2,
            22.56,
            [-181_890.0, 101_480_000.0, -4_544_000_000.0, 0.0],
        ),
        (
            "20p0",
            20.0,
            22.64,
            [-174_680.0, 97_126_000.0, -3_884_000_000.0, 0.0],
        ),
        (
            "40p0",
            40.0,
            23.04,
            [-151_020.0, 83_928_000.0, -2_014_000_000.0, 0.0],
        ),
        (
            "60p0",
            60.0,
            23.63,
            [-117_590.0, 62_868_000.0, 1_550_000_000.0, 0.0],
        ),
        (
            "77p0",
            77.0,
            24.28,
            [
                -57_020.0,
                -1_197_000.0,
                25_851_000_000.0,
                -3_135_200_000_000.0,
            ],
        ),
    ];
    for (id, t, v0, a) in table {
        let c: Vec<f64> = a
            .iter()
            .zip([3, 5, 7, 9])
            .map(|(a, n)| a / (10.0 * f64::powi(v0, n)))
            .collect();
        let model = OddInversePower::new(v0 * factor, c[0], c[1], c[2], c[3]).unwrap();
        let record = material
            .record(&format!("argon_fcc_anderson_swenson_1975_{id}k"))
            .unwrap();
        for molar in [17.5, 19.0, 21.0, v0] {
            let v = molar * factor;
            let expected: f64 = a
                .iter()
                .zip([3, 5, 7, 9])
                .map(|(a, n)| a / f64::powi(molar, n) / 10.0)
                .sum();
            let p = model.pressure(v).unwrap();
            assert!((p - expected).abs() < 1e-13);
            assert!((record.pressure(v, t).unwrap() - expected).abs() < 1e-13);
            assert!((model.volume(p).unwrap() - v).abs() < 1e-8);
            let h = v * 1e-5;
            let numerical =
                -v * (model.pressure(v + h).unwrap() - model.pressure(v - h).unwrap()) / (2.0 * h);
            assert!((model.bulk_modulus(v).unwrap() - numerical).abs() < 1e-6);
        }
    }
    assert!(OddInversePower::new(0.0, 1.0, 0.0, 0.0, 0.0).is_err());
}
