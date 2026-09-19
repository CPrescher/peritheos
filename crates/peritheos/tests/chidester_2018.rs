#[test]
fn tho2_experimental_records_preserve_native_interchange_metadata() {
    for name in ["thorianite", "tho2_cotunnite"] {
        let path = std::path::Path::new(env!("CARGO_MANIFEST_DIR"))
            .join(format!("../../peritheos/data/materials/{name}.eosmat"));
        if !path.exists() {
            continue;
        }
        let loaded = peritheos::load_eosmat(path).unwrap();
        let restored = peritheos::load_eosmat_str(&loaded.to_json().unwrap()).unwrap();
        assert_eq!(loaded.document, restored.document);
        assert_eq!(loaded.document["formula_units_per_cell"], 4);
        assert_eq!(
            loaded.document["datasets"][0]["columns"][2]["role"],
            "uncertainty"
        );
    }
}
