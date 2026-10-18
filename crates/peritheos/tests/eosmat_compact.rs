use peritheos::{export_eosmat_v3, material_from_value};
use serde_json::{json, Value};

fn sample() -> Value {
    json!({"format_version":4,"identifier":"example","formula":"X",
        "shared_metadata": {
            "source": {"authors":["Example"],"year":2026,"source":"Fictional teaching source"},
            "alias": {"shared_ref":"source"},
            "unused": {"note":"Preserve unused supporting evidence"}
        },
        "eos_records":[{"identifier":"example_bm3","model":"BM3","temperature_ref":300,
            "parameters":{"V0":10,"K0":160,"K0_prime":4},
            "reference":{"shared_ref":"alias"},
            "evidence":{"reproduction":{"assessment":{"status":"parity","scope":"parameter_agreement",
                "basis":"Synthetic test only","qualifications":["No scientific validation claimed"]}}}}],
        "application_extension":{"shared_ref":"literal_not_a_reference"}})
}

#[test]
fn local_metadata_is_explicit_preserved_and_numerically_equivalent() {
    let source = sample();
    let native = material_from_value(source.clone()).unwrap();
    assert_eq!(native.document, source);
    let legacy = export_eosmat_v3(&source).unwrap();
    assert_eq!(
        legacy["eos_records"][0]["reference"],
        source["shared_metadata"]["source"]
    );
    assert_eq!(
        legacy["shared_metadata"]["unused"],
        source["shared_metadata"]["unused"]
    );
    assert_eq!(
        legacy["application_extension"],
        source["application_extension"]
    );
    let old = material_from_value(legacy).unwrap();
    assert_eq!(
        old.eos_records[0].pressure(9., 300.).unwrap(),
        native.eos_records[0].pressure(9., 300.).unwrap()
    );
}

#[test]
fn missing_cyclic_or_overridden_references_fail_including_unused_definitions() {
    for (pool, reference) in [
        (json!({}), json!({"shared_ref":"missing"})),
        (json!({"a":{"shared_ref":"a"}}), json!({"shared_ref":"a"})),
        (
            json!({"a":{"shared_ref":"b"},"b":{"shared_ref":"a"}}),
            json!({"shared_ref":"a"}),
        ),
        (
            json!({"a":{}}),
            json!({"shared_ref":"a","title":"Hidden override"}),
        ),
        (json!({"a":{}}), json!({"shared_ref":["a"]})),
        (json!({"bad-name":{}}), json!("Ordinary reference")),
        (json!({"a":"not an object"}), json!("Ordinary reference")),
    ] {
        let mut document = sample();
        document["shared_metadata"] = pool;
        document["eos_records"][0]["reference"] = reference;
        assert!(material_from_value(document).is_err());
    }
}

#[test]
fn assessment_keeps_reproduction_separate_from_source_validation() {
    for (key, value) in [
        ("status", json!("primary_source_validated")),
        ("scope", json!("unknown")),
        ("basis", json!("")),
        ("qualifications", json!([null])),
    ] {
        let mut document = sample();
        document["eos_records"][0]["evidence"]["reproduction"]["assessment"][key] = value;
        assert!(material_from_value(document).is_err());
    }
}

#[test]
fn documented_thermal_and_curated_examples_use_the_same_native_reader() {
    let directory = std::path::Path::new(env!("CARGO_MANIFEST_DIR")).join("../../docs/examples");
    if !directory.exists() {
        return;
    }
    for name in [
        "thermal-material.eosmat",
        "curated-material.eosmat",
        "shared-material.eosmat",
    ] {
        let value: Value =
            serde_json::from_slice(&std::fs::read(directory.join(name)).unwrap()).unwrap();
        let native = material_from_value(value.clone()).unwrap();
        let legacy = material_from_value(export_eosmat_v3(&value).unwrap()).unwrap();
        for (a, b) in native.eos_records.iter().zip(&legacy.eos_records) {
            assert_eq!(
                a.pressure(70., 300.).unwrap(),
                b.pressure(70., 300.).unwrap()
            );
        }
    }
}
