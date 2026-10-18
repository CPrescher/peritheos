use peritheos::{export_eosmat_v3, material_from_value};
use serde_json::{json, Value};

fn sample() -> Value {
    json!({"format_version":4,"identifier":"example","formula":"X",
    "eos_records":[{"identifier":"example_bm3","model":"BM3","temperature_ref":300,
        "parameters":{"V0":10,"K0":160,"K0_prime":4},
        "evidence":{"reproduction":{"assessment":{
            "status":"not_reproduced","scope":"parameter_agreement",
            "basis":"Synthetic comparison for parser tests",
                "report_url":"https://example.org/report/#record-example_bm3",
                "report_path":"docs/report.md#record-example_bm3",
            "parameters":[{"parameter":"K0","published":160,"published_error":null,
                "refit":180.12345678901235,"refit_error":2.5,"similar":false}],
            "fit":{"fixed_parameters":["V0","K0_prime"]}
        }}}}]})
}

#[test]
fn failed_reconstruction_values_survive_exports_without_changing_the_eos() {
    let source = sample();
    let native = material_from_value(source.clone()).unwrap();
    let legacy = export_eosmat_v3(&source).unwrap();
    assert_eq!(
        legacy["eos_records"][0]["reproduction"],
        source["eos_records"][0]["evidence"]["reproduction"]
    );
    assert_eq!(legacy["eos_records"][0]["eos"]["parameters"]["K0"], 160);
    let restored = material_from_value(legacy).unwrap();
    let mut without_assessment = source;
    without_assessment["eos_records"][0]
        .as_object_mut()
        .unwrap()
        .remove("evidence");
    let original = material_from_value(without_assessment).unwrap();
    for volume in [9., 8., 7.] {
        let expected = original.eos_records[0].pressure(volume, 300.).unwrap();
        assert_eq!(
            native.eos_records[0].pressure(volume, 300.).unwrap(),
            expected
        );
        assert_eq!(
            restored.eos_records[0].pressure(volume, 300.).unwrap(),
            expected
        );
    }
}

#[test]
fn invalid_reconstructed_parameter_tables_are_rejected() {
    for parameters in [
        json!({}),
        json!([{"parameter":"K0","published":160}]),
        json!([{"parameter":"","published":160,"refit":180}]),
        json!([{"parameter":"K0","published":true,"refit":180}]),
        json!([{"parameter":"K0","published":160,"refit":180,"refit_error":-1}]),
        json!([{"parameter":"K0","published":160,"refit":180},
               {"parameter":"K0","published":160,"refit":180}]),
    ] {
        let mut source = sample();
        source["eos_records"][0]["evidence"]["reproduction"]["assessment"]["parameters"] =
            parameters;
        assert!(material_from_value(source).is_err());
    }
}

#[test]
fn optional_report_links_require_nonempty_strings() {
    for key in ["report_url", "report_path"] {
        for value in [json!(null), json!(""), json!(123), json!([])] {
            let mut source = sample();
            source["eos_records"][0]["evidence"]["reproduction"]["assessment"][key] = value;
            assert!(material_from_value(source).is_err());
        }
    }
}
