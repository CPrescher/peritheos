//! EOSMAT v5 record `check` headlines and audit sidecars.
//!
//! The bundled-catalog tests mirror `tests/test_eosmat_audit_split.py`: the
//! Rust headline of every record must equal the stored `check` block and the
//! derivation from the merged audit, and loading the record file without its
//! sidecar must not change any evaluated pressure.

use std::fs;
use std::path::{Path, PathBuf};

use peritheos::{
    derive_eosmat_check, export_eosmat_v3, load_eosmat, load_eosmat_str,
    load_eosmat_str_with_audit, load_eosmat_without_audit, material_from_value, merge_eosmat_audit,
    validate_eosmat_document, EosmatErrorKind, Material, Reproduction, AUDIT_RECORD_KEYS,
};
use serde_json::{json, Value};

fn crate_directory() -> &'static Path {
    Path::new(env!("CARGO_MANIFEST_DIR"))
}

fn catalog_directory() -> Option<PathBuf> {
    let path = crate_directory().join("../../peritheos/data/materials");
    if path.is_dir() {
        Some(path)
    } else {
        eprintln!("skipping bundled-catalog test: the catalog is not packaged with the crate");
        None
    }
}

fn bundled_paths(directory: &Path) -> Vec<PathBuf> {
    let mut paths: Vec<PathBuf> = fs::read_dir(directory)
        .unwrap()
        .map(|entry| entry.unwrap().path())
        .filter(|path| {
            path.extension()
                .is_some_and(|extension| extension == "eosmat")
        })
        .collect();
    paths.sort();
    paths
}

fn without_audit_link(check: &Value) -> Value {
    let mut check = check.clone();
    check.as_object_mut().unwrap().remove("audit");
    check
}

/// Evaluate pressures at fixed fractions of the reference volume.
///
/// Errors (for example outside a Hugoniot branch domain) are compared by
/// their debug form, so a changed error is also detected.
fn pressure_samples(material: &Material) -> Vec<String> {
    let mut samples = Vec::new();
    for record in &material.eos_records {
        for fraction in [1.0, 0.95, 0.85] {
            let volume = record.reference_volume() * fraction;
            let result = record.pressure(volume, record.reference_temperature);
            samples.push(format!("{} {fraction} {result:?}", record.identifier));
        }
    }
    samples
}

fn assert_headline_matches_check(material: &Material, path: &Path) {
    for record in &material.eos_records {
        let check = record
            .check()
            .unwrap_or_else(|| panic!("{}: {} has no check", path.display(), record.identifier));
        assert_eq!(
            Some(record.reproduction().as_str()),
            check["reproduction"].as_str(),
            "{}: {}",
            path.display(),
            record.identifier
        );
        assert_eq!(
            record.warning(),
            check["warning"].as_str(),
            "{}: {}",
            path.display(),
            record.identifier
        );
    }
}

#[test]
fn bundled_headlines_equal_stored_checks_and_the_merged_audit_derivation() {
    let Some(directory) = catalog_directory() else {
        return;
    };
    let paths = bundled_paths(&directory);
    assert!(paths.len() > 100, "unexpectedly small catalog");
    let mut records = 0;
    for path in &paths {
        let merged = load_eosmat(path).unwrap_or_else(|error| {
            panic!("{} failed to load with its audit: {error}", path.display())
        });
        let record_file = load_eosmat_without_audit(path).unwrap_or_else(|error| {
            panic!(
                "{} failed to load without its audit: {error}",
                path.display()
            )
        });
        assert_headline_matches_check(&merged, path);
        assert_headline_matches_check(&record_file, path);
        for (complete, alone) in merged.eos_records.iter().zip(&record_file.eos_records) {
            records += 1;
            let check = complete.check().unwrap();
            assert_eq!(
                complete.audit_reference(),
                Some(
                    format!(
                        "../audits/{}.audit.json#{}",
                        merged.identifier, complete.identifier
                    )
                    .as_str()
                ),
            );
            assert_eq!(
                derive_eosmat_check(&complete.document, Some(&merged.document)),
                without_audit_link(check),
                "{}: {} check differs from its merged evidence",
                path.display(),
                complete.identifier
            );
            assert!(!complete.audit_is_external());
            assert!(complete.audit().is_some());
            assert_eq!(alone.reproduction(), complete.reproduction());
            assert_eq!(alone.warning(), complete.warning());
            assert_eq!(alone.check(), Some(check));
            // Every bundled record file names a sidecar entry; records with
            // no audit-only block have no sidecar entry and stay complete.
            assert_eq!(alone.audit_is_external(), alone.audit().is_none());
        }
        assert_eq!(
            pressure_samples(&merged),
            pressure_samples(&record_file),
            "{}: pressures depend on whether the audit is merged",
            path.display()
        );
    }
    assert!(records > 300, "unexpectedly few bundled EOS records");
}

#[test]
fn bundled_v3_exports_restore_evidence_and_drop_check() {
    let Some(directory) = catalog_directory() else {
        return;
    };
    for path in bundled_paths(&directory) {
        let material = load_eosmat(&path).unwrap();
        let legacy = export_eosmat_v3(&material.document).unwrap();
        let exported = material_from_value(legacy.clone()).unwrap();
        for (record, source) in legacy["eos_records"]
            .as_array()
            .unwrap()
            .iter()
            .zip(&material.eos_records)
        {
            assert!(record.get("check").is_none(), "{}", path.display());
            assert!(record.get("evidence").is_none(), "{}", path.display());
            assert_eq!(
                record["scientific_validation"]["status"],
                source.check().unwrap()["source"],
                "{}: {}",
                path.display(),
                source.identifier
            );
        }
        // Without check, the v3 headline is derived from the inline evidence.
        for (old, new) in exported.eos_records.iter().zip(&material.eos_records) {
            assert_eq!(old.reproduction(), new.reproduction(), "{}", path.display());
            assert_eq!(old.warning(), new.warning(), "{}", path.display());
        }
        assert_eq!(pressure_samples(&exported), pressure_samples(&material));
    }
}

#[test]
fn packaged_fixtures_merge_sidecars_relative_to_the_record_file() {
    let fixtures = crate_directory().join("tests/data");
    for name in [
        "brucite",
        "diopside",
        "forsterite",
        "ferrous_ringwoodite",
        "nacl_b2",
        "palladium",
        "rbcl",
        "rbbr",
        "rbi",
    ] {
        let path = fixtures.join(format!("{name}.eosmat"));
        let material = load_eosmat(&path).unwrap();
        assert_headline_matches_check(&material, &path);
        assert!(material
            .eos_records
            .iter()
            .all(|record| !record.audit_is_external()));
        if let Some(directory) = catalog_directory() {
            assert_eq!(
                fs::read(&path).unwrap(),
                fs::read(directory.join(format!("{name}.eosmat"))).unwrap(),
                "{name} fixture differs from the catalog record file"
            );
            assert_eq!(
                fs::read(crate_directory().join(format!("tests/audits/{name}.audit.json")))
                    .unwrap(),
                fs::read(directory.join(format!("../audits/{name}.audit.json"))).unwrap(),
                "{name} fixture sidecar differs from the catalog sidecar"
            );
        }
    }
}

#[test]
fn v4_fixture_derives_the_headlines_stored_in_the_v5_catalog() {
    // tests/data/argon_fcc.eosmat is the catalog's v4 compatibility export:
    // full evidence inline and no check block.
    let path = crate_directory().join("tests/data/argon_fcc.eosmat");
    let legacy = load_eosmat(&path).unwrap();
    assert_eq!(legacy.document["format_version"], 4);
    assert!(legacy
        .eos_records
        .iter()
        .all(|record| record.check().is_none()));
    let Some(directory) = catalog_directory() else {
        return;
    };
    let current = load_eosmat(directory.join("argon_fcc.eosmat")).unwrap();
    assert_eq!(legacy.eos_records.len(), current.eos_records.len());
    for (old, new) in legacy.eos_records.iter().zip(&current.eos_records) {
        assert_eq!(old.identifier, new.identifier);
        assert_eq!(old.reproduction(), new.reproduction(), "{}", old.identifier);
        assert_eq!(old.warning(), new.warning(), "{}", old.identifier);
        assert_eq!(
            derive_eosmat_check(&old.document, Some(&legacy.document)),
            without_audit_link(new.check().unwrap()),
            "{}",
            old.identifier
        );
    }
    assert_eq!(pressure_samples(&legacy), pressure_samples(&current));
}

fn standalone(extra: &Value) -> Value {
    let mut record = json!({
        "identifier": "example_bm3",
        "model": "BM3",
        "temperature_ref": 300,
        "parameters": {"V0": 74.7, "K0": 160, "K0_prime": 4},
        "reference": "Fictional example; no experimental validation."
    });
    record
        .as_object_mut()
        .unwrap()
        .extend(extra.as_object().unwrap().clone());
    json!({
        "format_version": 5,
        "identifier": "example",
        "formula": "MgO",
        "eos_records": [record]
    })
}

fn evidence(status: &str) -> Value {
    json!({
        "scientific_validation": {
            "status": "primary_source_validated",
            "audit_date": "2026-01-02",
            "primary_source_check": {"finding": "Table 2 values match."}
        },
        "reproduction": {
            "assessment": {
                "status": status,
                "scope": "parameter_agreement",
                "basis": "Independent refit of Table 1."
            }
        }
    })
}

fn invalid_reason(document: &Value) -> String {
    let error = validate_eosmat_document(document).unwrap_err();
    assert_eq!(error.kind(), EosmatErrorKind::InvalidDocument, "{error}");
    error.to_string()
}

#[test]
fn standalone_files_need_no_check_and_no_assessment_is_untested() {
    let plain = material_from_value(standalone(&json!({}))).unwrap();
    let record = &plain.eos_records[0];
    assert_eq!(record.reproduction(), Reproduction::Untested);
    assert_eq!(record.warning(), None);
    assert!(record.check().is_none());
    assert!(!record.audit_is_external());
    assert_eq!(record.audit(), Some(serde_json::Map::new()));

    for (status, expected) in [
        ("parity", Reproduction::Reproduced),
        ("similar", Reproduction::Close),
        ("not_reproduced", Reproduction::Disagrees),
        ("not_refittable", Reproduction::Untested),
        ("source_reconstruction", Reproduction::Untested),
        ("not_assessed", Reproduction::Untested),
    ] {
        let inline = material_from_value(standalone(&json!({"evidence": evidence(status)})));
        let record = &inline.unwrap().eos_records[0];
        assert_eq!(record.reproduction(), expected, "{status}");
        assert_eq!(record.audit().unwrap()["evidence"], evidence(status));
    }
}

#[test]
fn inline_check_must_agree_with_inline_evidence() {
    let document = standalone(&json!({"evidence": evidence("parity")}));
    let derived = derive_eosmat_check(&document["eos_records"][0], Some(&document));
    assert_eq!(
        derived,
        json!({
            "reproduction": "reproduced",
            "warning": null,
            "finding": "Table 2 values match.",
            "source": "primary_source_validated",
            "source_checked_on": "2026-01-02"
        })
    );
    let consistent = standalone(&json!({"evidence": evidence("parity"), "check": derived}));
    let material = material_from_value(consistent.clone()).unwrap();
    assert_eq!(
        material.eos_records[0].reproduction(),
        Reproduction::Reproduced
    );
    assert!(export_eosmat_v3(&consistent).unwrap()["eos_records"][0]
        .get("check")
        .is_none());

    let contradiction = standalone(&json!({
        "evidence": evidence("parity"),
        "check": {"reproduction": "disagrees"}
    }));
    assert!(invalid_reason(&contradiction).contains("check.reproduction"));
    assert!(material_from_value(contradiction).is_err());

    let source = standalone(&json!({
        "evidence": evidence("similar"),
        "check": {"reproduction": "close", "source": "deferred"}
    }));
    assert!(invalid_reason(&source).contains("check.source contradicts"));
}

#[test]
fn check_block_structure_is_validated() {
    for (check, reason) in [
        (json!([]), "check must be an object"),
        (json!({"status": "parity"}), "unknown fields"),
        (
            json!({"reproduction": "parity"}),
            "check.reproduction must be one of",
        ),
        (json!({"source": "validated"}), "check.source is invalid"),
        (
            json!({"warning": " "}),
            "check.warning must be null or text",
        ),
        (json!({"warning": 3}), "check.warning must be null or text"),
        (
            json!({"audit": ""}),
            "check.audit must be a non-empty string",
        ),
        (
            json!({"finding": 1}),
            "check.finding must be a non-empty string",
        ),
    ] {
        let document = standalone(&json!({ "check": check }));
        assert!(invalid_reason(&document).contains(reason), "{check}");
    }
    let document = standalone(&json!({"check": {"reproduction": "close", "warning": null}}));
    let record = &material_from_value(document).unwrap().eos_records[0];
    assert_eq!(record.reproduction(), Reproduction::Close);
    assert_eq!(record.warning(), None);
}

fn split_example() -> (Value, Value) {
    let mut audit_evidence = evidence("not_reproduced");
    audit_evidence["scientific_validation"]["usage_recommendation"] = json!("not_recommended");
    let record_file = standalone(&json!({
        "check": {
            "reproduction": "disagrees",
            "warning": "Not recommended for use.",
            "finding": "Table 2 values match.",
            "source": "primary_source_validated",
            "source_checked_on": "2026-01-02",
            "audit": "../audits/example.audit.json#example_bm3"
        }
    }));
    let sidecar = json!({
        "format": "peritheos.eosmat-audit",
        "format_version": 5,
        "material": "example",
        "records": {"example_bm3": {
            "evidence": audit_evidence,
            "scope_history": ["Fictional scope note."]
        }}
    });
    (record_file, sidecar)
}

#[test]
fn record_files_load_without_a_sidecar_and_merge_with_one() {
    let (record_file, sidecar) = split_example();
    let text = serde_json::to_string(&record_file).unwrap();

    let alone = load_eosmat_str(&text).unwrap();
    let record = &alone.eos_records[0];
    assert!(record.audit_is_external());
    assert_eq!(record.audit(), None);
    assert_eq!(record.reproduction(), Reproduction::Disagrees);
    assert_eq!(record.warning(), Some("Not recommended for use."));
    assert_eq!(
        record.audit_reference(),
        Some("../audits/example.audit.json#example_bm3")
    );
    // The v3 export of an unmerged record keeps check so the audit link is
    // not lost; its scientific_validation status is the stored headline.
    let legacy = export_eosmat_v3(&alone.document).unwrap();
    assert_eq!(
        legacy["eos_records"][0]["check"]["audit"],
        "../audits/example.audit.json#example_bm3"
    );
    assert_eq!(
        legacy["eos_records"][0]["scientific_validation"],
        json!({"status": "primary_source_validated"})
    );

    let merged = load_eosmat_str_with_audit(&text, &sidecar.to_string()).unwrap();
    let record = &merged.eos_records[0];
    assert!(!record.audit_is_external());
    let audit = record.audit().unwrap();
    assert_eq!(
        audit.keys().collect::<Vec<_>>(),
        ["evidence", "scope_history"]
    );
    let assessment = &record.document["evidence"]["reproduction"]["assessment"];
    assert_eq!(
        assessment["report_url"],
        "https://peritheos.readthedocs.io/en/latest/reconstruction-reports/example/#record-example_bm3"
    );
    assert_eq!(
        assessment["report_path"],
        "docs/reconstruction-reports/example.md#record-example_bm3"
    );
    assert_eq!(assessment["report"], "docs/data/primary-eos-refits.json");
    assert_eq!(
        derive_eosmat_check(&record.document, Some(&merged.document)),
        without_audit_link(record.check().unwrap())
    );
    assert_eq!(pressure_samples(&alone), pressure_samples(&merged));
    let legacy = export_eosmat_v3(&merged.document).unwrap();
    assert!(legacy["eos_records"][0].get("check").is_none());
    assert_eq!(
        legacy["eos_records"][0]["scientific_validation"],
        sidecar["records"]["example_bm3"]["evidence"]["scientific_validation"]
    );
    assert!(AUDIT_RECORD_KEYS.contains(&"scope_history"));
}

#[test]
fn audit_dependent_checks_are_deferred_until_the_sidecar_is_merged() {
    let (mut record_file, sidecar) = split_example();
    record_file["eos_records"][0]["record_kind"] = json!("refit");
    // fit_provenance is evidence, so a refit record file is valid on its own...
    let alone = load_eosmat_str(&record_file.to_string()).unwrap();
    assert!(alone.eos_records[0].audit_is_external());
    // ...but the complete document must provide it.
    let merged = merge_eosmat_audit(&record_file, &sidecar).unwrap();
    assert!(invalid_reason(&merged).contains("refit records require fit_provenance"));
    // Without check.audit the record is standalone and must carry it inline.
    let mut inline = record_file;
    inline["eos_records"][0]["check"]
        .as_object_mut()
        .unwrap()
        .remove("audit");
    assert!(invalid_reason(&inline).contains("refit records require fit_provenance"));
}

#[test]
fn sidecar_merge_rejects_mismatches() {
    let (record_file, sidecar) = split_example();
    let mut other = sidecar.clone();
    other["material"] = json!("mgo");
    assert!(merge_eosmat_audit(&record_file, &other)
        .unwrap_err()
        .to_string()
        .contains("expected \"example\""));
    let mut format = sidecar.clone();
    format["format_version"] = json!(4);
    assert!(merge_eosmat_audit(&record_file, &format).is_err());
    let mut unknown = sidecar.clone();
    unknown["records"]["missing_record"] = json!({});
    assert!(merge_eosmat_audit(&record_file, &unknown)
        .unwrap_err()
        .to_string()
        .contains("unknown records"));
    let mut extra = sidecar.clone();
    extra["records"]["example_bm3"]["eos"] = json!({});
    assert!(merge_eosmat_audit(&record_file, &extra)
        .unwrap_err()
        .to_string()
        .contains("may contain only"));
    let mut inline = record_file.clone();
    inline["eos_records"][0]["scope_history"] = json!(["Inline note."]);
    assert!(merge_eosmat_audit(&inline, &sidecar)
        .unwrap_err()
        .to_string()
        .contains("both inline and in its audit sidecar"));
    // A contradicting check is detected once the evidence is merged.
    let mut contradiction = record_file;
    contradiction["eos_records"][0]["check"]["reproduction"] = json!("reproduced");
    let merged = merge_eosmat_audit(&contradiction, &sidecar).unwrap();
    assert!(invalid_reason(&merged).contains("check.reproduction"));
}

#[test]
fn path_loads_resolve_sidecars_and_report_missing_ones() {
    let (record_file, sidecar) = split_example();
    let root = std::env::temp_dir().join(format!(
        "peritheos-eosmat-audit-{}-{:?}",
        std::process::id(),
        std::thread::current().id()
    ));
    let materials = root.join("materials");
    let audits = root.join("audits");
    fs::create_dir_all(&materials).unwrap();
    fs::create_dir_all(&audits).unwrap();
    let path = materials.join("example.eosmat");
    fs::write(&path, serde_json::to_string_pretty(&record_file).unwrap()).unwrap();

    let error = load_eosmat(&path).unwrap_err();
    assert_eq!(error.kind(), EosmatErrorKind::MissingAudit);
    assert_eq!(error.code(), "eosmat.audit_missing");
    assert_eq!(error.record_identifier(), Some("example_bm3"));
    assert!(error.to_string().contains("example.audit.json"));
    assert!(std::error::Error::source(&error).is_some());
    // Without the sidecar the record file is still readable as authored.
    assert!(load_eosmat_without_audit(&path).unwrap().eos_records[0].audit_is_external());

    fs::write(audits.join("example.audit.json"), sidecar.to_string()).unwrap();
    let material = load_eosmat(&path).unwrap();
    assert_eq!(
        material.eos_records[0].audit().unwrap()["scope_history"],
        json!(["Fictional scope note."])
    );

    let mut wrong_fragment = record_file;
    wrong_fragment["eos_records"][0]["check"]["audit"] =
        json!("../audits/example.audit.json#another_record");
    fs::write(&path, wrong_fragment.to_string()).unwrap();
    assert!(load_eosmat(&path)
        .unwrap_err()
        .to_string()
        .contains("does not name record"));
    fs::remove_dir_all(root).unwrap();
}

#[test]
fn coded_restrictions_and_holds_form_the_derived_warning() {
    let mut restricted = evidence("parity");
    restricted["reproduction"]["do_not_use"] = json!(true);
    restricted["scientific_validation"]["usage_recommendation"] =
        json!("consult_source_audit_before_quantitative_use");
    let mut document = standalone(&json!({
        "evidence": restricted,
        "catalog_access": "explicit_selection"
    }));
    document["recommendation_hold"] = json!({
        "status": "on_hold",
        "scope": "equilibrium_default",
        "records": ["example_bm3"],
        "reason": "Fictional hold."
    });
    let expected = "Do not use. Consult the source audit before quantitative use. \
                    Default selection on hold pending review. \
                    Explicit selection only; not loaded by default.";
    assert_eq!(
        derive_eosmat_check(&document["eos_records"][0], Some(&document))["warning"],
        expected
    );

    let mut prose = evidence("parity");
    prose["scientific_validation"]["usage_recommendation"] = json!("Use only below 10 GPa.");
    let record = standalone(&json!({ "evidence": prose }));
    assert_eq!(
        derive_eosmat_check(&record["eos_records"][0], None)["warning"],
        Value::Null
    );
    let mut uncoded = evidence("parity");
    uncoded["scientific_validation"]["usage_recommendation"] = json!("diagnostic_only");
    let record = standalone(&json!({ "evidence": uncoded }));
    assert_eq!(
        derive_eosmat_check(&record["eos_records"][0], None)["warning"],
        "Diagnostic only."
    );
}
