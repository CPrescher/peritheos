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
    validate_eosmat_document, Compared, EosmatErrorKind, Material, Reproduction, SourceCheck,
    AUDIT_RECORD_KEYS, RECORD_NOTE_MAX_LENGTH, WARNING_REASON_MAX_LENGTH, WHY_MAX_LENGTH,
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
        let context = format!("{}: {}", path.display(), record.identifier);
        assert_eq!(
            Some(record.reproduction().as_str()),
            check["reproduction"].as_str(),
            "{context}"
        );
        assert_eq!(
            Some(record.compared().as_str()),
            check["compared"].as_str(),
            "{context}"
        );
        assert_eq!(record.why(), check["why"].as_str(), "{context}");
        assert_eq!(record.warning(), check["warning"].as_str(), "{context}");
        assert_eq!(
            record.warning_reason(),
            check["warning_reason"].as_str(),
            "{context}"
        );
        assert_eq!(
            record.source_check().map(SourceCheck::as_str),
            check["source_check"].as_str(),
            "{context}"
        );
        assert_eq!(
            record.source_note(),
            check["source_note"].as_str(),
            "{context}"
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
            assert_eq!(alone.compared(), complete.compared());
            assert_eq!(alone.why(), complete.why());
            assert_eq!(alone.warning(), complete.warning());
            assert_eq!(alone.source_check(), complete.source_check());
            assert_eq!(alone.source_note(), complete.source_note());
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
                record["scientific_validation"]["status"].as_str(),
                source.source_check().map(SourceCheck::status_code),
                "{}: {}",
                path.display(),
                source.identifier
            );
            // Historical readers see the default selection on every record.
            assert_eq!(
                record["default"].as_bool(),
                Some(source.is_default),
                "{}: {}",
                path.display(),
                source.identifier
            );
        }
        // Without check, the v3 headline is derived from the inline evidence.
        for (old, new) in exported.eos_records.iter().zip(&material.eos_records) {
            assert_eq!(old.reproduction(), new.reproduction(), "{}", path.display());
            assert_eq!(old.compared(), new.compared(), "{}", path.display());
            assert_eq!(old.why(), new.why(), "{}", path.display());
            assert_eq!(old.warning(), new.warning(), "{}", path.display());
            assert_eq!(old.source_check(), new.source_check(), "{}", path.display());
            assert_eq!(old.source_note(), new.source_note(), "{}", path.display());
            assert_eq!(old.is_default, new.is_default, "{}", path.display());
        }
        assert_eq!(
            exported.default_record().map(|record| &record.identifier),
            material.default_record().map(|record| &record.identifier),
            "{}",
            path.display()
        );
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
        assert_eq!(old.compared(), new.compared(), "{}", old.identifier);
        assert_eq!(old.why(), new.why(), "{}", old.identifier);
        assert_eq!(old.warning(), new.warning(), "{}", old.identifier);
        assert_eq!(old.is_default, new.is_default, "{}", old.identifier);
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
            "compared": "published_fit",
            "warning": null,
            "source_check": "verified",
            "source_checked_on": "2026-01-02",
            "source_note": "Table 2 values match."
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
        "check": {"reproduction": "close", "source_check": "deferred"}
    }));
    assert!(invalid_reason(&source).contains("check.source_check \"deferred\" contradicts"));

    let mut independent = evidence("parity");
    independent["reproduction"]["assessment"]["scope"] = json!("independent_refit");
    let compared = standalone(&json!({
        "evidence": independent,
        "check": {"compared": "published_fit"}
    }));
    assert!(invalid_reason(&compared).contains("check.compared \"published_fit\" contradicts"));
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
        (
            json!({"source_check": "validated"}),
            "check.source_check must be one of",
        ),
        (
            json!({"compared": "paper"}),
            "check.compared must be one of",
        ),
        (
            json!({"reproduction": "reproduced", "compared": "none"}),
            "cannot accompany",
        ),
        (
            json!({"reproduction": "untested", "compared": "published_fit"}),
            "cannot accompany",
        ),
        (
            json!({"reproduction": "reproduced", "why": "x"}),
            "only given for reproduction close",
        ),
        (
            json!({"why": "x".repeat(WHY_MAX_LENGTH + 1)}),
            "at most 160",
        ),
        (
            json!({"finding": "Table 2 values match."}),
            "check.finding is no longer read; EOSMAT v5 names this field check.source_note",
        ),
        (
            json!({"source": "primary_source_validated"}),
            "check.source is no longer read; EOSMAT v5 names this field check.source_check \
             (primary_source_validated -> verified, pending_primary_source_check -> pending, \
             not_reproduced -> unresolved, deferred -> deferred)",
        ),
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
            json!({"source_note": 1}),
            "check.source_note must be a non-empty string",
        ),
    ] {
        let document = standalone(&json!({ "check": check }));
        assert!(invalid_reason(&document).contains(reason), "{check}");
    }
    let document = standalone(&json!({"check": {
        "reproduction": "untested",
        "compared": "source_equation",
        "warning": null
    }}));
    let record = &material_from_value(document).unwrap().eos_records[0];
    assert_eq!(record.reproduction(), Reproduction::Untested);
    assert_eq!(record.compared(), Compared::SourceEquation);
    assert_eq!(record.warning(), None);
}

fn split_example() -> (Value, Value) {
    let mut audit_evidence = evidence("not_reproduced");
    audit_evidence["scientific_validation"]["usage_recommendation"] = json!("not_recommended");
    let record_file = standalone(&json!({
        "check": {
            "reproduction": "disagrees",
            "compared": "published_fit",
            "warning": "Not recommended for use.",
            "source_check": "verified",
            "source_checked_on": "2026-01-02",
            "source_note": "Table 2 values match.",
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
    assert_eq!(record.compared(), Compared::PublishedFit);
    assert_eq!(record.why(), None);
    assert_eq!(record.warning(), Some("Not recommended for use."));
    assert_eq!(record.source_check(), Some(SourceCheck::Verified));
    assert_eq!(record.source_note(), Some("Table 2 values match."));
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
                    Never chosen automatically; select it explicitly.";
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

fn assessed(status: Option<&str>, scope: &str, rows: Option<Value>) -> Value {
    let mut evidence = evidence("parity");
    match status {
        Some(status) => {
            let assessment = &mut evidence["reproduction"]["assessment"];
            assessment["status"] = json!(status);
            assessment["scope"] = json!(scope);
            if let Some(rows) = rows {
                assessment["parameters"] = rows;
            }
        }
        None => {
            evidence.as_object_mut().unwrap().remove("reproduction");
        }
    }
    standalone(&json!({ "evidence": evidence }))
}

#[allow(clippy::needless_pass_by_value)] // Called with literal json! values.
fn row(parameter: &str, published: Value, refit: f64, flags: &Value) -> Value {
    let mut row = json!({
        "parameter": parameter,
        "published": published,
        "published_error": null,
        "refit": refit,
        "refit_error": null
    });
    row.as_object_mut()
        .unwrap()
        .extend(flags.as_object().unwrap().clone());
    row
}

fn why(document: &Value) -> Option<String> {
    let derived = derive_eosmat_check(&document["eos_records"][0], Some(document));
    let loaded = material_from_value(document.clone()).unwrap();
    let record = &loaded.eos_records[0];
    // Without a stored check the record derives what derive_eosmat_check writes.
    assert_eq!(record.why(), derived["why"].as_str());
    assert_eq!(
        Some(record.compared().as_str()),
        derived["compared"].as_str()
    );
    record.why().map(str::to_owned)
}

#[test]
fn compared_is_a_fixed_mapping_of_status_and_scope() {
    for (status, scope, expected) in [
        (
            Some("parity"),
            "parameter_agreement",
            Compared::PublishedFit,
        ),
        (
            Some("similar"),
            "parameter_agreement",
            Compared::PublishedFit,
        ),
        (
            Some("not_reproduced"),
            "author_workflow",
            Compared::PublishedFit,
        ),
        (
            Some("parity"),
            "independent_refit",
            Compared::PeritheosRefit,
        ),
        (Some("parity"), "source_equation", Compared::SourceEquation),
        (
            Some("source_reconstruction"),
            "source_equation",
            Compared::SourceEquation,
        ),
        (
            Some("not_refittable"),
            "parameter_agreement",
            Compared::Nothing,
        ),
        (Some("not_refittable"), "source_equation", Compared::Nothing),
        (
            Some("not_assessed"),
            "parameter_agreement",
            Compared::Nothing,
        ),
        (None, "", Compared::Nothing),
    ] {
        let document = assessed(status, scope, None);
        let record = &material_from_value(document.clone()).unwrap().eos_records[0];
        assert_eq!(record.compared(), expected, "{status:?} {scope}");
        assert_eq!(record.source_check(), Some(SourceCheck::Verified));
        assert_eq!(record.source_note(), Some("Table 2 values match."));
        // The derived block is a valid stored check for the same evidence.
        let mut checked = document.clone();
        checked["eos_records"][0]["check"] =
            derive_eosmat_check(&document["eos_records"][0], Some(&document));
        validate_eosmat_document(&checked).unwrap();
    }
    for code in [
        "primary_source_validated",
        "pending_primary_source_check",
        "not_reproduced",
        "deferred",
    ] {
        let check = SourceCheck::from_status_code(code).unwrap();
        assert_eq!(check.status_code(), code);
        assert_eq!(SourceCheck::from_check(check.as_str()), Some(check));
    }
}

#[test]
#[allow(clippy::too_many_lines)]
fn why_names_failed_coefficients_from_the_stored_table() {
    // The same cases as tests/test_eosmat_v5_check_fields.py.
    let close = assessed(
        Some("similar"),
        "parameter_agreement",
        Some(json!([
            row(
                "V0",
                json!(162.4),
                162.41,
                &json!({"similar": true, "within_combined_2sigma": true})
            ),
            row(
                "K0",
                json!(256),
                263.7307,
                &json!({"similar": true, "within_combined_2sigma": false})
            ),
        ])),
    );
    assert_eq!(
        why(&close).as_deref(),
        Some("K0: refit 263.7 vs published 256 GPa (+3.0%), outside the combined error.")
    );
    let disagrees = assessed(
        Some("not_reproduced"),
        "parameter_agreement",
        Some(json!([{
            "parameter": "rt_eos.K0_prime",
            "published": 4.91,
            "published_sigma": 0.12,
            "refit": 4.183_456,
            "refit_sigma": null,
            "similar": false,
            "within_combined_2sigma": false,
            "within_reported_error": false
        }])),
    );
    assert_eq!(
        why(&disagrees).as_deref(),
        Some(
            "K0': refit 4.183 vs published 4.91 (-14.8%), outside the combined error, \
             beyond the agreement limit and outside the published error."
        )
    );
    let untested = assessed(
        Some("similar"),
        "parameter_agreement",
        Some(json!([
            row(
                "K0",
                json!(160),
                161.0,
                &json!({"similar": true, "within_combined_2sigma": null})
            ),
            row(
                "K0_prime",
                json!(4),
                4.1,
                &json!({"similar": true, "within_combined_2sigma": null})
            ),
        ])),
    );
    assert_eq!(
        why(&untested).as_deref(),
        Some(
            "K0 and K0' have no error to test against the combined-error rule; \
             all within the agreement limit."
        )
    );
    let single = assessed(
        Some("similar"),
        "parameter_agreement",
        Some(json!([row(
            "V0",
            json!(74.7),
            74.71,
            &json!({"similar": true, "within_combined_2sigma": null})
        )])),
    );
    assert_eq!(
        why(&single).as_deref(),
        Some(
            "V0 has no error to test against the combined-error rule; \
             all within the agreement limit."
        )
    );
    // Published values print as stored, without a trailing ".0".
    let integral = assessed(
        Some("similar"),
        "parameter_agreement",
        Some(json!([row(
            "theta0",
            json!(600.0),
            766.2634,
            &json!({"similar": false})
        )])),
    );
    assert_eq!(
        why(&integral).as_deref(),
        Some("theta0: refit 766.3 vs published 600 K (+27.7%), beyond the agreement limit.")
    );
    // Exponent notation follows the stored spelling's significant digits.
    let small = assessed(
        Some("similar"),
        "parameter_agreement",
        Some(json!([row(
            "alpha1",
            json!(8.7e-9),
            5.6123e-9,
            &json!({"similar": false})
        )])),
    );
    assert_eq!(
        why(&small).as_deref(),
        Some("alpha1: refit 5.61e-09 vs published 8.7e-09 (-35.5%), beyond the agreement limit.")
    );
    // More than three failures: the three largest relative differences are
    // named and the rest counted.
    let four = assessed(
        Some("not_reproduced"),
        "parameter_agreement",
        Some(json!([
            row("V0", json!(74.7), 75.0, &json!({"similar": false})),
            row("K0", json!(160), 200.0, &json!({"similar": false})),
            row("K0_prime", json!(4), 2.0, &json!({"similar": false})),
            row("gamma0", json!(1.5), 1.6, &json!({"similar": false})),
        ])),
    );
    assert_eq!(
        why(&four).as_deref(),
        Some("K0' (-50.0%), K0 (+25.0%), gamma0 (+6.7%) and 1 more: beyond the agreement limit.")
    );
    // Nothing failed and nothing untested, or a parity headline: no statement.
    let passing = assessed(
        Some("not_reproduced"),
        "parameter_agreement",
        Some(json!([row(
            "K0",
            json!(160),
            161.0,
            &json!({"similar": true, "within_combined_2sigma": true})
        )])),
    );
    assert_eq!(why(&passing), None);
    let parity = assessed(
        Some("parity"),
        "parameter_agreement",
        Some(
            close["eos_records"][0]["evidence"]["reproduction"]["assessment"]["parameters"].clone(),
        ),
    );
    assert_eq!(why(&parity), None);
    assert_eq!(
        why(&assessed(Some("similar"), "parameter_agreement", None)),
        None
    );
    let many: Vec<Value> = (0..12)
        .map(|index| {
            row(
                &format!("coefficient_{index}"),
                json!(1.25),
                1.5 + f64::from(index),
                &json!({"similar": false}),
            )
        })
        .collect();
    let many = assessed(
        Some("not_reproduced"),
        "parameter_agreement",
        Some(json!(many)),
    );
    assert!(why(&many).is_none_or(|text| text.chars().count() <= WHY_MAX_LENGTH));
}

fn two_records(map: Option<Value>, first: &Value, second: &Value) -> Value {
    let mut document = standalone(&json!({}));
    let base = document["eos_records"][0].clone();
    let mut a = base.clone();
    a.as_object_mut()
        .unwrap()
        .extend(first.as_object().unwrap().clone());
    let mut b = base;
    b["identifier"] = json!("example_vinet");
    b["model"] = json!("Vinet");
    b.as_object_mut()
        .unwrap()
        .extend(second.as_object().unwrap().clone());
    document["eos_records"] = json!([a, b]);
    if let Some(map) = map {
        document["default_records"] = map;
    }
    document
}

#[test]
fn default_records_map_selects_exactly_like_record_flags() {
    let mapped = two_records(
        Some(json!({"equilibrium": "example_vinet"})),
        &json!({}),
        &json!({}),
    );
    let flagged = two_records(None, &json!({}), &json!({"default": true}));
    let by_map = material_from_value(mapped.clone()).unwrap();
    let by_flag = material_from_value(flagged).unwrap();
    assert_eq!(by_map.default_record().unwrap().identifier, "example_vinet");
    assert_eq!(
        by_map
            .eos_records
            .iter()
            .map(|r| r.is_default)
            .collect::<Vec<_>>(),
        by_flag
            .eos_records
            .iter()
            .map(|r| r.is_default)
            .collect::<Vec<_>>()
    );
    // The retained document keeps the map as authored.
    assert_eq!(
        by_map.document["default_records"],
        json!({"equilibrium": "example_vinet"})
    );
    // A v3 export writes an explicit flag on every record instead.
    let legacy = export_eosmat_v3(&mapped).unwrap();
    assert!(legacy.get("default_records").is_none());
    assert_eq!(legacy["eos_records"][0]["default"], false);
    assert!(legacy["eos_records"][0].get("default_for").is_none());
    assert_eq!(legacy["eos_records"][1]["default"], true);
    assert_eq!(legacy["eos_records"][1]["default_for"], "equilibrium");
    // Both forms in one file are accepted only when they agree.
    validate_eosmat_document(&two_records(
        Some(json!({"equilibrium": "example_vinet"})),
        &json!({}),
        &json!({"default_for": "equilibrium"}),
    ))
    .unwrap();
    let map = Some(json!({"equilibrium": "example_vinet"}));
    for (document, message) in [
        (
            two_records(map.clone(), &json!({"default": true}), &json!({})),
            "flagged as the default equilibrium record, but default_records names 'example_vinet'",
        ),
        (
            two_records(map.clone(), &json!({}), &json!({"default": false})),
            "flagged default: false",
        ),
        (
            two_records(
                Some(json!({"equilibrium": "missing"})),
                &json!({}),
                &json!({}),
            ),
            "unknown record",
        ),
        (
            two_records(
                Some(json!({"hugoniot": "example_vinet"})),
                &json!({}),
                &json!({}),
            ),
            "which is not a hugoniot record",
        ),
        (
            two_records(
                Some(json!({"best": "example_vinet"})),
                &json!({}),
                &json!({}),
            ),
            "keys must be",
        ),
        (
            two_records(Some(json!({})), &json!({}), &json!({})),
            "non-empty object",
        ),
    ] {
        let reason = invalid_reason(&document);
        assert!(reason.contains(message), "{reason}");
        assert!(material_from_value(document).is_err());
    }
    let mut v4 = mapped;
    v4["format_version"] = json!(4);
    assert!(invalid_reason(&v4).contains("requires EOSMAT format_version 5"));
}

#[test]
fn bundled_default_selection_uses_the_material_map() {
    let Some(directory) = catalog_directory() else {
        return;
    };
    let mut mapped = 0;
    for path in bundled_paths(&directory) {
        let stored: Value = serde_json::from_slice(&fs::read(&path).unwrap()).unwrap();
        let material = load_eosmat_without_audit(&path).unwrap();
        for record in stored["eos_records"].as_array().unwrap() {
            assert!(record.get("default").is_none(), "{}", path.display());
            assert!(record.get("default_for").is_none(), "{}", path.display());
        }
        let Some(map) = stored.get("default_records") else {
            assert!(
                material.eos_records.iter().all(|record| !record.is_default),
                "{}",
                path.display()
            );
            continue;
        };
        mapped += 1;
        let named: Vec<&str> = map
            .as_object()
            .unwrap()
            .values()
            .map(|value| value.as_str().unwrap())
            .collect();
        for record in &material.eos_records {
            assert_eq!(
                record.is_default,
                named.contains(&record.identifier.as_str()),
                "{}: {}",
                path.display(),
                record.identifier
            );
        }
        if let Some(identifier) = map.get("equilibrium").and_then(Value::as_str) {
            if let Some(default) = material.default_equilibrium_record() {
                assert_eq!(default.identifier, identifier, "{}", path.display());
            }
        }
        if let Some(identifier) = map.get("hugoniot").and_then(Value::as_str) {
            assert_eq!(
                material
                    .default_hugoniot_record()
                    .unwrap()
                    .record
                    .identifier,
                identifier,
                "{}",
                path.display()
            );
        }
    }
    assert!(mapped > 50, "unexpectedly few default_records maps");
}

#[test]
fn superseded_record_text_is_an_audit_block_read_like_the_others() {
    let block = json!({
        "label": "Example (2020), BM3 (diagnostic only)",
        "notes": "Original wording, kept verbatim.",
        "validity_notes": ["Original range note."],
        "thermal_notes": null,
        "superseded_on": "2026-11-16"
    });
    assert!(AUDIT_RECORD_KEYS.contains(&"superseded_record_text"));
    let (record_file, mut sidecar) = split_example();
    sidecar["records"]["example_bm3"]["superseded_record_text"] = block.clone();
    let merged =
        load_eosmat_str_with_audit(&record_file.to_string(), &sidecar.to_string()).unwrap();
    let audit = merged.eos_records[0].audit().unwrap();
    assert_eq!(audit["superseded_record_text"], block);
    // Compatibility exports carry it inline, like the other audit blocks.
    let legacy = export_eosmat_v3(&merged.document).unwrap();
    assert_eq!(legacy["eos_records"][0]["superseded_record_text"], block);
    // A standalone file may keep it inline.
    let inline =
        material_from_value(standalone(&json!({"superseded_record_text": block}))).unwrap();
    assert_eq!(
        inline.eos_records[0].audit().unwrap()["superseded_record_text"],
        block
    );
    let mut both = record_file.clone();
    both["eos_records"][0]["superseded_record_text"] = block.clone();
    assert!(merge_eosmat_audit(&both, &sidecar)
        .unwrap_err()
        .to_string()
        .contains("both inline and in its audit sidecar"));
    for (invalid, message) in [
        (json!({"notes": "x"}), "superseded_on must be an ISO date"),
        (
            json!({"notes": "x", "superseded_on": "2026-02-30"}),
            "superseded_on must be an ISO date",
        ),
        (
            json!({"note": "x", "superseded_on": "2026-11-16"}),
            "unknown fields",
        ),
        (
            json!({"superseded_on": "2026-11-16"}),
            "names no superseded field",
        ),
        (
            json!({"validity_notes": "x", "superseded_on": "2026-11-16"}),
            "validity_notes must be an array of strings or null",
        ),
        (
            json!({"label": 3, "superseded_on": "2026-11-16"}),
            "label must be a string or null",
        ),
    ] {
        let mut wrong = sidecar.clone();
        wrong["records"]["example_bm3"]["superseded_record_text"] = invalid.clone();
        let error = merge_eosmat_audit(&record_file, &wrong).unwrap_err();
        assert!(error.to_string().contains(message), "{error}");
        let error = material_from_value(standalone(&json!({"superseded_record_text": invalid})))
            .unwrap_err();
        assert!(error.to_string().contains(message), "{error}");
    }
}

#[test]
fn bundled_superseded_record_text_matches_the_sidecar() {
    // Python reads the same sidecars (tests/test_record_text.py); both must
    // expose each stored block unchanged through the record's audit.
    let Some(directory) = catalog_directory() else {
        return;
    };
    let mut found = 0;
    for path in bundled_paths(&directory) {
        let material = load_eosmat(&path).unwrap();
        let sidecar_path = directory
            .join("../audits")
            .join(format!("{}.audit.json", material.identifier));
        let sidecar: Value =
            serde_json::from_str(&fs::read_to_string(sidecar_path).unwrap()).unwrap();
        for record in &material.eos_records {
            let stored = sidecar["records"]
                .get(&record.identifier)
                .and_then(|entry| entry.get("superseded_record_text"));
            let audit = record.audit().unwrap();
            assert_eq!(
                audit.get("superseded_record_text"),
                stored,
                "{}",
                record.identifier
            );
            found += usize::from(stored.is_some());
        }
    }
    assert!(found > 0, "no bundled record keeps superseded text");
    // The packaged Rust fixture reads its own copy of the Pd sidecar.
    let fixture = load_eosmat(crate_directory().join("tests/data/palladium.eosmat")).unwrap();
    let record = fixture
        .eos_records
        .iter()
        .find(|record| record.identifier == "palladium_fedotenko_2020_bm2_2")
        .unwrap();
    assert_eq!(
        record.audit().unwrap()["superseded_record_text"]["label"],
        "Fedotenko et al. (2020), constrained BM2 (poor fit; diagnostic only)"
    );
    assert_eq!(record.label, "Fedotenko et al. (2020), constrained BM2");
}

/// A restricted record whose evidence carries both curated check texts.
fn curated(note: Option<&str>, reason: Option<&str>, restricted: bool) -> Value {
    let mut evidence = evidence("parity");
    if restricted {
        evidence["scientific_validation"]["usage_recommendation"] = json!("qualified_use");
    }
    if let Some(note) = note {
        evidence["scientific_validation"]["primary_source_check"]["record_note"] = json!(note);
    }
    if let Some(reason) = reason {
        evidence["scientific_validation"]["warning_reason"] = json!(reason);
    }
    standalone(&json!({ "evidence": evidence }))
}

#[test]
fn record_note_and_warning_reason_are_copied_into_check() {
    let note = "Table 2 V0 and K0 match; V0 converted from cm3/mol.";
    let reason = "The fit omits the two highest-pressure points.";
    let document = curated(Some(note), Some(reason), true);
    let derived = derive_eosmat_check(&document["eos_records"][0], Some(&document));
    assert_eq!(derived["source_note"], note);
    assert_eq!(derived["warning"], "Qualified use only; see the audit.");
    assert_eq!(derived["warning_reason"], reason);
    let record = material_from_value(document.clone()).unwrap().eos_records[0].clone();
    assert_eq!(record.source_note(), Some(note));
    assert_eq!(record.warning_reason(), Some(reason));

    // Stored values: the check block Python writes for this record (key order
    // as formatted) reads back unchanged and agrees with the derivation.
    let mut stored = document.clone();
    stored["eos_records"][0]["check"] = json!({
        "reproduction": "reproduced",
        "compared": "published_fit",
        "warning": "Qualified use only; see the audit.",
        "warning_reason": reason,
        "source_check": "verified",
        "source_checked_on": "2026-01-02",
        "source_note": note
    });
    assert_eq!(stored["eos_records"][0]["check"], derived);
    let record = material_from_value(stored).unwrap().eos_records[0].clone();
    assert_eq!(record.warning_reason(), Some(reason));
    assert_eq!(record.source_note(), Some(note));

    // Absent: the verbatim finding, and no reason.
    let plain = curated(None, None, true);
    let derived = derive_eosmat_check(&plain["eos_records"][0], None);
    assert_eq!(derived["source_note"], "Table 2 values match.");
    assert!(derived.get("warning_reason").is_none());
    let record = material_from_value(plain).unwrap().eos_records[0].clone();
    assert_eq!(record.source_note(), Some("Table 2 values match."));
    assert_eq!(record.warning_reason(), None);
}

#[test]
fn curated_check_texts_are_validated() {
    let long_note = "x".repeat(RECORD_NOTE_MAX_LENGTH + 1);
    let long_reason = "x".repeat(WARNING_REASON_MAX_LENGTH + 1);
    let cases = [
        (
            curated(Some(&long_note), None, true),
            "primary_source_check.record_note must be at most 300 characters",
        ),
        (
            curated(None, Some(&long_reason), true),
            "warning_reason must be at most 140 characters",
        ),
        (
            curated(Some(" "), None, true),
            "record_note must be a non-empty string",
        ),
        (
            curated(None, Some(""), true),
            "warning_reason must be a non-empty string",
        ),
        (
            curated(None, Some("Fictional reason."), false),
            "warning_reason is given but the record has no warning",
        ),
    ];
    for (document, message) in cases {
        let reason = invalid_reason(&document);
        assert!(reason.contains(message), "{reason}");
    }
    // Within the limits is accepted.
    let note = "x".repeat(RECORD_NOTE_MAX_LENGTH);
    let reason = "x".repeat(WARNING_REASON_MAX_LENGTH);
    validate_eosmat_document(&curated(Some(&note), Some(&reason), true)).unwrap();

    // The check field itself: a reason needs a warning and obeys the limit.
    for (check, message) in [
        (
            json!({"warning": null, "warning_reason": "Fictional reason."}),
            "check.warning_reason is only given with a warning",
        ),
        (
            json!({"warning_reason": "Fictional reason."}),
            "check.warning_reason is only given with a warning",
        ),
        (
            json!({"warning": "Do not use.", "warning_reason": long_reason}),
            "check.warning_reason must be at most 140 characters",
        ),
        (
            json!({"warning": "Do not use.", "warning_reason": ""}),
            "check.warning_reason must be a non-empty string",
        ),
    ] {
        let reason = invalid_reason(&standalone(&json!({ "check": check })));
        assert!(reason.contains(message), "{reason}");
    }
    let record = material_from_value(standalone(&json!({
        "check": {"warning": "Do not use.", "warning_reason": "Fictional reason."}
    })))
    .unwrap()
    .eos_records[0]
        .clone();
    assert_eq!(record.warning(), Some("Do not use."));
    assert_eq!(record.warning_reason(), Some("Fictional reason."));
}

/// Mirrors `test_stored_copies_must_equal_the_inline_evidence` in Python.
#[test]
fn stored_check_texts_must_equal_the_evidence() {
    let note = "Table 2 V0 and K0 match; V0 converted from cm3/mol.";
    let reason = "The fit omits the two highest-pressure points.";
    let with_check = |document: &Value, key: &str, value: Option<&str>| {
        let mut document = document.clone();
        let mut check = derive_eosmat_check(&document["eos_records"][0], Some(&document));
        if let Some(value) = value {
            check[key] = json!(value);
        }
        document["eos_records"][0]["check"] = check;
        document
    };
    for document in [
        curated(Some(note), Some(reason), true),
        curated(None, None, true),
    ] {
        material_from_value(with_check(&document, "source_note", None)).unwrap();
    }
    for (document, key, value) in [
        (
            curated(Some(note), Some(reason), true),
            "source_note",
            "Table 2 values match.",
        ),
        (curated(None, Some(reason), true), "source_note", note),
        (
            curated(Some(note), Some(reason), true),
            "warning_reason",
            "Other reason.",
        ),
        (curated(Some(note), None, true), "warning_reason", reason),
    ] {
        let document = with_check(&document, key, Some(value));
        let message = invalid_reason(&document);
        assert!(
            message.contains(&format!("check.{key} differs from the evidence")),
            "{message}"
        );
    }
    // A standalone file without evidence keeps its own texts.
    let mut own = curated(None, None, true);
    own["eos_records"][0]
        .as_object_mut()
        .unwrap()
        .remove("evidence");
    own["eos_records"][0]["check"] = json!({
        "warning": "Do not use.",
        "warning_reason": "My own reason.",
        "source_note": "My own note."
    });
    let record = material_from_value(own).unwrap().eos_records[0].clone();
    assert_eq!(record.source_note(), Some("My own note."));
}

/// Blank means Unicode `White_Space` only, as in Python; limits count characters.
#[test]
fn curated_texts_trim_white_space_and_count_characters() {
    for blank in [" ", "\t\n", "\u{a0}", "\u{2003}\u{3000}"] {
        let mut plain = curated(None, None, true);
        plain["eos_records"][0]
            .as_object_mut()
            .unwrap()
            .remove("evidence");
        for check in [
            json!({"warning": "Do not use.", "warning_reason": blank}),
            json!({"source_note": blank}),
        ] {
            plain["eos_records"][0]["check"] = check;
            let message = invalid_reason(&plain);
            assert!(message.contains("must be a non-empty string"), "{message}");
        }
        for document in [
            curated(Some(blank), None, true),
            curated(None, Some(blank), true),
        ] {
            let message = invalid_reason(&document);
            assert!(message.contains("must be a non-empty string"), "{message}");
        }
    }
    // U+001F is not White_Space: kept, as by Python.
    validate_eosmat_document(&curated(Some("\u{1f}"), None, true)).unwrap();

    // Two-byte characters: the limits are characters, not bytes.
    let note = "é".repeat(RECORD_NOTE_MAX_LENGTH);
    let reason = "Å".repeat(WARNING_REASON_MAX_LENGTH);
    assert_eq!(note.len(), 2 * RECORD_NOTE_MAX_LENGTH);
    let document = curated(Some(&note), Some(&reason), true);
    validate_eosmat_document(&document).unwrap();
    let mut stored = document.clone();
    stored["eos_records"][0]["check"] =
        derive_eosmat_check(&document["eos_records"][0], Some(&document));
    let record = material_from_value(stored).unwrap().eos_records[0].clone();
    assert_eq!(record.warning_reason(), Some(reason.as_str()));
    let longer = format!("{note}é");
    assert!(invalid_reason(&curated(Some(&longer), None, true))
        .contains("record_note must be at most 300 characters"));
    let longer = format!("{reason}Å");
    assert!(invalid_reason(&curated(None, Some(&longer), true))
        .contains("warning_reason must be at most 140 characters"));
    let mut plain = curated(None, None, true);
    plain["eos_records"][0]
        .as_object_mut()
        .unwrap()
        .remove("evidence");
    plain["eos_records"][0]["check"] = json!({"warning": "Do not use.", "warning_reason": reason});
    validate_eosmat_document(&plain).unwrap();
    plain["eos_records"][0]["check"]["warning_reason"] = json!(longer);
    assert!(invalid_reason(&plain).contains("check.warning_reason must be at most 140"));
}
