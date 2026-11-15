//! EOSMAT v5 record `check` headlines and per-material audit sidecars.
//!
//! This mirrors `peritheos/_eosmat_audit.py`. A bundled material is stored as
//! a record file (`<id>.eosmat`, whose records end with a short `check` block)
//! and an audit sidecar (`<id>.audit.json`) holding each record's evidence,
//! detailed `pressure_calibration` and other audit-only blocks. Merging the
//! sidecar restores the complete inline v5 document. Nothing here reclassifies
//! a record: `check.reproduction` is a fixed mapping of the stored assessment
//! status and `check.warning` is built only from stored coded restrictions.

use std::collections::{BTreeMap, BTreeSet};
use std::fmt::{self, Display, Formatter};
use std::path::Path;

use serde_json::{Map, Value};

use super::{invalid_document, EosmatError};

/// Format discriminator of an EOSMAT v5 audit sidecar.
pub const EOSMAT_AUDIT_FORMAT: &str = "peritheos.eosmat-audit";

/// Record keys that describe how a record was checked, not what it is.
///
/// Bundled record files store these in the audit sidecar; every other record
/// key stays in the record file.
pub const AUDIT_RECORD_KEYS: &[&str] = &[
    "evidence",
    "pressure_calibration",
    "fit_input_inventory",
    "network_reconstruction_audit",
    "reproduction_audit",
    "numerical_reproduction",
    "validation_selection",
    "migration_corrections",
    "scope_history",
    "author_input_followup",
];

const SOURCE_STATUSES: &[&str] = &[
    "primary_source_validated",
    "pending_primary_source_check",
    "not_reproduced",
    "deferred",
];
const CHECK_KEYS: &[&str] = &[
    "reproduction",
    "warning",
    "finding",
    "source",
    "source_checked_on",
    "audit",
];

const REPORT_LEDGER: &str = "docs/data/primary-eos-refits.json";

const DO_NOT_USE_WARNING: &str = "Do not use.";
const HOLD_WARNING: &str = "Default selection on hold pending review.";
const EXPLICIT_SELECTION_WARNING: &str = "Explicit selection only; not loaded by default.";

/// Headline outcome of a record's reconstruction assessment (`check.reproduction`).
///
/// This is a fixed mapping of `evidence.reproduction.assessment.status`:
/// `parity` is [`Reproduction::Reproduced`], `similar` is
/// [`Reproduction::Close`], `not_reproduced` is [`Reproduction::Disagrees`],
/// and every other or absent status is [`Reproduction::Untested`]. It is not
/// a physical validation or a recommendation for use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Reproduction {
    /// A refit of the source data reproduced the published values (`parity`).
    Reproduced,
    /// A refit agreed within the stated similarity scope (`similar`).
    Close,
    /// A refit did not reproduce the published values (`not_reproduced`).
    Disagrees,
    /// No reproduction was assessed, or it could not be attempted.
    Untested,
}

impl Reproduction {
    /// The EOSMAT `check.reproduction` spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Reproduced => "reproduced",
            Self::Close => "close",
            Self::Disagrees => "disagrees",
            Self::Untested => "untested",
        }
    }

    /// Parse a `check.reproduction` value.
    #[must_use]
    pub fn from_check(value: &str) -> Option<Self> {
        match value {
            "reproduced" => Some(Self::Reproduced),
            "close" => Some(Self::Close),
            "disagrees" => Some(Self::Disagrees),
            "untested" => Some(Self::Untested),
            _ => None,
        }
    }

    /// Map a stored `evidence.reproduction.assessment.status`.
    ///
    /// Unknown statuses are [`Reproduction::Untested`]; nothing is inferred.
    #[must_use]
    pub fn from_assessment_status(status: &str) -> Self {
        match status {
            "parity" => Self::Reproduced,
            "similar" => Self::Close,
            "not_reproduced" => Self::Disagrees,
            _ => Self::Untested,
        }
    }
}

impl Display for Reproduction {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// Evidence object of a native record, or the flattened v3 record itself.
fn evidence(record: &Value) -> &Value {
    match record.get("evidence") {
        Some(value) if value.is_object() => value,
        _ => record,
    }
}

fn section<'a>(value: &'a Value, key: &str) -> Option<&'a Map<String, Value>> {
    value.get(key).and_then(Value::as_object)
}

fn assessment(record: &Value) -> Option<&Map<String, Value>> {
    evidence(record)
        .get("reproduction")
        .and_then(|reproduction| section(reproduction, "assessment"))
}

fn assessment_status(record: &Value) -> Option<String> {
    // Python maps `str(status)`; only string statuses can match a known value.
    assessment(record)
        .and_then(|assessment| assessment.get("status"))
        .map(|status| match status {
            Value::String(text) => text.clone(),
            other => other.to_string(),
        })
}

/// Map the stored assessment status of a record; anything else is untested.
pub(super) fn reproduction_headline(record: &Value) -> Reproduction {
    assessment_status(record).map_or(Reproduction::Untested, |status| {
        Reproduction::from_assessment_status(&status)
    })
}

fn hold_covers(hold: &Value, identifier: &str) -> bool {
    if hold.get("status").and_then(Value::as_str) != Some("on_hold") {
        return false;
    }
    match hold.get("records") {
        Some(Value::Array(records)) => records.iter().any(|record| record == identifier),
        _ => hold.get("record").and_then(Value::as_str) == Some(identifier),
    }
}

fn record_identifier(record: &Value) -> String {
    match record.get("identifier") {
        Some(Value::String(text)) => text.clone(),
        Some(other) => other.to_string(),
        None => "None".to_owned(),
    }
}

/// Whether a usage recommendation is a coded lower-snake-case value.
fn is_coded(recommendation: &str) -> bool {
    let letters: String = recommendation.chars().filter(|c| *c != '_').collect();
    !letters.is_empty()
        && letters.chars().all(char::is_alphabetic)
        && recommendation.chars().any(char::is_lowercase)
        && !recommendation.chars().any(char::is_uppercase)
}

fn recommendation_warning(recommendation: &str) -> String {
    match recommendation {
        "not_recommended" => "Not recommended for use.".to_owned(),
        "not_recommended_for_quantitative_use" => {
            "Not recommended for quantitative use.".to_owned()
        }
        "consult_source_audit_before_quantitative_use" => {
            "Consult the source audit before quantitative use.".to_owned()
        }
        "qualified_use" => "Qualified use only; see the audit.".to_owned(),
        "on_hold_not_recommended_pending_source_mapping_review" => {
            "On hold; not recommended pending source mapping review.".to_owned()
        }
        other => {
            let words = other.replace('_', " ");
            let mut characters = words.chars();
            let first = characters
                .next()
                .map(|c| c.to_uppercase().collect::<String>())
                .unwrap_or_default();
            format!("{first}{}.", characters.as_str())
        }
    }
}

/// Join the fixed phrases for stored, coded usage restrictions of a record.
///
/// Sources, in order: `evidence.reproduction.do_not_use` (true), a coded
/// `scientific_validation.usage_recommendation`, a material-level
/// `recommendation_hold` naming the record, and `catalog_access:
/// explicit_selection`. Prose fields are never interpreted.
pub(super) fn usage_warning(record: &Value, material: Option<&Value>) -> Option<String> {
    let evidence = evidence(record);
    let mut phrases: Vec<String> = Vec::new();
    if section(evidence, "reproduction")
        .and_then(|reproduction| reproduction.get("do_not_use"))
        .and_then(Value::as_bool)
        == Some(true)
    {
        phrases.push(DO_NOT_USE_WARNING.to_owned());
    }
    let recommendation = section(evidence, "scientific_validation")
        .and_then(|validation| validation.get("usage_recommendation"))
        .and_then(Value::as_str)
        .unwrap_or("");
    let coded = is_coded(recommendation);
    if coded {
        phrases.push(recommendation_warning(recommendation));
    }
    let on_hold_code = coded && recommendation.starts_with("on_hold");
    if let Some(material) = material {
        if !on_hold_code
            && material
                .get("recommendation_hold")
                .is_some_and(|hold| hold_covers(hold, &record_identifier(record)))
        {
            phrases.push(HOLD_WARNING.to_owned());
        }
    }
    if record.get("catalog_access").and_then(Value::as_str) == Some("explicit_selection") {
        phrases.push(EXPLICIT_SELECTION_WARNING.to_owned());
    }
    if phrases.is_empty() {
        None
    } else {
        Some(phrases.join(" "))
    }
}

/// Derive a record's `check` block from its inline evidence.
pub(super) fn check_block(record: &Value, material: Option<&Value>, audit: Option<&str>) -> Value {
    let validation = section(evidence(record), "scientific_validation");
    let mut result = Map::new();
    result.insert(
        "reproduction".into(),
        Value::String(reproduction_headline(record).as_str().to_owned()),
    );
    result.insert(
        "warning".into(),
        usage_warning(record, material).map_or(Value::Null, Value::String),
    );
    if let Some(validation) = validation {
        if let Some(finding) = validation
            .get("primary_source_check")
            .and_then(Value::as_object)
            .and_then(|check| check.get("finding"))
            .and_then(Value::as_str)
            .filter(|finding| !finding.is_empty())
        {
            result.insert("finding".into(), Value::String(finding.to_owned()));
        }
        if let Some(status) = validation.get("status") {
            result.insert("source".into(), status.clone());
        }
        if let Some(date) = validation.get("audit_date") {
            result.insert("source_checked_on".into(), date.clone());
        }
    }
    if let Some(audit) = audit {
        result.insert("audit".into(), Value::String(audit.to_owned()));
    }
    Value::Object(result)
}

/// Effective reproduction headline and warning of a record.
///
/// A stored `check` value wins; otherwise the value is derived from inline
/// evidence. No assessment means [`Reproduction::Untested`].
pub(super) fn headline(record: &Value, material: Option<&Value>) -> (Reproduction, Option<String>) {
    let check = record.get("check").and_then(Value::as_object);
    let reproduction = check
        .and_then(|check| check.get("reproduction"))
        .and_then(Value::as_str)
        .filter(|value| !value.is_empty())
        .and_then(Reproduction::from_check)
        .unwrap_or_else(|| reproduction_headline(record));
    let warning = match check.and_then(|check| check.get("warning")) {
        Some(value) => value.as_str().map(str::to_owned),
        None => usage_warning(record, material),
    };
    (reproduction, warning)
}

/// Validate the structure of an optional record `check` block.
pub(super) fn validate_check(check: &Value, location: &str) -> Result<(), EosmatError> {
    let check = check
        .as_object()
        .ok_or_else(|| invalid_document(format!("{location}.check must be an object")))?;
    let unknown: BTreeSet<&str> = check
        .keys()
        .map(String::as_str)
        .filter(|key| !CHECK_KEYS.contains(key))
        .collect();
    if !unknown.is_empty() {
        return Err(invalid_document(format!(
            "{location}.check has unknown fields {unknown:?}"
        )));
    }
    if let Some(value) = check.get("reproduction") {
        if value.as_str().and_then(Reproduction::from_check).is_none() {
            return Err(invalid_document(format!(
                "{location}.check.reproduction must be one of reproduced, close, disagrees, untested"
            )));
        }
    }
    if let Some(value) = check.get("source") {
        if !value.as_str().is_some_and(|v| SOURCE_STATUSES.contains(&v)) {
            return Err(invalid_document(format!(
                "{location}.check.source is invalid"
            )));
        }
    }
    if let Some(value) = check.get("warning") {
        if !value.is_null() && value.as_str().is_none_or(|text| text.trim().is_empty()) {
            return Err(invalid_document(format!(
                "{location}.check.warning must be null or text"
            )));
        }
    }
    for key in ["finding", "source_checked_on", "audit"] {
        if let Some(value) = check.get(key) {
            if value.as_str().is_none_or(str::is_empty) {
                return Err(invalid_document(format!(
                    "{location}.check.{key} must be a non-empty string"
                )));
            }
        }
    }
    Ok(())
}

/// A `check` headline must agree with the inline evidence it summarizes.
pub(super) fn validate_check_consistency(
    record: &Value,
    location: &str,
) -> Result<(), EosmatError> {
    let Some(check) = record.get("check").and_then(Value::as_object) else {
        return Ok(());
    };
    if let (Some(stored), Some(status)) = (check.get("reproduction"), assessment_status(record)) {
        let expected = match status.as_str() {
            "parity"
            | "similar"
            | "not_reproduced"
            | "not_refittable"
            | "source_reconstruction"
            | "not_assessed" => Some(Reproduction::from_assessment_status(&status).as_str()),
            _ => None,
        };
        if stored.as_str() != expected {
            return Err(invalid_document(format!(
                "{location}.check.reproduction {stored} contradicts assessment status {status:?}"
            )));
        }
    }
    if let (Some(source), Some(status)) = (
        check.get("source"),
        section(evidence(record), "scientific_validation")
            .and_then(|validation| validation.get("status")),
    ) {
        if source != status {
            return Err(invalid_document(format!(
                "{location}.check.source contradicts scientific_validation.status"
            )));
        }
    }
    Ok(())
}

/// Whether a record names a sidecar entry and carries no inline audit.
pub(super) fn is_external(record: &Value) -> bool {
    record
        .get("check")
        .and_then(|check| check.get("audit"))
        .is_some_and(Value::is_string)
        && !AUDIT_RECORD_KEYS
            .iter()
            .any(|key| record.get(*key).is_some())
}

/// Reconstruction-report links implied by the material and record identifiers.
fn derived_report_links(material: &str, record: &str) -> [(&'static str, String); 3] {
    [
        (
            "report_url",
            format!(
                "https://peritheos.readthedocs.io/en/latest/reconstruction-reports/{material}/#record-{record}"
            ),
        ),
        (
            "report_path",
            format!("docs/reconstruction-reports/{material}.md#record-{record}"),
        ),
        ("report", REPORT_LEDGER.to_owned()),
    ]
}

/// Regenerate the identifier-derived report links of a bundled entry.
fn expand_entry(entry: &Map<String, Value>, material: &str, record: &str) -> Map<String, Value> {
    let mut result = entry.clone();
    if let Some(assessment) = result
        .get_mut("evidence")
        .and_then(Value::as_object_mut)
        .and_then(|evidence| evidence.get_mut("reproduction"))
        .and_then(Value::as_object_mut)
        .and_then(|reproduction| reproduction.get_mut("assessment"))
        .and_then(Value::as_object_mut)
        .filter(|assessment| !assessment.is_empty())
    {
        for (key, value) in derived_report_links(material, record) {
            assessment
                .entry(key.to_owned())
                .or_insert(Value::String(value));
        }
    }
    result
}

fn identifier_text(value: Option<&Value>) -> String {
    match value {
        Some(Value::String(text)) => text.clone(),
        Some(other) => other.to_string(),
        None => "None".to_owned(),
    }
}

/// Validate the envelope and entry keys of an audit sidecar.
pub(super) fn validate_sidecar(sidecar: &Value, material: &str) -> Result<(), EosmatError> {
    let object = sidecar
        .as_object()
        .ok_or_else(|| invalid_document("an EOSMAT audit sidecar must be a JSON object"))?;
    if object.get("format").and_then(Value::as_str) != Some(EOSMAT_AUDIT_FORMAT)
        || object.get("format_version").and_then(Value::as_u64) != Some(5)
    {
        return Err(invalid_document(format!(
            "audit sidecar for {material:?} must declare format {EOSMAT_AUDIT_FORMAT:?} and format_version 5"
        )));
    }
    if object.get("material").and_then(Value::as_str) != Some(material) {
        return Err(invalid_document(format!(
            "audit sidecar names material {}, expected {material:?}",
            object.get("material").unwrap_or(&Value::Null)
        )));
    }
    let records = object
        .get("records")
        .and_then(Value::as_object)
        .ok_or_else(|| invalid_document("audit sidecar records must be an object"))?;
    for (identifier, entry) in records {
        if !entry.as_object().is_some_and(|entry| {
            entry
                .keys()
                .all(|key| AUDIT_RECORD_KEYS.contains(&key.as_str()))
        }) {
            return Err(invalid_document(format!(
                "audit sidecar entry {identifier:?} may contain only {}",
                AUDIT_RECORD_KEYS.join(", ")
            )));
        }
    }
    Ok(())
}

/// Return the complete inline document: record file plus sidecar entries.
pub(super) fn merge(document: &Value, sidecar: &Value) -> Result<Value, EosmatError> {
    let material = identifier_text(document.get("identifier"));
    validate_sidecar(sidecar, &material)?;
    let entries = sidecar["records"]
        .as_object()
        .expect("validated sidecar records");
    let mut result = document.clone();
    let records = result
        .get_mut("eos_records")
        .and_then(Value::as_array_mut)
        .ok_or_else(|| invalid_document("eos_records must be an array"))?;
    let mut seen = BTreeSet::new();
    for record in records {
        let identifier = identifier_text(record.get("identifier"));
        let Some(entry) = entries.get(&identifier).and_then(Value::as_object) else {
            continue;
        };
        let record = record
            .as_object_mut()
            .ok_or_else(|| invalid_document("EOS record must be an object"))?;
        let collision: BTreeSet<&str> = entry
            .keys()
            .filter(|key| record.contains_key(*key))
            .map(String::as_str)
            .collect();
        if !collision.is_empty() {
            return Err(invalid_document(format!(
                "record {identifier:?} holds {collision:?} both inline and in its audit sidecar"
            )));
        }
        record.extend(expand_entry(entry, &material, &identifier));
        seen.insert(identifier);
    }
    let unknown: BTreeSet<&str> = entries
        .keys()
        .filter(|key| !seen.contains(*key))
        .map(String::as_str)
        .collect();
    if !unknown.is_empty() {
        return Err(invalid_document(format!(
            "audit sidecar for {material:?} names unknown records {unknown:?}"
        )));
    }
    Ok(result)
}

/// Merge the sidecars named by `check.audit` references, relative to `base`.
pub(super) fn merge_referenced(document: &Value, base: &Path) -> Result<Value, EosmatError> {
    let Some(records) = document.get("eos_records").and_then(Value::as_array) else {
        return Ok(document.clone());
    };
    let mut files: BTreeMap<String, Value> = BTreeMap::new();
    for record in records.iter().filter(|record| is_external(record)) {
        let reference = record["check"]["audit"].as_str().unwrap_or_default();
        let (relative, fragment) = reference.split_once('#').unwrap_or((reference, ""));
        let identifier = identifier_text(record.get("identifier"));
        if !fragment.is_empty() && fragment != identifier {
            return Err(invalid_document(format!(
                "check.audit {reference:?} does not name record {identifier:?}"
            )));
        }
        if files.contains_key(relative) {
            continue;
        }
        let path = base.join(relative);
        let text = std::fs::read_to_string(&path).map_err(|source| {
            if source.kind() == std::io::ErrorKind::NotFound {
                EosmatError::MissingAudit {
                    path: path.clone(),
                    record: identifier.clone(),
                    source,
                }
            } else {
                EosmatError::Io(source)
            }
        })?;
        files.insert(relative.to_owned(), serde_json::from_str(&text)?);
    }
    let material = identifier_text(document.get("identifier"));
    let present: BTreeSet<String> = records
        .iter()
        .map(|record| identifier_text(record.get("identifier")))
        .collect();
    let mut result = document.clone();
    for sidecar in files.values() {
        validate_sidecar(sidecar, &material)?;
        let mut selected = sidecar.clone();
        selected["records"]
            .as_object_mut()
            .expect("validated sidecar records")
            .retain(|key, _| present.contains(key));
        result = merge(&result, &selected)?;
    }
    Ok(result)
}
