//! Native v4 decoding. The scientific constructors share historical validation.

use serde_json::{json, Map, Value};

#[path = "shared_metadata.rs"]
mod shared_metadata;

use super::{audit, component_model_identifier, invalid_document, EosmatError, RawComponent};

const EVIDENCE_FIELDS: &[&str] = &[
    "scientific_validation",
    "parameter_provenance",
    "source_lineage",
    "fit_provenance",
    "fit_protocol",
    "reproduction",
    "audit_corrections",
];

const EVIDENCE_RESERVED: &[&str] = &[
    "derivation",
    "notes",
    "aliases",
    "identifier",
    "label",
    "record_kind",
    "determination_method",
    "reference",
    "default",
    "default_for",
    "equation_kind",
    "loading_path",
    "branch_kind",
    "initial_state",
    "volume_basis",
    "branch_domain",
    "derived_from_record",
    "pressure_scale_variant_of",
    "eos",
    "thermal",
    "parameter_errors",
    "fixed_parameters",
    "parameter_error_confidence",
    "parameter_error_volume_basis",
    "parameter_covariance",
    "volume",
    "validity",
    "experimental_pressure_range_gpa",
    "pressure_range_status",
    "experimental_temperature_range_k",
    "temperature_ref",
    "temperature_ref_provenance",
    "pressure_calibration",
    "catalog_access",
    "model",
    "parameters",
    "equation_options",
    "evidence",
    "check",
];

fn object(value: &Value, label: &str) -> Result<Map<String, Value>, EosmatError> {
    value
        .as_object()
        .cloned()
        .ok_or_else(|| invalid_document(format!("{label} must be an object")))
}

fn model(kind: &Value, thermal: bool) -> Result<String, EosmatError> {
    let kind = kind
        .as_str()
        .ok_or_else(|| invalid_document("model must be a supported name"))?;
    let raw: RawComponent = serde_json::from_value(json!({"type": kind}))?;
    component_model_identifier(&raw, thermal)
        .map(str::to_owned)
        .map_err(invalid_document)
}

type Parameters = (Map<String, Value>, Map<String, Value>, Vec<Value>);

fn parameters(value: &Value) -> Result<Parameters, EosmatError> {
    let input = object(value, "parameters")?;
    let (mut values, mut errors, mut fixed) = (Map::new(), Map::new(), Vec::new());
    for (name, value) in input {
        if let Some(annotation) = value.as_object() {
            if annotation.is_empty()
                || annotation
                    .keys()
                    .any(|k| !["value", "error", "fixed"].contains(&k.as_str()))
            {
                return Err(invalid_document(format!(
                    "invalid parameter annotation for {name}"
                )));
            }
            if let Some(value) = annotation.get("value") {
                values.insert(name.clone(), value.clone());
            }
            if let Some(value) = annotation.get("error") {
                errors.insert(name.clone(), value.clone());
            }
            if let Some(value) = annotation.get("fixed") {
                match value.as_bool() {
                    Some(true) => fixed.push(Value::String(name.clone())),
                    Some(false) if annotation.len() == 1 => {
                        return Err(invalid_document("parameter annotation has no information"))
                    }
                    Some(false) => (),
                    None => return Err(invalid_document("parameter fixed flag must be boolean")),
                }
            }
        } else {
            values.insert(name, value);
        }
    }
    Ok((values, errors, fixed))
}

/// Default-selection categories of `default_records`.
const DEFAULT_CATEGORIES: [&str; 2] = ["equilibrium", "hugoniot"];

/// `hugoniot` for a Hugoniot model, otherwise `equilibrium`.
pub(super) fn default_category(record: &Map<String, Value>) -> &'static str {
    let model = match record.get("eos").and_then(Value::as_object) {
        Some(eos) => eos.get("type"),
        None => record.get("model"),
    };
    if model.and_then(Value::as_str) == Some("LinearUsUpHugoniot") {
        "hugoniot"
    } else {
        "equilibrium"
    }
}

fn flags_default(record: &Map<String, Value>) -> bool {
    record.get("default") == Some(&Value::Bool(true))
        || record
            .get("default_for")
            .is_some_and(|value| !value.is_null())
}

fn quoted(value: &Value) -> String {
    value
        .as_str()
        .map_or_else(|| value.to_string(), |text| format!("'{text}'"))
}

/// Write a material `default_records` map as per-record flags, in place.
///
/// The named record gets `default_for: <category>` and the map is removed.
/// Per-record flags that are already present must agree with the map: a
/// named record flagged `default: false` and any other record flagged as a
/// default are errors, never silently resolved. Mirrors
/// `expand_default_records` in `peritheos/_eosmat_v4.py`.
fn expand_default_records(document: &mut Map<String, Value>) -> Result<(), EosmatError> {
    let Some(mapping) = document.remove("default_records") else {
        return Ok(());
    };
    let mapping = mapping
        .as_object()
        .filter(|mapping| !mapping.is_empty())
        .cloned()
        .ok_or_else(|| invalid_document("default_records must be a non-empty object"))?;
    let unknown: Vec<&str> = mapping
        .keys()
        .map(String::as_str)
        .filter(|key| !DEFAULT_CATEGORIES.contains(key))
        .collect();
    if !unknown.is_empty() {
        let mut unknown = unknown;
        unknown.sort_unstable();
        return Err(invalid_document(format!(
            "default_records keys must be equilibrium or hugoniot; got {unknown:?}"
        )));
    }
    let records = document
        .get_mut("eos_records")
        .and_then(Value::as_array_mut)
        .ok_or_else(|| invalid_document("eos_records must be an array"))?;
    for (category, identifier) in &mapping {
        let Some(identifier) = identifier.as_str().filter(|text| !text.is_empty()) else {
            return Err(invalid_document(format!(
                "default_records.{category} must be a record identifier"
            )));
        };
        let record = records
            .iter()
            .filter_map(Value::as_object)
            .find(|record| record.get("identifier").and_then(Value::as_str) == Some(identifier))
            .ok_or_else(|| {
                invalid_document(format!(
                    "default_records.{category} names unknown record '{identifier}'"
                ))
            })?;
        if default_category(record) != category {
            return Err(invalid_document(format!(
                "default_records.{category} names '{identifier}', which is not a {category} record"
            )));
        }
    }
    for record in records.iter_mut().filter_map(Value::as_object_mut) {
        let category = default_category(record);
        let identifier = record.get("identifier").cloned().unwrap_or(Value::Null);
        let named = mapping.get(category);
        if named.is_some_and(|named| *named == identifier) {
            if record.get("default") == Some(&Value::Bool(false)) {
                return Err(invalid_document(format!(
                    "default_records.{category} names {}, which is flagged default: false",
                    quoted(&identifier)
                )));
            }
            record.insert("default_for".into(), json!(category));
        } else if flags_default(record) {
            return Err(invalid_document(format!(
                "Record {} is flagged as the default {category} record, but default_records names {}",
                quoted(&identifier),
                named.map_or_else(|| format!("no {category} record"), quoted)
            )));
        }
    }
    Ok(())
}

/// Give every record an explicit `default` boolean, for v3 exports of v5.
///
/// The default of each category gets `default: true` and `default_for:
/// <category>`, every other record `default: false`, so historical readers
/// see the selection on every record.
pub(super) fn write_explicit_default_flags(document: &mut Value) {
    let Some(records) = document
        .get_mut("eos_records")
        .and_then(Value::as_array_mut)
    else {
        return;
    };
    for record in records.iter_mut().filter_map(Value::as_object_mut) {
        let category = default_category(record);
        let selected = record.get("default") == Some(&Value::Bool(true))
            || record.get("default_for").and_then(Value::as_str) == Some(category);
        record.insert("default".into(), Value::Bool(selected));
        if selected {
            record.insert("default_for".into(), json!(category));
        }
    }
}

#[allow(clippy::too_many_lines)]
pub(super) fn expand(document: &Value) -> Result<Value, EosmatError> {
    let parameter_units = document["format_version"].as_u64() == Some(5);
    let expanded_metadata = shared_metadata::expand(document)?;
    let mut result = object(&expanded_metadata, "document")?;
    if let Some(block) = result.get("superseded_material_text") {
        audit::validate_superseded_material_text(block, "material")?;
    }
    if result.contains_key("default_records") && !parameter_units {
        return Err(invalid_document(
            "default_records requires EOSMAT format_version 5",
        ));
    }
    expand_default_records(&mut result)?;
    if result.contains_key("format") {
        return Err(invalid_document("EOSMAT v4 uses format_version alone"));
    }
    let identifier = result
        .get("identifier")
        .and_then(Value::as_str)
        .filter(|v| !v.is_empty())
        .ok_or_else(|| invalid_document("EOSMAT v4 requires identifier"))?
        .to_owned();
    result.insert("format".into(), json!("peritheos.material"));
    result.insert("format_version".into(), json!(3));
    result.entry("name").or_insert(json!(identifier));
    result.entry("units").or_insert(
        json!({"pressure":"GPa","temperature":"K","volume":"angstrom^3/conventional_unit_cell"}),
    );
    let records = result
        .get("eos_records")
        .and_then(Value::as_array)
        .ok_or_else(|| invalid_document("eos_records must be an array"))?;
    let mut expanded = Vec::new();
    for source in records {
        let mut record = object(source, "EOS record")?;
        if parameter_units
            && (record.contains_key("parameter_error_volume_basis")
                || record
                    .get("parameter_covariance")
                    .is_some_and(|c| c.get("volume_basis").is_some()))
        {
            return Err(invalid_document(
                "EOSMAT v5 uncertainties always use parameter units; remove legacy unit flags",
            ));
        }
        if ["eos", "parameter_errors", "fixed_parameters"]
            .iter()
            .chain(EVIDENCE_FIELDS)
            .any(|key| record.contains_key(*key))
        {
            return Err(invalid_document(
                "EOSMAT v4 record mixes native and v3 fields",
            ));
        }
        let kind = record.remove("model").unwrap_or(Value::Null);
        let canonical = model(&kind, false)?;
        let (p, e, f) = parameters(&record.remove("parameters").unwrap_or(Value::Null))?;
        let mut options = object(
            &record.remove("equation_options").unwrap_or(json!({})),
            "equation_options",
        )?;
        if ["type", "model", "parameters"]
            .iter()
            .any(|key| options.contains_key(*key))
        {
            return Err(invalid_document(
                "equation_options cannot replace model or parameters",
            ));
        }
        options.insert("type".into(), kind);
        options.insert("model".into(), json!(canonical));
        options.insert("parameters".into(), Value::Object(p));
        record.insert("eos".into(), Value::Object(options));
        record.insert("parameter_errors".into(), Value::Object(e));
        record.insert("fixed_parameters".into(), Value::Array(f));
        if let Some(check) = source.get("check") {
            let location = format!(
                "eos_records[{}]",
                source.get("identifier").unwrap_or(&Value::Null)
            );
            audit::validate_check(check, &location)?;
            audit::validate_check_consistency(source, &location)?;
        }
        if let Some(block) = source.get("superseded_record_text") {
            let location = format!(
                "eos_records[{}]",
                source.get("identifier").unwrap_or(&Value::Null)
            );
            audit::validate_superseded_record_text(block, &location)?;
        }
        audit::validate_check_text(
            source,
            Some(document),
            &format!(
                "eos_records[{}]",
                source.get("identifier").unwrap_or(&Value::Null)
            ),
        )?;
        let evidence = object(&record.remove("evidence").unwrap_or(json!({})), "evidence")?;
        if evidence
            .keys()
            .any(|key| record.contains_key(key) || EVIDENCE_RESERVED.contains(&key.as_str()))
        {
            return Err(invalid_document(
                "evidence fields collide with EOS record fields",
            ));
        }
        record.extend(evidence);
        let label = record.get("identifier").cloned().unwrap_or(json!(""));
        record.entry("label").or_insert(label);
        record
            .entry("reference")
            .or_insert(json!("User-supplied parameters"));
        if !record.contains_key("scientific_validation") {
            // A record file whose evidence lives in an audit sidecar keeps the
            // stored source-check status in its headline (a fixed renaming of
            // the sidecar code); nothing is inferred.
            let status = record
                .get("check")
                .and_then(|check| check.get("source_check"))
                .and_then(Value::as_str)
                .and_then(audit::SourceCheck::from_check)
                .map_or(
                    "pending_primary_source_check",
                    audit::SourceCheck::status_code,
                );
            record.insert("scientific_validation".into(), json!({ "status": status }));
        }
        if let Some(value) = record.get("thermal") {
            let mut thermal = object(value, "thermal")?;
            if ["type", "parameter_errors", "fixed_parameters"]
                .iter()
                .any(|key| thermal.contains_key(*key))
            {
                return Err(invalid_document("native thermal component mixes v3 fields"));
            }
            let kind = thermal.get("model").cloned().unwrap_or(Value::Null);
            let canonical = model(&kind, true)?;
            let (p, e, f) = parameters(thermal.get("parameters").unwrap_or(&Value::Null))?;
            thermal.insert("type".into(), kind);
            thermal.insert("model".into(), json!(canonical));
            thermal.insert("parameters".into(), Value::Object(p));
            if !e.is_empty() {
                thermal.insert("parameter_errors".into(), Value::Object(e));
            }
            if !f.is_empty() {
                thermal.insert("fixed_parameters".into(), Value::Array(f));
            }
            record.insert("thermal".into(), Value::Object(thermal));
        }
        expanded.push(Value::Object(record));
    }
    result.insert("eos_records".into(), Value::Array(expanded));
    Ok(Value::Object(result))
}
