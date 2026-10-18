//! Native v4 decoding. The scientific constructors share historical validation.

use serde_json::{json, Map, Value};

#[path = "shared_metadata.rs"]
mod shared_metadata;

use super::{component_model_identifier, invalid_document, EosmatError, RawComponent};

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

pub(super) fn expand(document: &Value) -> Result<Value, EosmatError> {
    let expanded_metadata = shared_metadata::expand(document)?;
    let mut result = object(&expanded_metadata, "document")?;
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
        record
            .entry("scientific_validation")
            .or_insert(json!({"status":"pending_primary_source_check"}));
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
