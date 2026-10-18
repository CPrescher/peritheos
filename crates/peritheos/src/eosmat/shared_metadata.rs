//! File-local metadata references, never inheritance or numerical substitution.
use serde_json::{Map, Value};

use super::{invalid_document, object, EosmatError};

fn resolve(
    value: &Value,
    pool: &Map<String, Value>,
    stack: &mut Vec<String>,
    nodes: &mut usize,
) -> Result<Value, EosmatError> {
    *nodes += 1;
    if (!pool.is_empty() && *nodes > 100_000) || stack.len() > 32 {
        return Err(invalid_document("shared_metadata expansion limit exceeded"));
    }
    match value {
        Value::Object(fields) => {
            if let Some(reference) = fields
                .get("shared_ref")
                .filter(|_| !pool.is_empty() || fields.len() == 1)
            {
                if fields.len() != 1 {
                    return Err(invalid_document(
                        "shared_ref cannot have sibling fields or overrides",
                    ));
                }
                let name = reference
                    .as_str()
                    .ok_or_else(|| invalid_document("shared_ref must be a name"))?;
                let target = pool.get(name).ok_or_else(|| {
                    invalid_document(format!("unknown shared_metadata reference {name}"))
                })?;
                if stack.iter().any(|key| key == name) {
                    return Err(invalid_document(format!(
                        "cyclic shared_metadata reference {name}"
                    )));
                }
                stack.push(name.to_owned());
                let resolved = resolve(target, pool, stack, nodes)?;
                stack.pop();
                return Ok(resolved);
            }
            let mut result = Map::new();
            for (key, item) in fields {
                result.insert(key.clone(), resolve(item, pool, stack, nodes)?);
            }
            Ok(Value::Object(result))
        }
        Value::Array(items) => Ok(Value::Array(
            items
                .iter()
                .map(|v| resolve(v, pool, stack, nodes))
                .collect::<Result<_, _>>()?,
        )),
        _ => Ok(value.clone()),
    }
}

fn assessment(evidence: &Value) -> Result<(), EosmatError> {
    let Some(value) = evidence
        .get("reproduction")
        .and_then(|v| v.get("assessment"))
    else {
        return Ok(());
    };
    let fields = object(value, "evidence.reproduction.assessment")?;
    let status = fields.get("status").and_then(Value::as_str).unwrap_or("");
    let scope = fields.get("scope").and_then(Value::as_str).unwrap_or("");
    if ![
        "parity",
        "similar",
        "not_reproduced",
        "not_refittable",
        "source_reconstruction",
        "not_assessed",
    ]
    .contains(&status)
    {
        return Err(invalid_document(
            "invalid evidence.reproduction.assessment.status",
        ));
    }
    if ![
        "parameter_agreement",
        "author_workflow",
        "independent_refit",
        "source_equation",
    ]
    .contains(&scope)
    {
        return Err(invalid_document(
            "invalid evidence.reproduction.assessment.scope",
        ));
    }
    if status != "not_assessed" && !fields.contains_key("basis") {
        return Err(invalid_document(
            "evidence.reproduction.assessment.basis is required",
        ));
    }
    for key in ["basis", "report", "source_status", "basis_code"] {
        if let Some(value) = fields.get(key) {
            if value.as_str().is_none_or(|s| s.trim().is_empty()) {
                return Err(invalid_document(format!(
                    "assessment.{key} must be a nonempty string"
                )));
            }
        }
    }
    if let Some(value) = fields.get("qualifications") {
        let values = value
            .as_array()
            .ok_or_else(|| invalid_document("assessment.qualifications must be an array"))?;
        if values
            .iter()
            .any(|v| v.as_str().is_none_or(|s| s.trim().is_empty()))
        {
            return Err(invalid_document(
                "assessment.qualifications must contain nonempty strings",
            ));
        }
    }
    if fields.get("fit").is_some_and(|value| !value.is_object()) {
        return Err(invalid_document("assessment.fit must be an object"));
    }
    if let Some(value) = fields.get("parameters") {
        let rows = value
            .as_array()
            .ok_or_else(|| invalid_document("assessment.parameters must be an array"))?;
        let mut names = std::collections::BTreeSet::new();
        for row in rows {
            let row = object(row, "assessment.parameters entry")?;
            let name = row.get("parameter").and_then(Value::as_str).unwrap_or("");
            if name.trim().is_empty() || !names.insert(name.to_owned()) {
                return Err(invalid_document(
                    "assessment.parameters requires unique parameter names",
                ));
            }
            for key in ["published", "refit"] {
                if !row.contains_key(key) {
                    return Err(invalid_document(format!(
                        "assessment.parameters requires {key}"
                    )));
                }
            }
            for key in [
                "published",
                "refit",
                "published_error",
                "refit_error",
                "published_sigma",
                "refit_sigma",
            ] {
                if let Some(value) = row.get(key).filter(|v| !v.is_null()) {
                    let number = value.as_f64().ok_or_else(|| {
                        invalid_document(format!(
                            "assessment.parameters.{key} must be numeric or null"
                        ))
                    })?;
                    if (key.ends_with("_error") || key.ends_with("_sigma")) && number < 0.0 {
                        return Err(invalid_document(format!(
                            "assessment.parameters.{key} must be nonnegative"
                        )));
                    }
                }
            }
        }
    }
    Ok(())
}

pub(super) fn expand(document: &Value) -> Result<Value, EosmatError> {
    let mut result = object(document, "document")?;
    if let Some(datasets) = result.get("datasets").and_then(Value::as_array) {
        for dataset in datasets {
            if dataset
                .get("reference")
                .and_then(Value::as_object)
                .is_some_and(|v| {
                    v.contains_key("shared_ref")
                        && (v.len() == 1
                            || result
                                .get("shared_metadata")
                                .and_then(Value::as_object)
                                .is_some_and(|pool| !pool.is_empty()))
                })
            {
                return Err(invalid_document(
                    "dataset references must remain inline; shared_ref is metadata-only",
                ));
            }
        }
    }
    let original_pool = result.remove("shared_metadata");
    let pool = match &original_pool {
        Some(value) => object(value, "shared_metadata")?,
        None => Map::new(),
    };
    for (name, value) in &pool {
        if !name.split('_').all(|part| {
            !part.is_empty()
                && part
                    .bytes()
                    .all(|b| b.is_ascii_lowercase() || b.is_ascii_digit())
        }) || !value.is_object()
        {
            return Err(invalid_document(
                "shared_metadata requires lower-snake-case names and object values",
            ));
        }
    }
    let mut nodes = 0;
    let mut definitions = Map::new();
    for (name, value) in &pool {
        definitions.insert(
            name.clone(),
            resolve(value, &pool, &mut vec![name.clone()], &mut nodes)?,
        );
    }
    if original_pool.is_some() {
        result.insert("shared_metadata".into(), Value::Object(definitions));
    }
    if let Some(records) = result.get_mut("eos_records").and_then(Value::as_array_mut) {
        for record in records {
            if let Some(fields) = record.as_object_mut() {
                for key in ["reference", "evidence", "pressure_calibration"] {
                    if let Some(value) = fields.get(key) {
                        let expanded = resolve(value, &pool, &mut Vec::new(), &mut nodes)?;
                        fields.insert(key.into(), expanded);
                    }
                }
                if let Some(evidence) = fields.get("evidence") {
                    assessment(evidence)?;
                }
            }
        }
    }
    Ok(Value::Object(result))
}
