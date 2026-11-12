use serde_json::{json, Value};

use super::{invalid_document, is_molar_volume_model, EosmatError};

fn volume_scale(document: &Value, record: &Value) -> Result<f64, EosmatError> {
    let scale = if let Some(scale) = record["volume"]["public_to_model_scale"].as_f64() {
        scale
    } else if record["thermal"]["model"]
        .as_str()
        .is_some_and(is_molar_volume_model)
        || record["eos"]["model"].as_str() == Some("maltby_2024_published")
    {
        let z = document["formula_units_per_cell"].as_f64().ok_or_else(|| {
            invalid_document("uncertainty conversion requires formula_units_per_cell")
        })?;
        super::CELL_ANGSTROM3_TO_FORMULA_MOLAR_J_PER_BAR / z
    } else {
        1.0
    };
    if !scale.is_finite() || scale <= 0.0 {
        return Err(invalid_document(
            "uncertainty volume scale must be positive and finite",
        ));
    }
    Ok(scale)
}

#[allow(clippy::float_cmp)] // Identity conversion must retain literal JSON values.
fn scaled(value: &Value, factor: f64) -> Result<Value, EosmatError> {
    if factor == 1.0 {
        return Ok(value.clone());
    }
    let original = value
        .as_f64()
        .ok_or_else(|| invalid_document("uncertainty must be numeric"))?;
    let result = original * factor;
    if !result.is_finite() || original != 0.0 && result == 0.0 {
        return Err(invalid_document(
            "uncertainty unit conversion overflowed or underflowed",
        ));
    }
    Ok(json!(result))
}

/// Convert an expanded document to the historical model-unit uncertainty contract.
pub(super) fn to_model_units(document: &mut Value, from_v5: bool) -> Result<(), EosmatError> {
    let count = document["eos_records"].as_array().map_or(0, Vec::len);
    for index in 0..count {
        let record = &document["eos_records"][index];
        let error_public = from_v5 || record["parameter_error_volume_basis"] == "public";
        let covariance_public =
            from_v5 || record["parameter_covariance"]["volume_basis"] == "public";
        let error_present = record["parameter_errors"]
            .get("V0")
            .is_some_and(|v| !v.is_null());
        let order = record["parameter_covariance"]["parameter_order"].as_array();
        let volume_names = order.map(|names| {
            names
                .iter()
                .map(|v| matches!(v.as_str(), Some("V0" | "rt_eos.V0")))
                .collect::<Vec<_>>()
        });
        let covariance_has_volume = volume_names
            .as_ref()
            .is_some_and(|names| names.iter().any(|v| *v));
        let scale = if error_public && error_present || covariance_public && covariance_has_volume {
            volume_scale(document, record)?
        } else {
            1.0
        };
        let record = document["eos_records"][index]
            .as_object_mut()
            .ok_or_else(|| invalid_document("record must be an object"))?;
        record.remove("parameter_error_volume_basis");
        if error_public && error_present {
            let error = &mut record.get_mut("parameter_errors").unwrap()["V0"];
            *error = scaled(error, scale)?;
        }
        if let Some(covariance) = record
            .get_mut("parameter_covariance")
            .and_then(Value::as_object_mut)
        {
            covariance.remove("volume_basis");
            if covariance_public && covariance_has_volume {
                let names = volume_names.as_ref().unwrap();
                let matrix = covariance
                    .get_mut("matrix")
                    .and_then(Value::as_array_mut)
                    .ok_or_else(|| invalid_document("covariance matrix must be an array"))?;
                for (i, row) in matrix.iter_mut().enumerate() {
                    let row = row
                        .as_array_mut()
                        .ok_or_else(|| invalid_document("covariance row must be an array"))?;
                    for (j, value) in row.iter_mut().enumerate() {
                        let first = scaled(value, if names[i] { scale } else { 1.0 })?;
                        *value = scaled(&first, if names[j] { scale } else { 1.0 })?;
                    }
                }
            }
        }
    }
    Ok(())
}
