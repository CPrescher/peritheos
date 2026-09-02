use crate::{EosError, EosResult};

pub(crate) fn finite_parameter(value: f64, name: &'static str) -> EosResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EosError::InvalidParameter {
            name,
            reason: "must be finite",
        })
    }
}

pub(crate) fn positive_parameter(value: f64, name: &'static str) -> EosResult<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(EosError::InvalidParameter {
            name,
            reason: "must be positive and finite",
        })
    }
}

pub(crate) fn finite_state(value: f64, name: &'static str) -> EosResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EosError::InvalidState {
            name,
            reason: "must be finite",
        })
    }
}

pub(crate) fn positive_state(value: f64, name: &'static str) -> EosResult<f64> {
    if value.is_finite() && value > 0.0 {
        Ok(value)
    } else {
        Err(EosError::InvalidState {
            name,
            reason: "must be positive and finite",
        })
    }
}

pub(crate) fn finite_result(value: f64) -> EosResult<f64> {
    if value.is_finite() {
        Ok(value)
    } else {
        Err(EosError::NonFiniteResult)
    }
}
