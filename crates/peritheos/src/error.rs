use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// Errors produced by equation-of-state construction and evaluation.
#[derive(Clone, Debug, PartialEq)]
#[non_exhaustive]
pub enum EosError {
    /// A constructor parameter violates the model's domain.
    InvalidParameter {
        /// Public constructor parameter name.
        name: &'static str,
        /// Stable explanatory category, not a compatibility-stable sentence.
        reason: &'static str,
    },
    /// A requested pressure-volume-temperature state is invalid.
    InvalidState {
        /// Public state-variable name.
        name: &'static str,
        /// Stable explanatory category, not a compatibility-stable sentence.
        reason: &'static str,
    },
    /// The requested state is outside the invertible branch of the model.
    OutsideInvertibleRange,
    /// A positive physical bracket could not be constructed.
    BracketingFailed,
    /// A numerical solver did not satisfy its convergence requirements.
    ConvergenceFailed,
    /// Model evaluation produced a non-finite value.
    NonFiniteResult,
}

impl Display for EosError {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        match self {
            Self::InvalidParameter { name, reason } => {
                write!(formatter, "invalid parameter {name}: {reason}")
            }
            Self::InvalidState { name, reason } => {
                write!(formatter, "invalid state variable {name}: {reason}")
            }
            Self::OutsideInvertibleRange => {
                formatter.write_str("state is outside the model's invertible range")
            }
            Self::BracketingFailed => formatter.write_str("failed to bracket a root"),
            Self::ConvergenceFailed => formatter.write_str("numerical solver failed to converge"),
            Self::NonFiniteResult => formatter.write_str("model returned a non-finite result"),
        }
    }
}

impl Error for EosError {}
