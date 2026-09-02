use std::error::Error;
use std::fmt::{self, Display, Formatter};

/// Machine-readable category for an [`EosError`].
///
/// Variant names and [`EosError::code`] are stable API; display strings are
/// diagnostic text and may become more detailed.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
#[non_exhaustive]
pub enum EosErrorKind {
    InvalidParameter,
    InvalidState,
    OutsideInvertibleRange,
    BracketingFailed,
    ConvergenceFailed,
    NonFiniteResult,
}

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

impl EosError {
    /// Return the machine-readable error category.
    #[must_use]
    pub const fn kind(&self) -> EosErrorKind {
        match self {
            Self::InvalidParameter { .. } => EosErrorKind::InvalidParameter,
            Self::InvalidState { .. } => EosErrorKind::InvalidState,
            Self::OutsideInvertibleRange => EosErrorKind::OutsideInvertibleRange,
            Self::BracketingFailed => EosErrorKind::BracketingFailed,
            Self::ConvergenceFailed => EosErrorKind::ConvergenceFailed,
            Self::NonFiniteResult => EosErrorKind::NonFiniteResult,
        }
    }

    /// Return a stable, language-independent error code.
    #[must_use]
    pub const fn code(&self) -> &'static str {
        match self {
            Self::InvalidParameter { .. } => "eos.invalid_parameter",
            Self::InvalidState { .. } => "eos.invalid_state",
            Self::OutsideInvertibleRange => "eos.outside_invertible_range",
            Self::BracketingFailed => "eos.bracketing_failed",
            Self::ConvergenceFailed => "eos.convergence_failed",
            Self::NonFiniteResult => "eos.non_finite_result",
        }
    }

    /// Return the invalid public parameter or state-variable name, if known.
    #[must_use]
    pub const fn field(&self) -> Option<&'static str> {
        match self {
            Self::InvalidParameter { name, .. } | Self::InvalidState { name, .. } => Some(name),
            Self::OutsideInvertibleRange
            | Self::BracketingFailed
            | Self::ConvergenceFailed
            | Self::NonFiniteResult => None,
        }
    }

    /// Whether this error represents invalid input rather than numerical failure.
    #[must_use]
    pub const fn is_validation(&self) -> bool {
        matches!(
            self,
            Self::InvalidParameter { .. }
                | Self::InvalidState { .. }
                | Self::OutsideInvertibleRange
                | Self::BracketingFailed
        )
    }
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
