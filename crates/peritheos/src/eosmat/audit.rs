//! EOSMAT v5 record `check` headlines and per-material audit sidecars.
//!
//! This mirrors `peritheos/_eosmat_audit.py`. A bundled material is stored as
//! a record file (`<id>.eosmat`, whose records end with a short `check` block)
//! and an audit sidecar (`<id>.audit.json`) holding each record's evidence,
//! detailed `pressure_calibration` and other audit-only blocks. Merging the
//! sidecar restores the complete inline v5 document. Nothing here reclassifies
//! a record: `check.reproduction` is a fixed mapping of the stored assessment
//! status, `check.compared` a fixed mapping of its status and scope,
//! `check.why` is composed only from the stored coefficient comparison table,
//! `check.source_check` is a fixed renaming of the stored source-validation
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
    "superseded_record_text",
];

/// Free-text record fields whose replaced original text a record keeps in its
/// `superseded_record_text` audit block: `label`, `notes`, `validity_notes`
/// (`validity.notes`) and `thermal_notes` (`thermal.notes`). `null` means the
/// field was absent before.
const SUPERSEDED_TEXT_FIELDS: &[&str] = &["label", "notes", "validity_notes", "thermal_notes"];

const CHECK_KEYS: &[&str] = &[
    "reproduction",
    "compared",
    "why",
    "warning",
    "source_check",
    "source_checked_on",
    "source_note",
    "audit",
];
/// Earlier names of v5 `check` fields. They are rejected, never reinterpreted.
const RENAMED_CHECK_KEYS: &[(&str, &str)] =
    &[("finding", "source_note"), ("source", "source_check")];
/// Maximum length of a `check.why` line, in characters.
pub const WHY_MAX_LENGTH: usize = 160;

const REPORT_LEDGER: &str = "docs/data/primary-eos-refits.json";

const DO_NOT_USE_WARNING: &str = "Do not use.";
const HOLD_WARNING: &str = "Default selection on hold pending review.";
const EXPLICIT_SELECTION_WARNING: &str = "Never chosen automatically; select it explicitly.";

/// Headline outcome of a record's reconstruction assessment (`check.reproduction`).
///
/// This is a fixed mapping of `evidence.reproduction.assessment.status`:
/// `parity` is [`Reproduction::Reproduced`], `similar` is
/// [`Reproduction::Close`], `not_reproduced` is [`Reproduction::Disagrees`],
/// and every other or absent status is [`Reproduction::Untested`]. What was
/// compared is [`Compared`]. It is not a physical validation or a
/// recommendation for use.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Reproduction {
    /// The comparison named by [`Compared`] passed (`parity`).
    Reproduced,
    /// The comparison agreed within the stated similarity scope (`similar`).
    Close,
    /// The comparison did not agree (`not_reproduced`).
    Disagrees,
    /// No comparison gave a result, or none could be attempted.
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

/// What a record's [`Reproduction`] headline refers to (`check.compared`).
///
/// A fixed mapping of the stored assessment status and scope: no assessment,
/// or status `not_refittable` or `not_assessed`, is [`Compared::Nothing`];
/// status `source_reconstruction` is [`Compared::SourceEquation`]; otherwise
/// scope `parameter_agreement` or `author_workflow` is
/// [`Compared::PublishedFit`], `independent_refit` is
/// [`Compared::PeritheosRefit`] and `source_equation` is
/// [`Compared::SourceEquation`].
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum Compared {
    /// A refit of the source data compared with the published coefficients.
    PublishedFit,
    /// A rerun of a stored Peritheos refit under its own inputs; it says
    /// nothing about the original paper's coefficients.
    PeritheosRefit,
    /// The source equation or its calculated checkpoints were rebuilt.
    SourceEquation,
    /// Nothing could be compared (`none`).
    Nothing,
}

impl Compared {
    /// The EOSMAT `check.compared` spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::PublishedFit => "published_fit",
            Self::PeritheosRefit => "peritheos_refit",
            Self::SourceEquation => "source_equation",
            Self::Nothing => "none",
        }
    }

    /// Parse a `check.compared` value.
    #[must_use]
    pub fn from_check(value: &str) -> Option<Self> {
        match value {
            "published_fit" => Some(Self::PublishedFit),
            "peritheos_refit" => Some(Self::PeritheosRefit),
            "source_equation" => Some(Self::SourceEquation),
            "none" => Some(Self::Nothing),
            _ => None,
        }
    }
}

impl Display for Compared {
    fn fmt(&self, formatter: &mut Formatter<'_>) -> fmt::Result {
        formatter.write_str(self.as_str())
    }
}

/// Primary-source review of a record (`check.source_check`).
///
/// A fixed renaming of `scientific_validation.status`, whose original codes
/// the audit sidecar and v3 exports keep.
#[derive(Clone, Copy, Debug, PartialEq, Eq, Hash)]
#[non_exhaustive]
pub enum SourceCheck {
    /// The values were checked against the primary source
    /// (`primary_source_validated`).
    Verified,
    /// The primary-source check is outstanding (`pending_primary_source_check`).
    Pending,
    /// The check found a discrepancy that is not resolved (`not_reproduced`).
    Unresolved,
    /// The check was deliberately deferred (`deferred`).
    Deferred,
}

impl SourceCheck {
    const ALL: [Self; 4] = [
        Self::Verified,
        Self::Pending,
        Self::Unresolved,
        Self::Deferred,
    ];

    /// The EOSMAT `check.source_check` spelling.
    #[must_use]
    pub const fn as_str(self) -> &'static str {
        match self {
            Self::Verified => "verified",
            Self::Pending => "pending",
            Self::Unresolved => "unresolved",
            Self::Deferred => "deferred",
        }
    }

    /// The `scientific_validation.status` code this value renames.
    #[must_use]
    pub const fn status_code(self) -> &'static str {
        match self {
            Self::Verified => "primary_source_validated",
            Self::Pending => "pending_primary_source_check",
            Self::Unresolved => "not_reproduced",
            Self::Deferred => "deferred",
        }
    }

    /// Parse a `check.source_check` value.
    #[must_use]
    pub fn from_check(value: &str) -> Option<Self> {
        Self::ALL.into_iter().find(|check| check.as_str() == value)
    }

    /// Rename a stored `scientific_validation.status` code.
    #[must_use]
    pub fn from_status_code(code: &str) -> Option<Self> {
        Self::ALL
            .into_iter()
            .find(|check| check.status_code() == code)
    }
}

impl Display for SourceCheck {
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

/// Map the stored assessment status and scope (see [`Compared`]).
pub(super) fn compared_scope(record: &Value) -> Compared {
    let Some(assessment) = assessment(record) else {
        return Compared::Nothing;
    };
    let status = match assessment.get("status") {
        None | Some(Value::Null) => return Compared::Nothing,
        Some(Value::String(text)) => text.as_str(),
        // Python maps `str(status)`; no other value is a known status.
        Some(_) => "",
    };
    match status {
        "not_refittable" | "not_assessed" => Compared::Nothing,
        "source_reconstruction" => Compared::SourceEquation,
        _ => match assessment.get("scope").and_then(Value::as_str) {
            Some("parameter_agreement" | "author_workflow") => Compared::PublishedFit,
            Some("independent_refit") => Compared::PeritheosRefit,
            Some("source_equation") => Compared::SourceEquation,
            _ => Compared::Nothing,
        },
    }
}

fn why_unit(parameter: &str) -> &'static str {
    match parameter {
        "K0" => " GPa",
        "K0_double_prime" => " GPa^-1",
        "theta0" => " K",
        _ => "",
    }
}

fn parameter_text(row: &Map<String, Value>) -> String {
    match row.get("parameter") {
        Some(Value::String(text)) => text.clone(),
        Some(other) => other.to_string(),
        None => "None".to_owned(),
    }
}

fn parameter_name(name: &str) -> String {
    let name = name.strip_prefix("rt_eos.").unwrap_or(name);
    if let Some(base) = name.strip_suffix("_double_prime") {
        format!("{base}''")
    } else if let Some(base) = name.strip_suffix("_prime") {
        format!("{base}'")
    } else {
        name.to_owned()
    }
}

fn number(value: Option<&Value>) -> Option<f64> {
    value
        .and_then(Value::as_number)
        .and_then(serde_json::Number::as_f64)
}

/// Python's exponent spelling: a sign and at least two digits.
fn python_exponent(text: &str) -> String {
    let Some((mantissa, exponent)) = text.split_once('e') else {
        return text.to_owned();
    };
    let value: i32 = exponent.parse().unwrap_or_default();
    let sign = if value < 0 { '-' } else { '+' };
    format!("{mantissa}e{sign}{:02}", value.unsigned_abs())
}

/// A stored number as Python's `json.dumps` prints it, without a spurious
/// trailing `.0` for an integral float.
fn printed(value: &Value) -> String {
    let Some(number) = value.as_number() else {
        return value.to_string();
    };
    if !number.is_f64() {
        return number.to_string();
    }
    let float = number.as_f64().unwrap_or_default();
    if float.fract() == 0.0 && float.abs() < 1e15 {
        // Exact: an integral float below 1e15 converts without loss.
        #[allow(clippy::cast_possible_truncation)]
        return format!("{}", float as i64);
    }
    // Python's repr: shortest round-trip digits, scientific below 1e-4 and
    // from 1e16.
    let scientific = format!("{float:e}");
    let exponent: i32 = scientific
        .split_once('e')
        .and_then(|(_, exponent)| exponent.parse().ok())
        .unwrap_or_default();
    if (-4..16).contains(&exponent) {
        let text = format!("{float}");
        if text.contains('.') {
            text
        } else {
            format!("{text}.0")
        }
    } else {
        python_exponent(&scientific)
    }
}

/// Round to at most one digit more than the printed published value.
fn rounded_refit(refit: f64, published: &str) -> String {
    let lower = published.to_lowercase();
    if let Some((mantissa, _)) = lower.split_once('e') {
        let digits = mantissa
            .chars()
            .filter(|c| *c != '-' && *c != '.')
            .skip_while(|c| *c == '0')
            .count();
        return python_exponent(&format!("{refit:.digits$e}"));
    }
    let decimals = lower
        .split_once('.')
        .map_or(0, |(_, fraction)| fraction.len());
    let text = format!("{refit:.prec$}", prec = decimals + 1);
    let text = if text.contains('.') {
        text.trim_end_matches('0').trim_end_matches('.')
    } else {
        &text
    };
    if text == "-0" { "0" } else { text }.to_owned()
}

/// `a`, `a and b`, `a, b and c`.
fn join(items: &[String]) -> String {
    match items {
        [] => String::new(),
        [only] => only.clone(),
        [rest @ .., last] => format!("{} and {last}", rest.join(", ")),
    }
}

/// Phrases for the stored criterion flags that are `false` in one row.
///
/// `similar` is the audit's own agreement screen (the generic limits or the
/// criterion recorded in the assessment basis), so it is reported as "beyond
/// the agreement limit" without naming a limit.
fn failed_criteria(row: &Map<String, Value>) -> Vec<String> {
    [
        ("within_combined_2sigma", "outside the combined error"),
        ("similar", "beyond the agreement limit"),
        ("within_reported_error", "outside the published error"),
    ]
    .into_iter()
    .filter(|(flag, _)| row.get(*flag) == Some(&Value::Bool(false)))
    .map(|(_, phrase)| phrase.to_owned())
    .collect()
}

fn relative_change(row: &Map<String, Value>) -> Option<f64> {
    let refit = number(row.get("refit"))?;
    let base = number(row.get("published"))?;
    (base != 0.0).then(|| (refit - base) / base.abs())
}

fn why_row(row: &Map<String, Value>, numbers: bool) -> String {
    let parameter = parameter_text(row);
    let name = parameter_name(&parameter);
    let change = relative_change(row).map(|relative| format!("{:+.1}%", relative * 100.0));
    if !numbers {
        return change.map_or_else(|| name.clone(), |change| format!("{name} ({change})"));
    }
    let (Some(refit), Some(published)) = (
        number(row.get("refit")),
        row.get("published").filter(|value| value.is_number()),
    ) else {
        return name;
    };
    let unit = why_unit(parameter.strip_prefix("rt_eos.").unwrap_or(&parameter));
    let published = printed(published);
    let text = format!(
        "{name}: refit {} vs published {published}{unit}",
        rounded_refit(refit, &published)
    );
    match change {
        Some(change) => format!("{text} ({change})"),
        None => text,
    }
}

/// At most this many failing coefficients are named; the rest are counted.
const WHY_NAMED: usize = 3;

type Failing<'a> = (&'a Map<String, Value>, Vec<String>);

/// Candidate texts for the failing rows, most detailed first.
fn failing_lines(mut failing: Vec<Failing<'_>>) -> Vec<String> {
    let hidden = failing.len().saturating_sub(WHY_NAMED);
    if hidden > 0 {
        // Name the largest relative differences; rows without one come last.
        let size = |item: &Failing<'_>| relative_change(item.0).map_or(-1.0, f64::abs);
        failing.sort_by(|left, right| size(right).total_cmp(&size(left)));
        failing.truncate(WHY_NAMED);
    }
    let mut groups: Vec<(String, Vec<&Map<String, Value>>)> = Vec::new();
    for (row, failed) in &failing {
        let phrase = join(failed);
        match groups.iter_mut().find(|(key, _)| *key == phrase) {
            Some((_, rows)) => rows.push(row),
            None => groups.push((phrase, vec![row])),
        }
    }
    let more = format!("{hidden} more");
    let grouped = |name: &dyn Fn(&Map<String, Value>) -> String| {
        if hidden > 0 && groups.len() == 1 {
            let (failed, rows) = &groups[0];
            let mut names: Vec<String> = rows.iter().map(|row| name(row)).collect();
            names.push(more.clone());
            return format!("{}: {failed}", join(&names));
        }
        let text = groups
            .iter()
            .map(|(failed, rows)| {
                let names: Vec<String> = rows.iter().map(|row| name(row)).collect();
                format!("{}: {failed}", join(&names))
            })
            .collect::<Vec<_>>()
            .join("; ");
        if hidden > 0 {
            format!("{text}; and {more}")
        } else {
            text
        }
    };
    let mut candidates = Vec::new();
    if hidden == 0 {
        candidates.push(
            failing
                .iter()
                .map(|(row, failed)| format!("{}, {}", why_row(row, true), join(failed)))
                .collect::<Vec<_>>()
                .join("; "),
        );
    }
    candidates.push(grouped(&|row| why_row(row, false)));
    candidates.push(grouped(&|row| parameter_name(&parameter_text(row))));
    candidates
}

/// One generated line saying which coefficients kept a refit from parity.
///
/// Only for [`Reproduction::Close`] and [`Reproduction::Disagrees`]. Built
/// from the stored `assessment.parameters` rows and their criterion flags
/// (`within_combined_2sigma`, `similar`, `within_reported_error`); no prose
/// is interpreted. Published values print as stored (without a trailing
/// `.0`); the refit is rounded to at most one digit more. More than three
/// failing coefficients: the three largest relative differences are named
/// and the rest counted. `None` when the table supports no statement within
/// [`WHY_MAX_LENGTH`] characters. This mirrors `reproduction_why` in
/// `peritheos/_eosmat_audit.py`.
pub(super) fn reproduction_why(record: &Value) -> Option<String> {
    let headline = reproduction_headline(record);
    if !matches!(headline, Reproduction::Close | Reproduction::Disagrees) {
        return None;
    }
    let rows: Vec<&Map<String, Value>> = assessment(record)?
        .get("parameters")?
        .as_array()?
        .iter()
        .filter_map(Value::as_object)
        .filter(|row| row.contains_key("parameter"))
        .collect();
    let failing: Vec<Failing<'_>> = rows
        .iter()
        .map(|row| (*row, failed_criteria(row)))
        .filter(|(_, failed)| !failed.is_empty())
        .collect();
    if !failing.is_empty() {
        return failing_lines(failing)
            .into_iter()
            .find(|text| text.chars().count() < WHY_MAX_LENGTH)
            .map(|text| text + ".");
    }
    if headline != Reproduction::Close {
        return None;
    }
    let untested: Vec<String> = rows
        .iter()
        .filter(|row| row.get("within_combined_2sigma") == Some(&Value::Null))
        .map(|row| parameter_name(&parameter_text(row)))
        .collect();
    if untested.is_empty() {
        return None;
    }
    let verb = if untested.len() == 1 { "has" } else { "have" };
    let mut text = format!(
        "{} {verb} no error to test against the combined-error rule",
        join(&untested)
    );
    if rows
        .iter()
        .all(|row| row.get("similar") == Some(&Value::Bool(true)))
    {
        text.push_str("; all within the agreement limit");
    }
    text.push('.');
    (text.chars().count() <= WHY_MAX_LENGTH).then_some(text)
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

fn validation(record: &Value) -> Option<&Map<String, Value>> {
    section(evidence(record), "scientific_validation")
}

fn source_finding(validation: Option<&Map<String, Value>>) -> Option<&str> {
    validation?
        .get("primary_source_check")?
        .as_object()?
        .get("finding")?
        .as_str()
        .filter(|finding| !finding.is_empty())
}

/// Derive a record's `check` block from its inline evidence.
pub(super) fn check_block(record: &Value, material: Option<&Value>, audit: Option<&str>) -> Value {
    let validation = validation(record);
    let mut result = Map::new();
    result.insert(
        "reproduction".into(),
        Value::String(reproduction_headline(record).as_str().to_owned()),
    );
    result.insert(
        "compared".into(),
        Value::String(compared_scope(record).as_str().to_owned()),
    );
    if let Some(why) = reproduction_why(record) {
        result.insert("why".into(), Value::String(why));
    }
    result.insert(
        "warning".into(),
        usage_warning(record, material).map_or(Value::Null, Value::String),
    );
    if let Some(validation) = validation {
        if let Some(check) = validation
            .get("status")
            .and_then(Value::as_str)
            .and_then(SourceCheck::from_status_code)
        {
            result.insert("source_check".into(), Value::String(check.as_str().into()));
        }
        if let Some(date) = validation.get("audit_date") {
            result.insert("source_checked_on".into(), date.clone());
        }
    }
    if let Some(finding) = source_finding(validation) {
        result.insert("source_note".into(), Value::String(finding.to_owned()));
    }
    if let Some(audit) = audit {
        result.insert("audit".into(), Value::String(audit.to_owned()));
    }
    Value::Object(result)
}

/// Effective `check` fields of a record, stored or derived.
#[derive(Clone, Debug, PartialEq)]
pub(super) struct Headline {
    pub(super) reproduction: Reproduction,
    pub(super) compared: Compared,
    pub(super) why: Option<String>,
    pub(super) warning: Option<String>,
    pub(super) source_check: Option<SourceCheck>,
    pub(super) source_note: Option<String>,
}

impl Default for Headline {
    fn default() -> Self {
        Self {
            reproduction: Reproduction::Untested,
            compared: Compared::Nothing,
            why: None,
            warning: None,
            source_check: None,
            source_note: None,
        }
    }
}

/// Effective `check` fields of a record.
///
/// A stored `check` value wins; otherwise the value is derived from inline
/// evidence exactly as [`check_block`] would. No assessment means
/// [`Reproduction::Untested`] and [`Compared::Nothing`].
pub(super) fn headline(record: &Value, material: Option<&Value>) -> Headline {
    let check = record.get("check").and_then(Value::as_object);
    let stored = |key: &str| {
        check
            .and_then(|check| check.get(key))
            .and_then(Value::as_str)
            .filter(|value| !value.is_empty())
    };
    let validation = validation(record);
    Headline {
        reproduction: stored("reproduction")
            .and_then(Reproduction::from_check)
            .unwrap_or_else(|| reproduction_headline(record)),
        compared: stored("compared")
            .and_then(Compared::from_check)
            .unwrap_or_else(|| compared_scope(record)),
        why: match check.and_then(|check| check.get("why")) {
            Some(value) => value.as_str().map(str::to_owned),
            None => reproduction_why(record),
        },
        warning: match check.and_then(|check| check.get("warning")) {
            Some(value) => value.as_str().map(str::to_owned),
            None => usage_warning(record, material),
        },
        source_check: stored("source_check")
            .and_then(SourceCheck::from_check)
            .or_else(|| {
                validation
                    .and_then(|validation| validation.get("status"))
                    .and_then(Value::as_str)
                    .and_then(SourceCheck::from_status_code)
            }),
        source_note: stored("source_note")
            .or_else(|| source_finding(validation))
            .map(str::to_owned),
    }
}

/// Validate the structure of an optional record `check` block.
#[allow(clippy::too_many_lines)]
pub(super) fn validate_check(check: &Value, location: &str) -> Result<(), EosmatError> {
    let check = check
        .as_object()
        .ok_or_else(|| invalid_document(format!("{location}.check must be an object")))?;
    for (old, new) in RENAMED_CHECK_KEYS {
        if check.contains_key(*old) {
            let codes = if *old == "source" {
                format!(
                    " ({})",
                    SourceCheck::ALL
                        .map(|value| format!("{} -> {}", value.status_code(), value.as_str()))
                        .join(", ")
                )
            } else {
                String::new()
            };
            return Err(invalid_document(format!(
                "{location}.check.{old} is no longer read; EOSMAT v5 names this field check.{new}{codes}"
            )));
        }
    }
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
    let reproduction = match check.get("reproduction") {
        Some(value) => Some(value.as_str().and_then(Reproduction::from_check).ok_or_else(
            || {
                invalid_document(format!(
                    "{location}.check.reproduction must be one of reproduced, close, disagrees, untested"
                ))
            },
        )?),
        None => None,
    };
    let compared = match check.get("compared") {
        Some(value) => Some(value.as_str().and_then(Compared::from_check).ok_or_else(|| {
            invalid_document(format!(
                "{location}.check.compared must be one of published_fit, peritheos_refit, source_equation, none"
            ))
        })?),
        None => None,
    };
    if let (Some(reproduction), Some(compared)) = (reproduction, compared) {
        // Only an assessed comparison can give a result; `untested` has none
        // or only a rebuilt source equation.
        let allowed = if reproduction == Reproduction::Untested {
            matches!(compared, Compared::Nothing | Compared::SourceEquation)
        } else {
            compared != Compared::Nothing
        };
        if !allowed {
            return Err(invalid_document(format!(
                "{location}.check.compared {:?} cannot accompany reproduction {:?}",
                compared.as_str(),
                reproduction.as_str()
            )));
        }
    }
    if let Some(value) = check.get("source_check") {
        if value.as_str().and_then(SourceCheck::from_check).is_none() {
            return Err(invalid_document(format!(
                "{location}.check.source_check must be one of verified, pending, unresolved, deferred"
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
    for key in ["why", "source_note", "source_checked_on", "audit"] {
        if let Some(value) = check.get(key) {
            if value.as_str().is_none_or(str::is_empty) {
                return Err(invalid_document(format!(
                    "{location}.check.{key} must be a non-empty string"
                )));
            }
        }
    }
    if let Some(why) = check.get("why").and_then(Value::as_str) {
        if why.chars().count() > WHY_MAX_LENGTH {
            return Err(invalid_document(format!(
                "{location}.check.why must be at most {WHY_MAX_LENGTH} characters"
            )));
        }
        if !matches!(
            reproduction.unwrap_or(Reproduction::Close),
            Reproduction::Close | Reproduction::Disagrees
        ) {
            return Err(invalid_document(format!(
                "{location}.check.why is only given for reproduction close or disagrees"
            )));
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
    let status = assessment_status(record);
    if let (Some(stored), Some(status)) = (check.get("reproduction"), &status) {
        let expected = match status.as_str() {
            "parity"
            | "similar"
            | "not_reproduced"
            | "not_refittable"
            | "source_reconstruction"
            | "not_assessed" => Some(Reproduction::from_assessment_status(status).as_str()),
            _ => None,
        };
        if stored.as_str() != expected {
            return Err(invalid_document(format!(
                "{location}.check.reproduction {stored} contradicts assessment status {status:?}"
            )));
        }
    }
    if let (Some(stored), Some(status)) = (check.get("compared"), &status) {
        let expected = compared_scope(record);
        if stored.as_str() != Some(expected.as_str()) {
            let scope = assessment(record)
                .and_then(|assessment| assessment.get("scope"))
                .unwrap_or(&Value::Null);
            return Err(invalid_document(format!(
                "{location}.check.compared {stored} contradicts assessment status {status:?} and scope {scope} (expected {:?})",
                expected.as_str()
            )));
        }
    }
    if let (Some(stored), Some(status)) = (
        check.get("source_check"),
        validation(record).and_then(|validation| validation.get("status")),
    ) {
        let expected = status
            .as_str()
            .and_then(SourceCheck::from_status_code)
            .map(SourceCheck::as_str);
        if stored.as_str() != expected {
            return Err(invalid_document(format!(
                "{location}.check.source_check {stored} contradicts scientific_validation.status {status}"
            )));
        }
    }
    Ok(())
}

/// Validate a record's `superseded_record_text` audit block.
///
/// The block keeps the original value of every free-text field replaced by a
/// text revision (see [`SUPERSEDED_TEXT_FIELDS`]) and the ISO date
/// `superseded_on`. It never changes how the record is read.
pub(super) fn validate_superseded_record_text(
    value: &Value,
    location: &str,
) -> Result<(), EosmatError> {
    let place = format!("{location}.superseded_record_text");
    let object = value
        .as_object()
        .ok_or_else(|| invalid_document(format!("{place} must be an object")))?;
    let unknown: BTreeSet<&str> = object
        .keys()
        .map(String::as_str)
        .filter(|key| *key != "superseded_on" && !SUPERSEDED_TEXT_FIELDS.contains(key))
        .collect();
    if !unknown.is_empty() {
        return Err(invalid_document(format!(
            "{place} has unknown fields {unknown:?}"
        )));
    }
    if !SUPERSEDED_TEXT_FIELDS
        .iter()
        .any(|key| object.contains_key(*key))
    {
        return Err(invalid_document(format!(
            "{place} names no superseded field"
        )));
    }
    if !object
        .get("superseded_on")
        .and_then(Value::as_str)
        .is_some_and(is_iso_date)
    {
        return Err(invalid_document(format!(
            "{place}.superseded_on must be an ISO date (YYYY-MM-DD)"
        )));
    }
    for key in ["label", "notes", "thermal_notes"] {
        if object
            .get(key)
            .is_some_and(|value| !value.is_null() && !value.is_string())
        {
            return Err(invalid_document(format!(
                "{place}.{key} must be a string or null"
            )));
        }
    }
    if object.get("validity_notes").is_some_and(|value| {
        !value.is_null()
            && !value
                .as_array()
                .is_some_and(|notes| notes.iter().all(Value::is_string))
    }) {
        return Err(invalid_document(format!(
            "{place}.validity_notes must be an array of strings or null"
        )));
    }
    Ok(())
}

/// Whether `text` is a calendar date written `YYYY-MM-DD`.
fn is_iso_date(text: &str) -> bool {
    let bytes = text.as_bytes();
    if bytes.len() != 10
        || bytes[4] != b'-'
        || bytes[7] != b'-'
        || !bytes
            .iter()
            .enumerate()
            .all(|(index, byte)| index == 4 || index == 7 || byte.is_ascii_digit())
    {
        return false;
    }
    let number = |range: std::ops::Range<usize>| text[range].parse::<u32>().unwrap_or(0);
    let (year, month, day) = (number(0..4), number(5..7), number(8..10));
    let leap = (year % 4 == 0 && year % 100 != 0) || year % 400 == 0;
    let days = match month {
        1 | 3 | 5 | 7 | 8 | 10 | 12 => 31,
        4 | 6 | 9 | 11 => 30,
        2 if leap => 29,
        2 => 28,
        _ => return false,
    };
    year >= 1 && (1..=days).contains(&day)
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
        if let Some(block) = entry.get("superseded_record_text") {
            validate_superseded_record_text(block, &format!("audit sidecar entry {identifier:?}"))?;
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
