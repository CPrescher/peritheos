# Peritheos project goals and working guidance

## What this project is for

Peritheos is a reusable scientific library for equations of state,
thermodynamics, fitting and uncertainty. Its Python and Rust interfaces and its
material database should support researchers doing their own calculations as
well as applications consuming the library. Peritheos Studio and Dioptas are
downstream applications, not the definition of the library's interface.

Make the scientific work understandable and reproducible for human users.
Machine readability is necessary, but a file or API is not successful merely
because an AI can navigate it. Prefer clear names, predictable organization,
explicit conventions and complete working examples over layers of abstraction.

## Human-readable EOSMAT

- Support both small standalone user files and complete curated database records
  in the same native format. Basic authoring must not require database
  registration, audit boilerplate or a generated companion representation.
- A reader should quickly find the material/phase, EOS model, coefficients,
  units, reference state, source and limits. Keep the usual path directly
  readable; advanced capabilities should be optional.
- Optimize for clarity, not just fewer bytes or lines. Removing repetition is
  useful only when it does not replace readable content with confusing
  indirection. Inline metadata remains a first-class choice.
- Keep scientific evidence available and traceable. A short summary cannot
  replace the original coefficients, observations, uncertainty semantics,
  calibration, fitting protocol or important qualifications.
- Preserve old files and downstream interoperability through supported readers
  and explicit compatibility exports. Document capability boundaries; never
  silently reinterpret an unsupported model or file feature.

## Scientific fidelity

- Separate measured observations, digitized observations, calculated checkpoints,
  published parameterizations and independent refits. Do not invent missing
  observations, errors, cell contents, calibration inputs or fitting choices.
- Preserve published values and source data when investigating a discrepancy.
  Keep any diagnostic reconstruction or independently requested refit separately
  identifiable. Numerical agreement does not authorize replacing a source record.
- State units, volume normalization, reference-state conventions, model laws,
  fixed parameters and uncertainty/confidence/covariance assumptions explicitly.
- Use the [EOS label definitions](docs/eos-labels.md) consistently.
- Distinguish structural validation, source verification, numerical parameter
  parity, author-workflow reproduction, physical validation and recommendation
  for use. Reconstruction findings need their actual comparison scope, basis
  and unresolved limitations; avoid a single misleading quality score.
- For catalog curation, follow [Adding materials and EOS](docs/adding-materials-and-eos.md)
  and the [EOS curation guidance](.agents/skills/peritheos-eos-curation/SKILL.md).
  Format or presentation changes must not alter scientific findings.

## Documentation and verification

- Treat documentation as part of implementation. Update the authoring guide,
  relevant format/API reference and complete examples when behavior changes.
  Label legacy examples and compatibility requirements clearly.
- End-user instructions must work with the installed package. Identify
  repository-only maintenance tools as such.
- Test the scientific or compatibility property that could break: source-value
  and array preservation, loading and round trips, independent numerical cases,
  Python/Rust consistency, invalid inputs and executable documentation as relevant.
  Do not equate passing tests with experimental validation.
- Use [Contributing](CONTRIBUTING.md) and [API stability](docs/api-stability.md)
  for maintained development checks and compatibility policy.

## Scope and shared-checkout care

Work on the requested repository and task. Do not let a library task drift into
Studio work or an unrelated scientific investigation. Preserve existing local
changes and untracked evidence. Do not regenerate catalogs, rewrite material
files or update scientific classifications merely as incidental cleanup.

Keep commits focused; stage only owned changes, including individual hunks when
files contain other work. Report what was tested and what was committed,
published or still unverified. Do not push, release or deploy unless requested.

These are enduring goals, not a backlog or authorization to continue unrelated
roadmap work. A user's explicit instructions determine the current task.
