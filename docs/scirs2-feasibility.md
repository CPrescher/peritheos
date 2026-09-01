# SciRS2 feasibility audit

This audit evaluates whether SciRS2 can replace the SciPy functionality used
by Peritheos. It is intentionally narrower than asking whether SciRS2 is a
general scientific-computing library: every decision below is tied to an
existing Peritheos code path and compatibility requirement.

## Audit scope and decision

The source-level audit was performed on 2026-08-30 against the independently
versioned SciRS2 0.6.5 crates published on crates.io:
`scirs2-optimize`, `scirs2-integrate`, `scirs2-sparse`, `scirs2-linalg`, and
`scirs2-stats`. All are Apache-2.0 and identify the
[SciRS2 repository](https://github.com/cool-japan/scirs) as their source. The
published manifests do not declare a Rust version, so their effective MSRV is
not a documented compatibility guarantee.

**Decision:** do not add SciRS2 to the production dependency graph at this
stage. Version 0.6.5 is not a drop-in foundation for Peritheos fitting, sparse
covariance, root finding, or production quadrature. Selected future SciRS2
components may be reconsidered behind Peritheos-owned interfaces after they
pass the same numerical and platform tests as any other candidate.

This is not a permanent rejection of the project. The audit is versioned, and
must be repeated before adopting a later release.

## Exact Peritheos requirements

The current implementation uses SciPy and NumPy for the following operations:

| Peritheos path | Required behavior |
|---|---|
| EOS inversion | A caller-supplied physical bracket, reliable scalar root convergence, and a post-solve residual check |
| Debye and Sokolova thermal terms | Accurate adaptive one-dimensional quadrature over smooth and limiting-regime integrands |
| Ordinary fitting | Nonlinear least squares, parameter bounds, column/Jacobian scaling, finite-difference Jacobians, robust loss, and complete diagnostics |
| Errors-in-variables fitting | All ordinary-fit features plus hundreds or thousands of positive latent state variables and a block-local Jacobian sparsity pattern |
| Fit covariance | Dense SVD pseudoinverse for ordinary fits; sparse Schur-complement profiling for latent-variable fits |
| Observation covariance | Batched Cholesky factorization and triangular solves |
| Uncertainty | Symmetric eigenvalues, normal quantiles, multivariate-normal sampling, and output covariance |

These requirements must work in combination. A separate bounded solver, a
separate robust solver, and a separate sparse solver do not satisfy the fitting
contract unless their capabilities compose without changing the objective or
reported statistics.

## Findings

### Nonlinear least squares: not suitable as the fitting backend

`scirs2-optimize` 0.6.5 contains useful individual APIs, but not the combined
solver required by Peritheos:

- `least_squares::least_squares` offers LM, TRF, and Dogbox methods and dense
  optional Jacobians. Its options have evaluation and convergence tolerances
  but no bounds, robust loss, residual scale, Jacobian sparsity, or equivalent
  of SciPy's `x_scale="jac"`.
- `least_squares::bounded::bounded_least_squares` is a separate dense path. It
  has bounds, but no robust loss or sparse-Jacobian input. Its no-bounds branch
  contains a source comment that calling the regular solver is still to be
  wired up.
- `least_squares::robust::robust_least_squares` is another separate dense path
  based on IRLS or steepest descent. It has no bounds or sparse-Jacobian input
  and forms dense normal equations.
- `least_squares::sparse::sparse_least_squares` initially accepts a dense
  Jacobian, converts it to an internal sparse representation, then recomputes a
  dense Jacobian on every sparse Gauss--Newton iteration. Its step divides the
  gradient by the diagonal of `J^T J`; it does not solve the coupled sparse
  least-squares system. The initially supplied sparse matrix is not used by
  the iterative solve.
- The common result type does not provide SciPy-equivalent `cost`,
  `optimality`, or shaped final-Jacobian semantics across all solver paths.
  Several paths leave the status code at its default value.

Peritheos currently passes bounds, `jac_sparsity`, `x_scale="jac"`, robust
`loss`, `f_scale`, and `max_nfev` together. Errors-in-variables fits also rely
on the final sparse Jacobian to profile latent observations out of the
parameter covariance. None of the SciRS2 0.6.5 entry points provides that
combined contract.

### Sparse linear algebra: unsuitable for large latent systems

The public `scirs2_sparse::linalg::spsolve` implementation converts its CSR
matrix to a dense matrix and then performs Gaussian elimination.
`sparse_direct_solve` delegates to it, including when symmetry or positive
definiteness hints are supplied. `sparse_lstsq` forms normal equations and then
uses the same dense-converting solve.

This behavior is unacceptable for the latent-information system in
Peritheos's covariance calculation: memory changes from proportional to the
number of nonzeros to quadratic in the number of latent variables. It would
turn a presently sparse workload into the exact scaling failure the Rust
migration is meant to avoid.

SciRS2 also exposes iterative sparse solvers, but substituting one is not
mechanical. Peritheos needs deterministic convergence criteria, handling of
rank loss and ill conditioning, multiple right-hand sides, and verified
covariance accuracy. Those properties require a dedicated prototype and
stress suite before any dependency decision.

### Scalar root finding: insufficient physical control

The optimizer's scalar root method is reached through the general vector-root
API. It starts from one point, searches around `x0` using a fixed initial
offset, expands at most ten times, and then combines Newton steps with a
bisection fallback. The public call does not accept Peritheos's already-known
physical volume or temperature bracket.

That is not equivalent to Peritheos's branch contract. In particular, EOS
inversion must deliberately select compression or the first expansion branch,
stay in the positive physical domain, and verify the pressure residual.
A small Peritheos-owned bracketed Brent-style solver is therefore the lower
risk option and avoids exposing optimizer types in the core crate.

### Quadrature: available but not yet production-grade for this use

`scirs2_integrate::quad` defaults to composite Simpson estimates on 10 and 20
subintervals and recursively bisects intervals. It is not a QUADPACK-style
Gauss--Kronrod implementation. Its optional fixed-Simpson path always reports
an approximate `1e-8` error and successful convergence after 1,001
evaluations, rather than deriving that error from the integrand.

The API may be adequate for some smooth integrals, but Peritheos evaluates a
Debye integral across small-argument, ordinary, and large-argument regimes and
the Sokolova volume integral over model-dependent functions. Adoption would
require comparison with existing SciPy values, high-precision reference
values, reversed integration limits, domain edges, evaluation counts, and
failure behavior. Given the narrow one-dimensional need, a tested specialized
implementation or a focused quadrature crate is preferable to adding the
broader SciRS2 dependency now.

### Dense linear algebra and statistics: capable pieces, weak dependency case

`scirs2-linalg` exposes Cholesky, SVD, pseudoinverse, and eigenvalue routines,
and `scirs2-stats` exposes a normal distribution with a percent-point
function. These cover names that Peritheos needs, but API presence is not
numerical parity evidence. They also bring much broader dependency surfaces
than the isolated operations justify.

Dense matrix operations should be selected together with the fitting solver
after rank-deficient covariance and cross-platform wheel prototypes. A normal
quantile and multivariate-normal sampling can be implemented with smaller,
focused dependencies once the RNG reproducibility policy is fixed.

## Capability matrix

| Capability | SciRS2 0.6.5 assessment | Migration decision |
|---|---|---|
| Analytical EOS formulas | Not needed | Implement directly in `peritheos-core` |
| Physical bracketed roots | Partial, wrong public control surface | Implement behind a private core interface |
| One-dimensional quadrature | Candidate only | Benchmark and validate focused alternatives first |
| Bounded robust least squares | Separate non-composable APIs | Do not use as the fitting backend |
| Sparse EIV least squares | Dense Jacobian regeneration and diagonal step | Reject for production fitting |
| Sparse direct solve | Densifies CSR | Reject for latent covariance |
| Dense Cholesky/SVD/eigen | API exists | Re-evaluate with the native fitting prototype |
| Normal quantile | API exists | Prefer a focused implementation/dependency |
| Random sampling | Broad alternatives exist | Select only after defining deterministic backend behavior |

## Architecture consequence

The migration keeps numerical ownership inside Peritheos:

1. `peritheos-core` will contain formulas, state validation, explicit branch
   selection, a private bracketed root solver, and only the quadrature needed
   by thermal models.
2. `peritheos-fit` will define its own residual, bounds, robust-loss,
   convergence, diagnostics, and covariance interfaces. Dense and sparse
   matrix backends remain replaceable implementation details.
3. The initial Python binding will port deterministic EOS evaluation before
   fitting. SciPy remains the fitting oracle and fallback until native solver
   parity is independently demonstrated.
4. No Rust public API will contain SciRS2 types. A later audited SciRS2 release
   can therefore be adopted for an individual component without an API break.

## Re-evaluation gates

A future SciRS2 release may be adopted only when a pinned prototype shows:

- a single solver path combining bounds, all five losses, Jacobian scaling,
  sparse finite-difference coloring or an equivalent analytic strategy, and
  complete convergence diagnostics;
- no dense materialization proportional to the square of the latent-variable
  count;
- covariance agreement on full-rank, ill-conditioned, and rank-deficient
  fixtures;
- explicit physical brackets for scalar roots;
- thermal-integral agreement throughout Peritheos's tested domain;
- successful Linux, macOS, and Windows builds at the workspace MSRV; and
- acceptable compile time, wheel size, runtime, and license inventory.

The relevant published source is available through the exact-version
[optimizer](https://docs.rs/crate/scirs2-optimize/0.6.5/source/),
[integrator](https://docs.rs/crate/scirs2-integrate/0.6.5/source/), and
[sparse](https://docs.rs/crate/scirs2-sparse/0.6.5/source/) crate archives.
