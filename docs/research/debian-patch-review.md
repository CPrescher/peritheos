# Debian packaging patch review

Reviewed 4 November 2026 against Peritheos PR 27 and Debian packaging revision
`41c5e1b75c63156eccfd6eab00510b380e2d533c`, whose changelog identifies upstream
0.17.1. The [original patch series](https://salsa.debian.org/python-team/packages/python-peritheos/-/tree/41c5e1b75c63156eccfd6eab00510b380e2d533c/debian/patches)
is by Roland Mas and remains the authority for the proposed edits.

| Patch | Upstream decision for 0.17.3 | Reason |
| --- | --- | --- |
| `0001-Update-versioned-deps.patch`: PyO3 0.29.2 → 0.28.2, rust-numpy 0.29 → 0.28 | Retain as a downstream adaptation | This selects Debian's available versions. No upstream build failure with the current versions was demonstrated, so there is no reason to downgrade the released binding stack. |
| Same patch: rand 0.8.5 → 0.10.3 | Defer upstream | Rand 0.10.3 declares Rust 1.85, while Peritheos supports Rust 1.83. This requires a deliberate compatibility-policy change rather than a patch-release dependency edit. |
| `0002-Fix-build-error.patch` | Incorporate, with attribution | An explicit `-2.0_f64` in the normal sampler makes its type clear. Current rand 0.8 builds already pass; the annotation changes no formula, dependency or intended sampling behavior. |
| `0003-Port-to-current-rand-crate.patch` | Keep with the downstream rand upgrade | `RngExt` and `random_range` are the rand 0.10 API. They cannot be applied on their own to the current rand 0.8 dependency. Revisit together with an upstream rand/MSRV migration. |
| `0004-Relax-tolerance-for-one-numerical-test.patch` | Incorporate, with an explanatory comment and attribution | Increase only the absolute saved-optimizer-SSE comparison from `1e-10` to `1e-9`. Keep the separate `1e-11` SSE-versus-current-residual check, source checksum, row selection, covariance consistency, runtime parity and scientific comparisons. |

The [rand migration guide](https://rust-random.github.io/book/update-0.10.html)
documents the trait changes. The downloaded
[rand 0.10.3 crate metadata](https://crates.io/crates/rand/0.10.3)
specifies Rust 1.85. The eight initial floating-point samples at seed 42 agreed
bit-for-bit in an isolated rand 0.8.8/0.10.3 probe; that small comparison does
not establish complete stream compatibility. Rand's
[`StdRng` documentation](https://docs.rs/rand/0.10.0/rand/rngs/struct.StdRng.html)
does not promise portable output across versions and platforms. No random
dependency is changed in 0.17.3.

## Numerical tolerance scope

The brucite report contains four diagnostic optimizer objectives between
0.234 and 2.232, with pressure and volume objectives kept separate. A `1e-9`
absolute allowance is approximately `4.5e-10` to `4.3e-9` relative to those
objectives. It accommodates numerical optimizer/backend variation; it does
not change observation errors, experimental confidence, fit acceptance or
recommendation status. The saved report and all observations remain unchanged.

The reported Debian failure was not reproduced locally. The existing tighter
test already passed on Python 3.9 with NumPy 2.0.2/SciPy 1.13.1 and Python 3.13
with NumPy 2.5.2/SciPy 1.18.1. Fresh objectives on the latter differed from the
saved values by at most `2.88e-12`; the former reproduced the saved objectives
exactly. The narrowly scoped adjustment is accepted on the maintainer's
portability report, rather than claimed as a locally reproduced failure.

The two adopted edits are checked by the brucite comparison suite on both
Python versions, the independent scientific-source fingerprint and native
tests on stable Rust and Rust 1.83. The complete release matrix must pass on
the final release revision before publication. Licensing is addressed
separately in the [installed code/data guide](https://github.com/CPrescher/peritheos/blob/main/peritheos/data/LICENSING.md);
this patch review does not establish third-party data redistribution rights.
