# Rust third-party notices

Peritheos code remains MIT. Each dependency retains the license expression below;
OR alternatives remain alternatives. Supplied texts include all available options,
COPYRIGHT/NOTICE files, nested notices and additional source-header attributions.

The notice set conservatively covers all normal/build dependencies, including
procedural macros and conditional packages not selected by the five wheel targets.
Build-tool inclusion does not claim that tool binaries are bundled. `manifest.json`
records exact crate archive checksums, target roles and every notice-file hash.

Official wheels use Rust 1.99.1 (58cf8f9c35b355199727627d07ae34a4e2846d67).
Its complete standard-library copyright document is in
[COPYRIGHT-library.html](rust-stdlib/COPYRIGHT-library.html). This includes
standard-library dependencies and build inputs without relicensing Peritheos.
Audited Rust source: https://github.com/rust-lang/rust/tree/1.99.1
Dependency source archives are available at the exact crates.io versions below;
Peritheos does not modify those dependency or standard-library sources.
Other-toolchain redistributors must refresh matching notices before distribution.

| Crate | Version | Declared license | Wheel-target selection |
| --- | --- | --- | --- |
| [autocfg](https://crates.io/crates/autocfg/1.5.1) | 1.5.1 | Apache-2.0 OR MIT | all five targets |
| [cfg-if](https://crates.io/crates/cfg-if/1.0.4) | 1.0.4 | MIT OR Apache-2.0 | all five targets |
| [crossbeam-deque](https://crates.io/crates/crossbeam-deque/0.8.7) | 0.8.7 | MIT OR Apache-2.0 | all five targets |
| [crossbeam-epoch](https://crates.io/crates/crossbeam-epoch/0.9.20) | 0.9.20 | MIT OR Apache-2.0 | all five targets |
| [crossbeam-utils](https://crates.io/crates/crossbeam-utils/0.8.22) | 0.8.22 | MIT OR Apache-2.0 | all five targets |
| [either](https://crates.io/crates/either/1.18.0) | 1.18.0 | MIT OR Apache-2.0 | all five targets |
| [getrandom](https://crates.io/crates/getrandom/0.2.17) | 0.2.17 | MIT OR Apache-2.0 | all five targets |
| [heck](https://crates.io/crates/heck/0.5.0) | 0.5.0 | MIT OR Apache-2.0 | all five targets |
| [itoa](https://crates.io/crates/itoa/1.0.18) | 1.0.18 | MIT OR Apache-2.0 | all five targets |
| [libc](https://crates.io/crates/libc/0.2.189) | 0.2.189 | MIT OR Apache-2.0 | all five targets |
| [matrixmultiply](https://crates.io/crates/matrixmultiply/0.3.11) | 0.3.11 | MIT/Apache-2.0 | all five targets |
| [memchr](https://crates.io/crates/memchr/2.8.3) | 2.8.3 | Unlicense OR MIT | all five targets |
| [ndarray](https://crates.io/crates/ndarray/0.17.2) | 0.17.2 | MIT OR Apache-2.0 | all five targets |
| [num-complex](https://crates.io/crates/num-complex/0.4.6) | 0.4.6 | MIT OR Apache-2.0 | all five targets |
| [num-integer](https://crates.io/crates/num-integer/0.1.47) | 0.1.47 | MIT OR Apache-2.0 | all five targets |
| [num-traits](https://crates.io/crates/num-traits/0.2.19) | 0.2.19 | MIT OR Apache-2.0 | all five targets |
| [numpy](https://crates.io/crates/numpy/0.29.0) | 0.29.0 | BSD-2-Clause | all five targets |
| [once_cell](https://crates.io/crates/once_cell/1.21.4) | 1.21.4 | MIT OR Apache-2.0 | all five targets |
| [portable-atomic](https://crates.io/crates/portable-atomic/1.15.0) | 1.15.0 | Apache-2.0 OR MIT | conditional extra; none of the five targets |
| [portable-atomic-util](https://crates.io/crates/portable-atomic-util/0.2.7) | 0.2.7 | Apache-2.0 OR MIT | conditional extra; none of the five targets |
| [ppv-lite86](https://crates.io/crates/ppv-lite86/0.2.21) | 0.2.21 | MIT OR Apache-2.0 | all five targets |
| [proc-macro2](https://crates.io/crates/proc-macro2/1.0.107) | 1.0.107 | MIT OR Apache-2.0 | all five targets |
| [pyo3](https://crates.io/crates/pyo3/0.29.2) | 0.29.2 | MIT OR Apache-2.0 | all five targets |
| [pyo3-build-config](https://crates.io/crates/pyo3-build-config/0.29.2) | 0.29.2 | MIT OR Apache-2.0 | all five targets |
| [pyo3-ffi](https://crates.io/crates/pyo3-ffi/0.29.2) | 0.29.2 | MIT OR Apache-2.0 | all five targets |
| [pyo3-macros](https://crates.io/crates/pyo3-macros/0.29.2) | 0.29.2 | MIT OR Apache-2.0 | all five targets |
| [pyo3-macros-backend](https://crates.io/crates/pyo3-macros-backend/0.29.2) | 0.29.2 | MIT OR Apache-2.0 | all five targets |
| [quote](https://crates.io/crates/quote/1.0.47) | 1.0.47 | MIT OR Apache-2.0 | all five targets |
| [rand](https://crates.io/crates/rand/0.8.8) | 0.8.8 | MIT OR Apache-2.0 | all five targets |
| [rand_chacha](https://crates.io/crates/rand_chacha/0.3.1) | 0.3.1 | MIT OR Apache-2.0 | all five targets |
| [rand_core](https://crates.io/crates/rand_core/0.6.4) | 0.6.4 | MIT OR Apache-2.0 | all five targets |
| [rawpointer](https://crates.io/crates/rawpointer/0.2.1) | 0.2.1 | MIT/Apache-2.0 | all five targets |
| [rayon](https://crates.io/crates/rayon/1.12.0) | 1.12.0 | MIT OR Apache-2.0 | all five targets |
| [rayon-core](https://crates.io/crates/rayon-core/1.13.0) | 1.13.0 | MIT OR Apache-2.0 | all five targets |
| [rustc-hash](https://crates.io/crates/rustc-hash/2.1.3) | 2.1.3 | Apache-2.0 OR MIT | all five targets |
| [serde](https://crates.io/crates/serde/1.0.229) | 1.0.229 | MIT OR Apache-2.0 | all five targets |
| [serde_core](https://crates.io/crates/serde_core/1.0.229) | 1.0.229 | MIT OR Apache-2.0 | all five targets |
| [serde_derive](https://crates.io/crates/serde_derive/1.0.229) | 1.0.229 | MIT OR Apache-2.0 | all five targets |
| [serde_json](https://crates.io/crates/serde_json/1.0.151) | 1.0.151 | MIT OR Apache-2.0 | all five targets |
| [syn](https://crates.io/crates/syn/2.0.119) | 2.0.119 | MIT OR Apache-2.0 | all five targets |
| [syn](https://crates.io/crates/syn/3.0.4) | 3.0.4 | MIT OR Apache-2.0 | all five targets |
| [target-lexicon](https://crates.io/crates/target-lexicon/0.13.5) | 0.13.5 | Apache-2.0 WITH LLVM-exception | all five targets |
| [unicode-ident](https://crates.io/crates/unicode-ident/1.0.24) | 1.0.24 | (MIT OR Apache-2.0) AND Unicode-3.0 | all five targets |
| [wasi](https://crates.io/crates/wasi/0.11.1+wasi-snapshot-preview1) | 0.11.1+wasi-snapshot-preview1 | Apache-2.0 WITH LLVM-exception OR Apache-2.0 OR MIT | conditional extra; none of the five targets |
| [zerocopy](https://crates.io/crates/zerocopy/0.8.56) | 0.8.56 | BSD-2-Clause OR Apache-2.0 OR MIT | all five targets |
| [zerocopy-derive](https://crates.io/crates/zerocopy-derive/0.8.56) | 0.8.56 | BSD-2-Clause OR Apache-2.0 OR MIT | conditional extra; none of the five targets |
| [zmij](https://crates.io/crates/zmij/1.0.23) | 1.0.23 | MIT | all five targets |

Source archive/license files are unchanged copies verified against original crate
archive checksums in Cargo.lock. Additional source comment notices are grouped with
their exact source locations. Do not omit these notices when redistributing the
corresponding compiled code. This is not a scientific-data license grant.
