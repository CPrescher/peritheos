# Rust third-party notice audit

Reviewed 4 November 2026 against Cargo.lock and Rust 1.99.1, commit
`58cf8f9c35b355199727627d07ae34a4e2846d67`.

The Python wheel statically incorporates Rust code into `_rust`; the Python
and Rust implementation licenses remain MIT. Dependency licenses retain their
own terms. Cargo's license-policy check establishes accepted expressions; it
does not by itself put the necessary copyright/license notices in packages.

## Included evidence and notices

The installed `peritheos/data/third_party/rust/` directory contains:

- `NOTICE.md`: readable dependency index and source-code availability links.
- `manifest.json`: all 47 locked normal/build packages, archive hashes, notice
  hashes and dependency roles for five supported wheel targets.
- `licenses/`: unmodified upstream LICENSE, COPYING, COPYRIGHT and NOTICE files,
  including nested license files. OR alternatives retain their original meaning.
- `source-notices/`: additional original source-comment attributions with file
  locations, including serde_json's lexical-derived Alexander Huszagh notice.
- `rust-stdlib/COPYRIGHT-library.html` and `stdlib.json`: the complete, unmodified
  official standard-library copyright document from the audited Rustup compiler,
  with version, commit, hash and source URL. The document includes standard-library
  dependencies and build inputs; including their notices does not claim those
  tools or the Rust compiler itself are distributed.

All these files are included in wheels and source distributions and declared
as `License-File` metadata. They are also readable without finding the wheel's
metadata directory:

```python
from importlib.resources import files

print(files("peritheos").joinpath("data/third_party/rust/NOTICE.md").read_text())
```

Original crate archives are checked against Cargo.lock. Every copied license
and scanned source file must match the corresponding original archive bytes.
Git preserves notice bytes without text conversion. Only the Cargo.lock
fingerprint normalizes CRLF to LF; notice and source payload hashes do not.

## Target selection and limits

Cargo metadata selects 43 third-party packages for each of:

- `aarch64-apple-darwin`
- `aarch64-unknown-linux-gnu`
- `x86_64-apple-darwin`
- `x86_64-pc-windows-msvc`
- `x86_64-unknown-linux-gnu`

The additional all-target conditional packages are portable-atomic,
portable-atomic-util, wasi and zerocopy-derive. Their notices are supplied
conservatively too. Target roles distinguish libraries from host build tools
and procedural macros. A selected package does not imply all its machine code
survives linking; build-tool notices are not assertions that tool binaries ship.
Dev-only dependencies are excluded. Generated-code inputs retain their notices.

The earlier complete 35-wheel archive inventory found one native extension per
wheel, with no additional native payloads. The local macOS extension dynamically
uses platform libSystem/libiconv. NumPy, SciPy and Python are separately installed,
not bundled by dependency declarations. The installed smoke check rejects extra
native libraries so that a future repaired wheel/DLL payload requires a new audit.
This review does not establish scientific-data redistribution rights.

## Reproducing and maintaining the audit

Repository-only checks:

```bash
python scripts/generate_rust_notices.py --check
```

To update after dependency changes, regenerate and review the resulting files:

```bash
python scripts/generate_rust_notices.py
```

When changing the official release compiler, install its Rust documentation,
refresh the standard-library notices, and update `RELEASE_RUST_TOOLCHAIN` in
the publication workflow to the same audited version:

```bash
rustup component add rust-docs
python scripts/generate_rust_notices.py --refresh-stdlib
python scripts/generate_rust_notices.py --check --verify-toolchain
```

Official wheel builds pin that compiler and check the bundle before Linux builds
or alongside host builds. Source/MSRV builds continue to support Rust 1.83;
pinning the publication compiler does not raise the library's minimum version.
Redistributors rebuilding with another compiler must supply/review its matching
standard-library notices. Toolchains without a library-only copyright document
require a separate review/import of their complete upstream copyright material;
the generator fails instead of substituting the 1.99.1 document. The committed
Rust 1.99.1 document is not a claim of notice coverage for arbitrary other
toolchain versions.

The full dependency-audit job reruns the notice check. Routine tests guard graph
roles, source/archive tampering and committed payload hashes. Artifact checks
verify package copies against license metadata and reject unexpected native files.
