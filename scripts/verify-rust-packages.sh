#!/usr/bin/env bash

set -euo pipefail

repository_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
verification_dir="$(mktemp -d "${TMPDIR:-/tmp}/peritheos-packages.XXXXXX")"
trap 'rm -rf -- "${verification_dir}"' EXIT

cd "${repository_dir}"
workspace_version="$(awk -F '"' '/^version = "/ { print $2; exit }' Cargo.toml)"
if [[ -z "${workspace_version}" ]]; then
  echo "Could not determine the Cargo workspace version." >&2
  exit 1
fi

# Cargo requires a registry dependency to exist before packaging a dependent
# crate. Prepare peritheos-fit with a local resolution patch, while leaving the
# generated archive's normalized dependency as a versioned registry dependency.
cargo package -p peritheos-core --locked
cargo package -p peritheos-fit --locked --no-verify \
  --config "patch.crates-io.peritheos-core.path=\"${repository_dir}/crates/peritheos-core\""

core_archive="${repository_dir}/target/package/peritheos-core-${workspace_version}.crate"
fit_archive="${repository_dir}/target/package/peritheos-fit-${workspace_version}.crate"
tar -xf "${core_archive}" -C "${verification_dir}"
tar -xf "${fit_archive}" -C "${verification_dir}"

packaged_core="${verification_dir}/peritheos-core-${workspace_version}"
packaged_fit="${verification_dir}/peritheos-fit-${workspace_version}"
cargo test --manifest-path "${packaged_core}/Cargo.toml" --locked
cargo test --manifest-path "${packaged_fit}/Cargo.toml" --locked \
  --config "patch.crates-io.peritheos-core.path=\"${packaged_core}\""
