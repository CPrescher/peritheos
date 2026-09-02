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

package_options=(--locked)
if [[ "${PERITHEOS_ALLOW_DIRTY:-0}" == "1" ]]; then
  package_options+=(--allow-dirty)
fi

cargo package -p peritheos "${package_options[@]}"

archive="${repository_dir}/target/package/peritheos-${workspace_version}.crate"
tar -xf "${archive}" -C "${verification_dir}"

packaged_crate="${verification_dir}/peritheos-${workspace_version}"
cargo test --manifest-path "${packaged_crate}/Cargo.toml" --locked
