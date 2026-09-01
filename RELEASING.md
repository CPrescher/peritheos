# Releasing Peritheos

Python releases are built by GitHub Actions and published to PyPI through
Trusted Publishing. No PyPI API token is stored in GitHub. The public
`peritheos-core` and `peritheos-fit` Rust crates are published separately to
crates.io; the private `peritheos-python` extension crate is never published.

## One-time setup

1. Create a GitHub environment named `pypi` in the `CPrescher/peritheos`
   repository. Adding required reviewers is recommended.
2. On PyPI, add a Trusted Publisher for these exact values:
   - PyPI project name: `peritheos`
   - GitHub owner: `CPrescher`
   - GitHub repository: `peritheos`
   - Workflow filename: `publish.yml`
   - Environment name: `pypi`
3. If the PyPI project does not exist yet, configure a pending publisher with
   the same values. The first successful publication will create the project.
4. Before the first Rust release, sign in to crates.io with GitHub, verify the
   account email, create a scoped API token, and authenticate `cargo`. The
   first publication claims the crate names and must be performed manually.
   Trusted Publishing can be configured for later Rust releases after each
   crate exists.

## Publish a release

1. Move the `Unreleased` changelog entries into a dated version section.
2. Update `peritheos.__version__`, the Cargo workspace version,
   `peritheos-fit`'s `peritheos-core` dependency, `Cargo.lock`, and
   `CITATION.cff` to the same version and release date.
3. Run the Python, Rust, documentation, wheel, source-distribution, and crate
   archive release gates. Merge the release change into `main` and ensure every
   CI job succeeds.
4. Confirm the hosted documentation reflects the release commit.
5. Create a matching annotated tag locally without pushing it yet, for example:

   ```bash
   git tag -a v0.6.0 -m "Release 0.6.0"
   ```

6. From the clean tagged commit, publish the public Rust crates in dependency
   order. Wait until `peritheos-core` is visible in the crates.io index before
   publishing `peritheos-fit`:

   ```bash
   cargo publish -p peritheos-core --locked
   cargo publish -p peritheos-fit --locked
   ```

   crates.io versions are permanent and cannot be overwritten. Do not publish
   `peritheos-python`; its manifest intentionally has `publish = false`.
7. Push the annotated tag to start the Python publication:

   ```bash
   git push origin v0.6.0
   ```

The `Publish to PyPI` workflow checks that the tag and package version match,
builds and smoke-tests the complete platform-specific native-wheel matrix,
builds and validates the source distribution, and then publishes the exact
uploaded artifacts from a separate OIDC-enabled job. The full test suite is a
required main-branch CI gate and must be green before the tag is pushed.
After PyPI accepts the artifacts, the workflow creates a GitHub Release with
the matching version section from `CHANGELOG.md` and attaches the same source
and wheel distributions. Re-running the workflow synchronizes an existing
release description with the changelog.

If an automated system pushes the tag and GitHub suppresses the push trigger,
dispatch the same workflow manually. It still checks out and validates the
existing release tag:

```bash
gh workflow run publish.yml --ref main -f version=0.6.0
```
