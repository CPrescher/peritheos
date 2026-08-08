# Releasing Peritheos

Releases are built by GitHub Actions and published to PyPI through Trusted
Publishing. No PyPI API token is stored in GitHub.

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

## Publish a release

1. Move the `Unreleased` changelog entries into a dated version section.
2. Update `peritheos.__version__` and `CITATION.cff` to the same version and
   release date.
3. Merge the release change into `main` and ensure every CI job succeeds.
4. Confirm the hosted documentation reflects the release commit.
5. Create and push a matching annotated tag, for example:

   ```bash
   git tag -a v0.1.0 -m "Release 0.1.0"
   git push origin v0.1.0
   ```

The `Publish to PyPI` workflow checks that the tag and package version match,
runs the tests, builds and validates the wheel and source distribution, and
then publishes the exact uploaded artifacts from a separate OIDC-enabled job.
After PyPI accepts the artifacts, the workflow creates a GitHub Release with
generated notes and attaches the same source and wheel distributions.

If an automated system pushes the tag and GitHub suppresses the push trigger,
dispatch the same workflow manually. It still checks out and validates the
existing release tag:

```bash
gh workflow run publish.yml --ref main -f version=0.1.0
```
