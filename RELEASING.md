# Releasing peritheos

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

1. Update `peritheos.__version__` in `peritheos/__init__.py`.
2. Merge the version change into `main` and ensure the CI workflow succeeds.
3. Create and push a matching annotated tag, for example:

   ```bash
   git tag -a v0.1.0 -m "Release 0.1.0"
   git push origin v0.1.0
   ```

The `Publish to PyPI` workflow checks that the tag and package version match,
runs the tests, builds and validates the wheel and source distribution, and
then publishes the exact uploaded artifacts from a separate OIDC-enabled job.
