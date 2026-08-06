from importlib.metadata import version

import peritheos


def test_runtime_version_matches_package_metadata():
    assert peritheos.__version__ == version("peritheos")
