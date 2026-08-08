import re
from importlib.metadata import version
from pathlib import Path

import peritheos


def test_runtime_version_matches_package_metadata():
    assert peritheos.__version__ == version("peritheos")


def test_citation_version_matches_runtime_version():
    citation = Path(__file__).parents[1].joinpath("CITATION.cff").read_text()
    match = re.search(r"^version:\s*(\S+)\s*$", citation, flags=re.MULTILINE)

    assert match is not None
    assert match.group(1) == peritheos.__version__
