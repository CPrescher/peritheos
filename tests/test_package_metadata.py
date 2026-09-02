import re
from importlib import resources
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


def test_package_declares_inline_typing_support():
    assert resources.files("peritheos").joinpath("py.typed").is_file()


def test_readme_catalog_count_matches_bundled_documents():
    readme = Path(__file__).parents[1].joinpath("README.md").read_text()
    record_count = sum(
        len(peritheos.get_material_document(identifier)["eos_records"])
        for identifier in peritheos.list_material_documents()
    )

    assert (
        f"{len(peritheos.list_material_documents())}-material/{record_count}-record"
        in readme
    )
