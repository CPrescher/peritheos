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


def test_readme_describes_the_catalog_without_hard_coded_counts():
    """Catalog totals are generated; the README links to them (AGENTS.md)."""
    readme = Path(__file__).parents[1].joinpath("README.md").read_text()
    documents = peritheos.list_material_documents()
    record_count = sum(
        len(peritheos.get_material_document(identifier)["eos_records"])
        for identifier in documents
    )
    words = set(readme.replace("/", " ").replace("-", " ").split())
    assert not {str(len(documents)), str(record_count)} & words
    assert "material documents and EOS records" in " ".join(readme.split())
    assert "docs/reconstruction-reports/index.md" in readme
