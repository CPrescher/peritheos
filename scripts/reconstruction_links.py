"""Stable public and repository-local destinations for bundled EOS assessments."""

from urllib.parse import quote

DOCS_URL = "https://peritheos.readthedocs.io/en/latest/"


def report_links(material: str, record: str) -> dict[str, str]:
    """Keep the same record anchor in hosted and locally built documentation."""
    page = f"reconstruction-reports/{quote(material, safe='')}/"
    anchor = "record-" + quote(record, safe="")
    return {
        "report_url": f"{DOCS_URL}{page}#{anchor}",
        "report_path": f"docs/reconstruction-reports/{material}.md#{anchor}",
    }
