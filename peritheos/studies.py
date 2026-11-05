"""Typed datasets from supporting studies without assigning an executable EOS."""

from __future__ import annotations

import copy
import json
from importlib import resources
from typing import Any

from peritheos.datasets import Dataset
from peritheos.errors import DatasetLookupError

__all__ = ["list_studies", "get_study_document", "get_study_dataset"]


def _documents() -> list[dict[str, Any]]:
    directory = resources.files("peritheos.data").joinpath("studies")
    if not directory.is_dir():
        return []
    documents = [
        json.loads(path.read_text(encoding="utf-8"))
        for path in sorted(directory.iterdir(), key=lambda p: p.name)
        if path.is_file() and path.name.endswith(".json")
    ]
    return [
        d
        for d in documents
        if d.get("datasets") and all("columns" in v for v in d["datasets"])
    ]


def list_studies() -> tuple[str, ...]:
    """Return supporting-study identifiers, separately from material/EOS discovery."""
    return tuple(sorted(document["identifier"] for document in _documents()))


def get_study_document(identifier: str) -> dict[str, Any]:
    """Return an independent metadata copy; observation values remain optional."""
    from peritheos._optional_data import study_dataset_declarations

    for original in _documents():
        if original["identifier"] == identifier:
            document = copy.deepcopy(original)
            optional = {
                d["identifier"]: d for d in study_dataset_declarations(identifier)
            }
            document["datasets"] = [
                optional.get(d["identifier"], d) for d in document.get("datasets", [])
            ]
            return document
    raise DatasetLookupError(
        f"Unknown supporting study {identifier!r}",
        code="dataset.study_not_found",
        operation="load_study",
        context={"study": identifier},
    )


def get_study_dataset(study: str, identifier: str) -> Dataset:
    """Load a study table with the ordinary Dataset checksum and version checks."""
    document = get_study_document(study)
    for declaration in document.get("datasets", []):
        if declaration["identifier"] == identifier:
            return Dataset.from_mapping(declaration)
    raise DatasetLookupError(
        f"Study {study!r} has no dataset {identifier!r}",
        operation="load_dataset",
        context={"study": study, "identifier": identifier},
    )
