"""Locate the optional, version-matched dataset distribution without network I/O."""

from __future__ import annotations

import importlib.util
import json
from importlib import import_module, resources

from peritheos.errors import DatasetError


def data_package():
    if importlib.util.find_spec("peritheos_datasets") is None:
        return None
    peritheos_datasets = import_module("peritheos_datasets")

    from peritheos import __version__

    if peritheos_datasets.__version__ != __version__:
        raise DatasetError(
            f"peritheos-datasets {peritheos_datasets.__version__} does not match "
            f"Peritheos {__version__}. Install matching versions with "
            '`python -m pip install --upgrade "peritheos[data]"`.',
            code="dataset.package_version_mismatch",
            operation="load_dataset",
        )
    return resources.files("peritheos_datasets")


def dataset_declarations(material):
    try:
        package = data_package()
    except DatasetError as error:
        if error.code == "dataset.package_version_mismatch":
            # EOS evaluation must remain independent from an incompatible
            # optional installation. Data access reports the mismatch itself.
            return []
        raise
    if package is None:
        return []
    path = package.joinpath("declarations", material + ".json")
    return json.loads(path.read_text(encoding="utf-8")) if path.is_file() else []


def unavailable_dataset(material, identifier, availability):
    package = data_package()
    url = availability["url"]
    if package is None:
        message = (
            f"Observation dataset {identifier!r} is not installed. "
            'Install selected datasets with `python -m pip install "peritheos[data]"`. '
            f"Other source datasets and their terms: {url}"
        )
        code = "dataset.package_not_installed"
    else:
        message = (
            f"Dataset {identifier!r} is not included in the optional distribution. "
            f"Its source evidence and terms are available at {url}."
        )
        code = "dataset.not_distributed"
    raise DatasetError(
        message,
        code=code,
        operation="load_dataset",
        context={"material": material, "identifier": identifier, "source_url": url},
    )
