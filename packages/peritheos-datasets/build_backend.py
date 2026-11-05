"""Repository data-package hooks; exported sdists use plain setuptools."""

import runpy
from pathlib import Path

_backend = runpy.run_path(
    str(Path(__file__).resolve().parents[2] / "peritheos_build.py")
)


def build_wheel(wheel_directory, config_settings=None, metadata_directory=None):
    return _backend["_build"](
        "build_wheel", wheel_directory, config_settings, metadata_directory, data=True
    )


def build_sdist(sdist_directory, config_settings=None):
    return _backend["_build"](
        "build_sdist", sdist_directory, config_settings, data=True
    )


def prepare_metadata_for_build_wheel(metadata_directory, config_settings=None):
    return _backend["_build"](
        "prepare_metadata_for_build_wheel",
        metadata_directory,
        config_settings,
        data=True,
    )


def get_requires_for_build_wheel(config_settings=None):
    return []


get_requires_for_build_sdist = get_requires_for_build_wheel
