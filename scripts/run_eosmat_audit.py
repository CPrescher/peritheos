"""Run a historical, hash-bound v3 audit against native library files.

Usage: ``python -m scripts.run_eosmat_audit AUDIT [--check | --register]``.
The original scientific script and its recorded fingerprints stay unchanged.
Only its material-file access uses a temporary v3 compatibility export.
"""

from __future__ import annotations

import argparse
import copy
import importlib
import json
import sys
import tempfile
from pathlib import Path
from types import ModuleType
from unittest.mock import patch

import peritheos
from scripts.eosmat_io import read_material, write_material

AUDITS = {
    "refit_hirose_fei_gold_speziale": "MATERIAL",
    "refit_miozzi_2020_speziale": "IRON_PATH",
    "refit_miozzi_2020_tange": "IRON_PATH",
}


def historical_reconstruction_view(document: dict) -> dict:
    """Exclude the later ledger summary from a historical full-record fingerprint.

    This projection is used only by the explicit historical-audit runner. It
    retains original evidence, coefficients and datasets. Ordinary library and
    compatibility exports always include the new assessment.
    """
    result = copy.deepcopy(document)
    for record in result.get("eos_records", []):
        reproduction = record.get("reproduction", {})
        assessment = reproduction.get("assessment", {})
        if (
            assessment.get("report") == "docs/data/primary-eos-refits.json"
            and "source_status" in assessment
        ):
            del reproduction["assessment"]
            if not reproduction:
                del record["reproduction"]
    return result


def run_audit(
    module: ModuleType, material_attribute: str, arguments: list[str]
) -> None:
    """Keep native files native, committing a material update only after success."""
    path = getattr(module, material_attribute)
    original = path.read_bytes()
    document = read_material(path)
    argv = sys.argv
    get_document = peritheos.get_material_document

    def historical_document(*args, **kwargs):
        return historical_reconstruction_view(get_document(*args, **kwargs))

    with tempfile.TemporaryDirectory(prefix="peritheos-eosmat-audit-") as directory:
        temporary = Path(directory) / path.name
        temporary.write_text(
            json.dumps(document, indent=1, ensure_ascii=False) + "\n", encoding="utf-8"
        )
        try:
            setattr(module, material_attribute, temporary)
            sys.argv = [module.__name__, *arguments]
            # The archived calibration hash predates the additive ledger summary.
            # Some scripts import this lookup locally, others at module scope.
            with (
                patch.object(peritheos, "get_material_document", historical_document),
                patch.object(
                    module, "get_material_document", historical_document, create=True
                ),
            ):
                module.main()
            updated = read_material(temporary)
        finally:
            setattr(module, material_attribute, path)
            sys.argv = argv
        if updated != document:
            if path.read_bytes() != original:
                raise RuntimeError(f"Material changed during audit: {path}")
            write_material(path, updated)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("audit", choices=sorted(AUDITS))
    args, arguments = parser.parse_known_args()
    module = importlib.import_module("scripts." + args.audit)
    run_audit(module, AUDITS[args.audit], arguments)


if __name__ == "__main__":
    main()
