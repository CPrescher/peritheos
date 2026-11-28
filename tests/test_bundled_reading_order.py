"""A reader of a bundled record file meets the EOS records before the datasets."""

from __future__ import annotations

import json
from pathlib import Path

MATERIALS = Path(__file__).resolve().parents[1] / "peritheos/data/materials"


def test_eos_records_come_before_datasets_in_every_bundled_file():
    late = []
    for path in sorted(MATERIALS.glob("*.eosmat")):
        keys = list(json.loads(path.read_bytes()))
        if "datasets" in keys and keys.index("datasets") < keys.index("eos_records"):
            late.append(path.stem)
    assert late == []
