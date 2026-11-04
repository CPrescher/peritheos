"""Load preserved numerical curve readouts without distributing source rasters."""

import hashlib
import json
from pathlib import Path

READOUTS = (
    Path(__file__).resolve().parents[1] / "docs/data/fei-source-curve-readouts.json"
)
READOUTS_SHA256 = "58c7a6835a4ecde80ef0ba4c1daf1a69bdb3a77aae7f479cce5fba4e4cc55935"


def curve_readouts(figure):
    content = READOUTS.read_bytes()
    if hashlib.sha256(content).hexdigest() != READOUTS_SHA256:
        raise ValueError(
            "Fei source-curve pixel readouts differ from the preserved audit"
        )
    return json.loads(content)[figure]
