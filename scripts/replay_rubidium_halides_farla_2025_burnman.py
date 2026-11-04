#!/usr/bin/env python3
"""Execute the hash-pinned historical CsCl module with optional BurnMan 1.2.0.

Example (isolated dependencies; does not change the project environment):
uv run --no-project --python 3.11 --with burnman==1.2.0 python \
    scripts/replay_rubidium_halides_farla_2025_burnman.py \
    --source /path/to/historical/Dewaele_2020.py --output /path/to/replay.json

The historical module remains external; its source hash is checked before
execution. No network requests are made by this script.
"""

import argparse
import csv
import hashlib
import importlib.util
import json
import platform
from pathlib import Path

import burnman
import numpy as np
import scipy

ROOT = Path(__file__).resolve().parents[1]
PREFIX = "rubidium-halides-farla-2025"


def replay(source):
    facts = json.loads(
        (ROOT / "peritheos/data/datasets" / f"{PREFIX}-author-inputs.json").read_text(
            encoding="utf-8"
        )
    )
    provenance = next(s for s in facts["sources"] if s["id"] == "cscl_module")
    digest = hashlib.sha256(source.read_bytes()).hexdigest()
    if digest != provenance["sha256"]:
        raise ValueError("Historical source checksum differs from the audited module")
    if burnman.__version__ != "1.2.0":
        raise ValueError("This replay requires the historically pinned BurnMan 1.2.0")
    spec = importlib.util.spec_from_file_location("farla_historical_cscl", source)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    phase = module.CsCl()
    na = facts["cscl_gui_defaults"]["avogadro_used_mol_inverse"]
    a0 = facts["cscl_gui_defaults"]["reference_lattice_a_angstrom"]
    phase.params["V_0"] = na * a0**3 * 1e-30
    with (ROOT / "peritheos/data/datasets" / f"{PREFIX}-rbcl-observations.csv").open(
        encoding="utf-8"
    ) as stream:
        rows = [r for r in csv.DictReader(stream) if r["source_table"] != "S1"]
    pressures = [
        float(
            phase.pressure_function(
                float(r["cscl_volume_a3"]) * na * 1e-30,
                float(r["temperature_k"]),
                phase.params,
            )
            * 1e-9
        )
        for r in rows
    ]
    return {
        "audit_date": "2026-10-16",
        "source": provenance,
        "burnman_version": burnman.__version__,
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "scipy_version": scipy.__version__,
        "rows": len(rows),
        "source_sha256": digest,
        "pressures_gpa": pressures,
        "qualification": "Historical code executed on printed V,T with recovered GUI defaults; original lattice-parameter inputs and saved session are absent.",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    args.output.write_text(
        json.dumps(replay(args.source), indent=2) + "\n", encoding="utf-8"
    )


if __name__ == "__main__":
    main()
