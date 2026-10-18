"""Migrate existing EOSMAT files to readable, inline v4 without scientific edits.

Run ``python -m scripts.migrate_eosmat_v4 [paths ...]``; no paths selects the
bundled library. ``--check`` reports pending migrations without writing files.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Any

from peritheos.eosmat import convert_eosmat
from scripts.eosmat_io import render_material

ROOT = Path(__file__).resolve().parents[1]


def scientific_content(document: dict[str, Any]) -> str:
    """Compare complete legacy content, allowing only empty thermal annotations.

    V4 puts errors/fixed flags on their coefficient. An empty optional container
    has no annotation to retain. All nonempty containers and array order count.
    """

    legacy = convert_eosmat(document, version=3)
    for record in legacy["eos_records"]:
        thermal = record.get("thermal", {})
        for key in ("parameter_errors", "fixed_parameters"):
            if key in thermal and not thermal[key]:
                thermal.pop(key)
    # JSON distinguishes 1 from 1.0 and preserves signed zero; Python equality
    # would accept those accidental changes. Object order alone is immaterial.
    return json.dumps(legacy, sort_keys=True, ensure_ascii=False, allow_nan=False)


def migrate_file(path: Path, *, check: bool = False) -> bool:
    """Return whether migration is needed; verify complete content before writing."""
    original = path.read_bytes()
    document = json.loads(original)
    rendered = render_material(document).encode("utf-8")
    if scientific_content(document) != scientific_content(json.loads(rendered)):
        raise ValueError(f"Migration changed scientific content: {path}")
    if original == rendered:
        return False
    if not check:
        if path.read_bytes() != original:
            raise RuntimeError(f"File changed while migrating: {path}")
        path.write_bytes(rendered)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    paths = args.paths or [ROOT / "peritheos/data/materials"]
    files = sorted(
        {
            file
            for path in paths
            for file in (path.glob("*.eosmat") if path.is_dir() else [path])
        }
    )
    changed = [path for path in files if migrate_file(path, check=args.check)]
    if args.check and changed:
        for path in changed:
            print(f"Needs v4 migration: {path}")
        return 1
    print(
        f"{'Checked' if args.check else 'Migrated'} {len(files)} files; {len(changed)} changed."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
