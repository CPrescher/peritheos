"""Format existing EOSMAT documents without altering their scientific contents.

Run ``python -m scripts.format_eosmat [paths ...]``. With no paths, format the
bundled materials and authoring template. ``--check`` is read-only.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from peritheos._eosmat_format import format_eosmat
from peritheos.eosmat import validate_eosmat_document

ROOT = Path(__file__).resolve().parents[1]


def format_file(path: Path, *, check: bool = False) -> bool:
    """Return whether a file needs formatting, preserving every JSON value."""
    original = path.read_bytes()
    document = json.loads(original)
    validate_eosmat_document(document)
    rendered = format_eosmat(document).encode("utf-8")
    if json.loads(rendered) != document:
        raise ValueError(f"Formatting changed JSON contents: {path}")
    if original == rendered:
        return False
    if not check:
        if path.read_bytes() != original:
            raise RuntimeError(f"File changed while formatting: {path}")
        path.write_bytes(rendered)
    return True


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("paths", nargs="*", type=Path)
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    paths = args.paths or [
        ROOT / "peritheos/data/materials",
        ROOT / "docs/examples/simple-material.eosmat",
    ]
    files: set[Path] = set()
    for path in paths:
        if path.is_dir():
            files.update(path.glob("*.eosmat"))
        else:
            files.add(path)
    changed = [path for path in sorted(files) if format_file(path, check=args.check)]
    if args.check and changed:
        for path in changed:
            print(f"Needs formatting: {path}")
        return 1
    print(
        f"{'Checked' if args.check else 'Formatted'} {len(files)} EOSMAT files; "
        f"{len(changed)} {'need changes' if args.check else 'changed'}."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
