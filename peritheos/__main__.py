"""Create, convert, validate and format installed EOSMAT material files."""

from __future__ import annotations

import argparse
import json
import math
from pathlib import Path
from typing import Sequence

from peritheos._eosmat_format import format_eosmat
from peritheos.authoring import create_eosmat
from peritheos.eosmat import convert_eosmat, load_eosmat


def _parameters(tokens: Sequence[str]) -> dict[str, float]:
    result = {}
    for token in tokens:
        key, separator, raw = token.partition("=")
        if not separator or not key:
            raise ValueError(f"Expected PARAMETER=VALUE, got {token!r}")
        if key in result:
            raise ValueError(f"Duplicate parameter {key!r}")
        value = float(raw)
        if not math.isfinite(value):
            raise ValueError(f"Parameter {key!r} must be finite")
        result[key] = value
    return result


def main(argv: Sequence[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m peritheos", description=__doc__)
    commands = parser.add_subparsers(dest="command", required=True)
    create = commands.add_parser(
        "create", help="Create a native standalone EOSMAT file"
    )
    create.add_argument(
        "path", type=Path, help="New .eosmat file; existing files are kept"
    )
    create.add_argument(
        "--id", dest="identifier", help="Material ID; default: file stem"
    )
    create.add_argument("--name", help="Human-readable name; default: material ID")
    create.add_argument("--formula", required=True, help="Material chemical formula")
    create.add_argument(
        "--model", required=True, help="Exact model name, e.g. BM3 or Vinet"
    )
    create.add_argument(
        "--parameters",
        nargs="+",
        required=True,
        help="PARAMETER=VALUE pairs; pressures in GPa, V0 in conventional-cell Å³",
    )
    create.add_argument(
        "--temperature", type=float, required=True, help="Reference T in K"
    )
    create.add_argument(
        "--reference", help="Source citation; default: user-supplied parameters"
    )
    create.add_argument(
        "--errors", nargs="+", default=[], help="Reported PARAMETER=ERROR pairs"
    )
    create.add_argument("--fixed", nargs="+", default=[], help="Parameters held fixed")
    create.add_argument("--z", type=float, help="Formula units per conventional cell")
    create.add_argument("--thermal-model", help="Optional exact thermal model name")
    create.add_argument(
        "--thermal-parameters", nargs="+", help="Thermal PARAMETER=VALUE pairs"
    )
    create.add_argument("--format-version", type=int, choices=[3, 4], default=4)
    convert = commands.add_parser(
        "convert", help="Migrate a complete EOSMAT file between v3 and v4"
    )
    convert.add_argument("source", type=Path)
    convert.add_argument(
        "path", type=Path, help="New destination; existing files are kept"
    )
    convert.add_argument("--format-version", type=int, choices=[3, 4], default=4)
    layout = convert.add_mutually_exclusive_group()
    layout.add_argument(
        "--compact",
        action="store_true",
        help="Share identical metadata inside native v4",
    )
    layout.add_argument(
        "--inline", action="store_true", help="Expand metadata for older v4 readers"
    )
    validate = commands.add_parser(
        "validate",
        help="Check file structure, local references and model declarations; not scientific validation",
    )
    validate.add_argument("paths", nargs="+", type=Path)
    formatter = commands.add_parser(
        "format", help="Format files in place without changing JSON values"
    )
    formatter.add_argument("paths", nargs="+", type=Path)
    formatter.add_argument(
        "--check", action="store_true", help="Read-only; exit 1 if formatting is needed"
    )
    args = parser.parse_args(argv)
    if args.command == "create" and bool(args.thermal_model) != bool(
        args.thermal_parameters
    ):
        parser.error(
            "--thermal-model and --thermal-parameters must be supplied together"
        )
    try:
        if args.command in {"validate", "format"}:
            prepared = []
            for path in args.paths:
                before = path.read_bytes()
                document = load_eosmat(path)
                if args.command == "validate":
                    print(
                        f"Valid EOSMAT v{document['format_version']}: {path} (structural validation only)"
                    )
                    continue
                rendered = format_eosmat(document).encode("utf-8")
                if json.loads(rendered) != document:
                    raise ValueError(f"Formatting changed JSON contents: {path}")
                prepared.append((path, before, rendered))
            if args.command == "validate":
                return 0
            changed = [
                (path, before, rendered)
                for path, before, rendered in prepared
                if before != rendered
            ]
            if args.check:
                for path, _, _ in changed:
                    print(f"Needs formatting: {path}")
                return int(bool(changed))
            # Validate every requested file and check ownership before any write.
            if any(path.read_bytes() != before for path, before, _ in prepared):
                raise ValueError(
                    "A source file changed while formatting; no files written"
                )
            for path, _, rendered in changed:
                path.write_bytes(rendered)
            print(f"Formatted {len(prepared)} EOSMAT files; {len(changed)} changed.")
            return 0
        if args.command == "convert":
            document = convert_eosmat(
                load_eosmat(args.source),
                version=args.format_version,
                compact=args.compact,
                inline=args.inline,
            )
        else:
            thermal = None
            if args.thermal_model:
                thermal = {
                    "model": args.thermal_model,
                    "parameters": _parameters(args.thermal_parameters),
                }
            document = create_eosmat(
                args.identifier if args.identifier is not None else args.path.stem,
                args.formula,
                args.model,
                _parameters(args.parameters),
                temperature_ref=args.temperature,
                name=args.name,
                reference=args.reference,
                parameter_errors=_parameters(args.errors),
                fixed_parameters=args.fixed,
                formula_units_per_cell=args.z,
                thermal=thermal,
                format_version=args.format_version,
            )
        payload = format_eosmat(document)
        with args.path.open("x", encoding="utf-8") as stream:
            stream.write(payload)
    except (ValueError, TypeError, OSError) as error:
        cause: BaseException = error
        while cause.__cause__ is not None:
            cause = cause.__cause__
        parser.error(str(error) if cause is error else f"{error}: {cause}")
    print(
        f"Wrote {args.path} (EOSMAT v{args.format_version}; scientific status preserved)."
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
