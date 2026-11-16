"""Show bundled materials and EOS records; create, convert, validate and format EOSMAT files."""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path
from typing import Any, Sequence

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


def _show(args: argparse.Namespace) -> str:
    from peritheos.catalog import get_eos_record, get_material
    from peritheos.errors import MaterialLookupError

    try:
        material = get_material(args.material)
    except MaterialLookupError:
        if args.record is not None:
            raise
        # A bare record identifier is accepted as well.
        try:
            record = get_eos_record(args.material)
        except MaterialLookupError:
            pass
        else:
            return _render(record, args)
        raise
    if args.record is None:
        return _render(material, args)
    return _render(material.get_eos_record(args.record), args)


def _render(item: Any, args: argparse.Namespace) -> str:
    overview = item.overview(include_audit=args.include_audit)
    if args.json:
        return json.dumps(overview, ensure_ascii=False, indent=2, allow_nan=False)
    return item.summary()


def _list(args: argparse.Namespace) -> str:
    from peritheos.catalog import list_materials, search_materials

    materials = search_materials(args.query) if args.query else list_materials()
    rows = []
    for material in materials:
        default = material.default_record()
        count = len(material.eos_records)
        rows.append(
            (
                material.identifier,
                material.formula,
                material.name or "",
                " ".join(material.phase.split()),
                f"{count} record{'s' if count != 1 else ''}",
                "-"
                if default is None
                else " ".join(default.name.split())
                + ("" if default.is_default else " (first record)"),
            )
        )
    if not rows:
        return f"No bundled material matches {args.query!r}."
    header = ("identifier", "formula", "name", "phase", "records", "default")
    widths = [
        min(40, max(len(row[column]) for row in (header, *rows))) for column in range(5)
    ]
    return "\n".join(
        "  ".join(
            [*(f"{cell[:width]:<{width}}" for cell, width in zip(row, widths)), row[5]]
        ).rstrip()
        for row in (header, *rows)
    )


def _catalog_command(args: argparse.Namespace, parser: argparse.ArgumentParser) -> int:
    from peritheos.errors import PeritheosError

    try:
        text = _show(args) if args.command == "show" else _list(args)
    except (PeritheosError, KeyError) as error:
        message = error.args[0] if error.args else str(error)
        print(f"{parser.prog} {args.command}: error: {message}", file=sys.stderr)
        return 2
    print(text)
    return 0


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
    create.add_argument("--format-version", type=int, choices=[3, 4, 5], default=5)
    convert = commands.add_parser(
        "convert", help="Migrate a complete EOSMAT file between v3, v4 and v5"
    )
    convert.add_argument("source", type=Path)
    convert.add_argument(
        "path", type=Path, help="New destination; existing files are kept"
    )
    convert.add_argument("--format-version", type=int, choices=[3, 4, 5], default=5)
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
    show = commands.add_parser(
        "show", help="Show a bundled material, or one of its EOS records"
    )
    show.add_argument(
        "material", help="Material identifier, alias, name or formula (e.g. MgO)"
    )
    show.add_argument("record", nargs="?", help="Optional EOS record identifier")
    show.add_argument(
        "--json", action="store_true", help="Print the overview as JSON instead"
    )
    show.add_argument(
        "--include-audit",
        action="store_true",
        help="Add the pressure scale from the audit sidecar (JSON and text)",
    )
    listing = commands.add_parser(
        "list", help="List bundled materials, optionally matching a search text"
    )
    listing.add_argument("query", nargs="?", help="Search text, e.g. a name or DOI")
    args = parser.parse_args(argv)
    if args.command in {"show", "list"}:
        return _catalog_command(args, parser)
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
                # Formatting keeps a record file's layout; validation also
                # requires referenced audit sidecars to exist and merge cleanly.
                document = load_eosmat(path, include_audit=args.command == "validate")
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
            # Conversion always writes a complete file: sidecar evidence that a
            # record file references is merged in.
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
        with args.path.open("x", encoding="utf-8", newline="\n") as stream:
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
