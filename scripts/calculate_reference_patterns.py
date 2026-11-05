"""Calculate reference powder arrays from the separately licensed input models.

Repository maintenance only: uv run --with gemmi==0.7.4 python
scripts/calculate_reference_patterns.py --check

The default prints a report; --write changes only the ten peak arrays. See
peritheos/data/calculated-pattern-inputs.json for assumptions and source limits.
"""

from __future__ import annotations

import argparse
import hashlib
import itertools
import json
import math
from collections import Counter
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
INPUTS = ROOT / "peritheos/data/calculated-pattern-inputs.json"


def expanded_model(entry, root=ROOT):
    import gemmi

    if gemmi.__version__ != "0.7.4":
        raise ValueError("Reproduction requires gemmi==0.7.4")
    tolerance = 5e-6
    if "resource" in entry:
        path = root / "peritheos/data" / entry["resource"]
        if hashlib.sha256(path.read_bytes()).hexdigest() != entry["sha256"]:
            raise ValueError("Pattern source checksum changed")
        structure = gemmi.read_small_structure(str(path))
        # Explicit CIF operations preserve the nonstandard majorite origin.
        operations = [gemmi.Op(op) for op in structure.symops]
        sites = [(s.element.name, s.fract.tolist(), s.occ) for s in structure.sites]
    else:
        operations = list(gemmi.SpaceGroup(entry["space_group"]).operations())
        sites = [
            (s["element"], [s[k] for k in ("x", "y", "z")], s["occupancy"])
            for s in entry["atom_sites"]
        ]
    if not operations:
        raise ValueError("Missing explicit source symmetry")
    atoms = []
    for element, xyz, occupancy in sites:
        # CIF decimal thirds otherwise split exact special-position orbits.
        xyz = [
            round(x * 12) / 12 if abs(x - round(x * 12) / 12) <= tolerance else x
            for x in xyz
        ]
        orbit = {
            tuple(round(x % 1, 7) % 1 for x in op.apply_to_xyz(xyz))
            for op in operations
        }
        for target, fraction in (
            entry["substitutions"].get(element, {element: 1}).items()
        ):
            atoms.extend(
                (target, position, occupancy * fraction) for position in sorted(orbit)
            )
    contents = Counter()
    for element, _, occupancy in atoms:
        contents[element] += occupancy
    return atoms, dict(sorted(contents.items()))


def calculate(entry):
    import gemmi

    atoms, contents = expanded_model(entry)
    axes = entry["lattice"]
    cell = gemmi.UnitCell(*(axes[k] for k in ("a", "b", "c", "alpha", "beta", "gamma")))
    limits = [math.ceil(axes[k] / 0.8) for k in ("a", "b", "c")]
    indices, distances = [], []
    for hkl in itertools.product(*(range(-n, n + 1) for n in limits)):
        if hkl == (0, 0, 0):
            continue
        d = cell.calculate_d(hkl)
        if 0.8 <= d <= 20:
            indices.append(hkl)
            distances.append(d)
    hkl = np.asarray(indices)
    d = np.asarray(distances)
    factors = np.zeros(len(indices), dtype=complex)
    # Full reciprocal sphere includes multiplicity and Friedel partners without
    # relying on a space-group setting lookup or a single representative's F.
    for element in sorted(contents):
        sites = [a for a in atoms if a[0] == element]
        xyz = np.asarray([a[1] for a in sites])
        occupancy = np.asarray([a[2] for a in sites])
        scattering = gemmi.Element("H" if element == "D" else element).it92
        form = np.asarray([scattering.calculate_sf(float(1 / (4 * x * x))) for x in d])
        factors += form * np.sum(
            np.exp(2j * np.pi * np.einsum("ij,kj->ik", hkl, xyz)) * occupancy, axis=1
        )
    intensity = abs(factors) ** 2
    groups = {}
    threshold = float(intensity.max()) * 1e-12
    for index, distance, value in zip(indices, d, intensity):
        if value > threshold:
            groups.setdefault(round(float(distance), 10), []).append(
                (index, float(value))
            )
    maximum = max(sum(v for _, v in family) for family in groups.values())
    peaks = []
    for distance, family in sorted(groups.items(), reverse=True):
        # Prefer an all-nonnegative representative, then a deterministic order.
        representative = min(
            (i for i, _ in family), key=lambda i: (sum(x < 0 for x in i), i)
        )
        strength = round(100 * sum(v for _, v in family) / maximum, 8)
        if strength > 0:
            peaks.append([*representative, distance, strength])
    return peaks, contents


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", action="store_true")
    mode.add_argument("--write", action="store_true")
    args = parser.parse_args()
    inputs = json.loads(INPUTS.read_bytes())
    for material, entry in inputs["materials"].items():
        peaks, contents = calculate(entry)
        path = ROOT / "peritheos/data/materials" / f"{material}.eosmat"
        document = json.loads(path.read_bytes())
        if args.check and document["peaks"] != peaks:
            raise SystemExit(f"Stale calculated pattern: {material}")
        if args.write:
            from scripts.eosmat_io import render_material

            document["peaks"] = peaks
            path.write_text(render_material(document), encoding="utf-8")
        print(f"{material}: {len(peaks)} calculated lines; cell contents {contents}")


if __name__ == "__main__":
    main()
