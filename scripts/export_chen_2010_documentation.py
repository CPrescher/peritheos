"""Export a local, network-independent comparison page for Studio.

Requires Markdown (included in the optional MkDocs documentation toolchain).
"""

import argparse
import hashlib
import json
import shutil
from pathlib import Path

from scripts.generate_paper_investigation_ledger import chen_comparison_section
from scripts.reproduce_chen_2010_argon import ROOT


def export(target: Path):
    import markdown

    target.mkdir(parents=True, exist_ok=True)
    files = [
        "docs/data/chen-2010-argon-comparison.png",
        "docs/data/chen-2010-argon-comparison.json",
        "docs/data/chen-2010-argon-refit.json",
        "docs/data/chen-2010-argon-refit.png",
        "docs/data/chen-2010-argon-reproduction.json",
        "docs/literature-reproductions/chen-2010-argon.md",
    ]
    hashes = {}
    for filename in files:
        source = ROOT / filename
        shutil.copyfile(source, target / source.name)
        hashes[filename] = hashlib.sha256(source.read_bytes()).hexdigest()
    text = chen_comparison_section()
    text = text.replace("(data/", "(").replace("(literature-reproductions/", "(")
    text += """
## Reference parameters

These are zero-pressure **coefficients of the stored equations**, not measured ambient solid properties.
The original-constant record derives them from the reported constraints at 2 GPa.

| Quantity | Printed-constant reconstruction: DO NOT USE | Peritheos refit |
|---|---:|---:|
| V0, Å³ / four-atom cell | 158.44581 | 157.05440 |
| K0, GPa | 2.57967 | 4.39886 |
| K0′ | 9.08160 | 4.62176 |
| Extrapolated rho0, g/cm³ | 1.67465 | 1.68948 |

## How the refit was obtained

All 80 distinct black-square positions in Figure 5 enter once. An errors-in-variables
least-squares objective permits density adjustment and uses the graphical pressure
and density halfwidths as relative coordinate weights. Original observations are
preserved. The paper's fitted curve, extrapolated zero-pressure density and printed
elastic constants do not enter this regression. The uncertainty convention and
cross-row correlations of the integrated densities are unknown, so no parameter
covariance, confidence interval or statistical standard errors are asserted.

Different weighting choices change predicted volumes by up to 0.62% in the observed
range; omitting individual pressure bins changes them by up to 1.67%. Zero-pressure
parameters are less stable. This is a separate approximation to the digitized data,
not a reconstruction of the complete original experiment or an independently
validated elastic equation of state.

Read the [complete source audit and fitting details](chen-2010-argon.md),
[original reconstruction diagnostics](chen-2010-argon-reproduction.json), and
[source publication](https://doi.org/10.1103/PhysRevB.81.144110).
"""
    body = markdown.markdown(text, extensions=["tables", "fenced_code"])
    html = (
        """<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<title>Chen 2010 argon: original constants and independent refit</title>
<style>
body{font:17px/1.65 system-ui,sans-serif;color:#263342;background:#f7f9fc;margin:0}
main{max-width:1120px;margin:32px auto;padding:30px;background:white;border:1px solid #dbe2ec;border-radius:12px}
h1,h2{line-height:1.25}h2{margin-top:2em}h2:first-child{margin-top:0}
img{width:100%;height:auto}table{width:100%;border-collapse:collapse;font-size:.92em;display:block;overflow:auto}
td,th{text-align:left;border-bottom:1px solid #dbe2ec;padding:10px}th{background:#eef2f8}
a{color:#195ca5}code{font-size:.85em;overflow-wrap:anywhere}strong{color:#792c25}
@media(max-width:600px){main{margin:0;padding:18px;border:0}body{font-size:15px}}
</style></head><body><main>"""
        + body
        + "</main></body></html>\n"
    )
    (target / "index.html").write_text(html, encoding="utf-8")
    (target / "source-manifest.json").write_text(
        json.dumps(
            {
                "source_files_sha256": hashes,
                "generator": "scripts/export_chen_2010_documentation.py",
            },
            indent=2,
        )
        + "\n"
    )


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    export(parser.parse_args().output)
