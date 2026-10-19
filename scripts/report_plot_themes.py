"""Light/dark report images using Material for MkDocs' palette selectors."""

from __future__ import annotations

import hashlib
import re
from pathlib import PurePosixPath

# Match the documentation's slate background; retain distinct source/refit hues.
# Transform paint only: paths, coordinates, labels, IDs and scientific values
# remain byte-for-byte identical to the light SVG.
DARK_COLORS = {
    "#ffffff": "#1e2129",  # figure and axes background
    "#000000": "#e6e8ef",  # ticks and spines
    "#245b8a": "#82baff",  # source curve/coefficients
    "#ae4c14": "#ffad73",  # reconstruction curve/coefficients
    "#4e5660": "#bdc5d3",  # observations and zero line
    "#aab2bc": "#8590a3",  # coordinate error bars
    "#aab4c0": "#8590a3",  # coefficient connectors
    "#c3cad2": "#687486",  # coefficient zero line
    "#e7eaee": "#3b414f",  # grid
}


def dark_path(path: str) -> str:
    source = PurePosixPath(path)
    return str(source.with_name(source.stem + "-dark" + source.suffix))


def themed_images(alt: str, path: str) -> list[str]:
    """Both images follow the active documentation palette, including its toggle."""
    return [
        f"![{alt}]({path}#only-light)",
        f"![{alt}]({dark_path(path)}#only-dark)",
    ]


def dark_svg(light: bytes) -> bytes:
    """Derive a deterministic dark SVG without reevaluating or rerendering data."""
    source_hash = hashlib.sha256(light).hexdigest()
    payload = light.decode("utf-8")
    # The light renderer's body hashes apply only to its original bytes.
    payload = re.sub(r"<!-- peritheos-(?:comparison|curves) .*? -->\n", "", payload)

    def paint(match: re.Match[str]) -> str:
        return re.sub(
            r"#[0-9a-fA-F]{6}\b",
            lambda color: DARK_COLORS.get(color[0].lower(), color[0]),
            match[0],
        )

    payload = re.sub(r'\b(?:style|fill|stroke)="[^"]*"', paint, payload)
    # Matplotlib omits the default black fill on text, so set it explicitly.
    payload = payload.replace(
        "</defs>",
        '<style type="text/css">text { fill: #e6e8ef; }</style>\n </defs>',
        1,
    )
    marker = f"<!-- peritheos-dark source={source_hash} -->\n"
    return payload.replace("?>\n", "?>\n" + marker, 1).encode("utf-8")
