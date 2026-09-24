"""Extract Figure 5 coordinates from the verified publisher PDF.

Requires optional PyMuPDF. PDF is deliberately external to the repository.
Repeated rendering commands are not counted as independent observations.
"""

import argparse
import csv
import hashlib
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
PDF_SHA256 = "17a0bcd30604bd81f125894a18a45a91255b235fb80643cddd4bcbc6e5bd01f1"
X0, X25 = 374.885986328125, 513.6329956054688
Y15, Y4 = 230.6690673828125, 81.11602783203125
MASS_PER_CELL = 4 * 39.948 / 0.602214076


def coordinates(x, y):
    return (x - X0) * 25 / (X25 - X0), 1.5 + (Y15 - y) * 2.5 / (Y15 - Y4)


def extract(pdf):
    import fitz

    if hashlib.sha256(pdf.read_bytes()).hexdigest() != PDF_SHA256:
        raise ValueError("This extraction is specific to the audited publisher PDF")
    with fitz.open(pdf) as doc:
        paths = doc[2].get_drawings()
    markers = []
    for index in range(297, 445):
        path = paths[index]
        assert path["fill"] == (0, 0, 0)
        for item in path["items"]:
            assert item[0] == "re"
            rect = item[1]
            assert 4.64 < rect.width < 4.65 and 4.64 < rect.height < 4.65
            x, y = (rect.x0 + rect.x1) / 2, (rect.y0 + rect.y1) / 2
            match = next(
                (m for m in markers if abs(m[0] - x) < 0.005 and abs(m[1] - y) < 0.005),
                None,
            )
            if match is None:
                markers.append([x, y, 1, index])
            else:
                match[2] += 1
    segments = [item for path in paths[291:297] for item in path["items"]]
    header = [
        "pressure_gpa",
        "density_g_cm3",
        "volume_a3",
        "temperature_k",
        "pressure_plot_halfwidth_gpa",
        "density_plot_halfwidth_g_cm3",
        "pdf_x_pt",
        "pdf_y_pt",
        "pdf_render_count",
        "pdf_first_path_index",
    ]
    rows = []
    for x, y, count, index in sorted(markers):
        p, rho = coordinates(x, y)
        pe, re = "", ""
        for kind, a, b in segments:
            assert kind == "l"
            if abs(a.x - x) < 0.01 and abs(b.x - x) < 0.01:
                if abs((a.y + b.y) / 2 - y) < 0.01 and abs(a.y - b.y) > 4.65:
                    re = abs(a.y - b.y) / 2 * 2.5 / (Y15 - Y4)
            if abs(a.y - y) < 0.01 and abs(b.y - y) < 0.01:
                if abs((a.x + b.x) / 2 - x) < 0.01 and abs(a.x - b.x) > 0.5:
                    pe = abs(a.x - b.x) / 2 * 25 / (X25 - X0)
        rows.append([p, rho, MASS_PER_CELL / rho, 290, pe, re, x, y, count, index])
    curve = paths[445]
    assert curve["type"] == "s" and len(curve["items"]) == 261
    points = [curve["items"][0][1]] + [i[2] for i in curve["items"]]
    curve_rows = []
    for x, y in points:
        p, rho = coordinates(x, y)
        curve_rows.append([p, rho, MASS_PER_CELL / rho, 290, x, y])
    return header, rows, header[:4] + ["pdf_x_pt", "pdf_y_pt"], curve_rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("pdf", type=Path)
    args = parser.parse_args()
    header, rows, curve_header, curve_rows = extract(args.pdf)
    for name, columns, data in [
        ("figure5-brillouin", header, rows),
        ("figure5-curve", curve_header, curve_rows),
    ]:
        path = ROOT / f"peritheos/data/datasets/argon-fcc-chen-2010-{name}.csv"
        with path.open("w", newline="") as stream:
            writer = csv.writer(stream, lineterminator="\n")
            writer.writerow(columns)
            for row in data:
                writer.writerow(
                    [f"{v:.8f}" if isinstance(v, float) else v for v in row]
                )
        print(f"{path}: {len(data)} rows")


if __name__ == "__main__":
    main()
