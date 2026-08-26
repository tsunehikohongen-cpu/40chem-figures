"""
Figure 2A -- PCA of the 40-chemical high-concentration samples (600 dpi).

Inputs (see README.md for column formats):
  --input-csv   DESeq-normalized expression matrix
                (default: data/40_Chem_DESeq_normalized.csv)
  --label-xlsx  Label A workbook, sheet "Labels A_B"
                (default: data/Table_Labels_A.xlsx)
Output:
  --out-dir     directory for the PNG/TIFF figure and score/loading CSVs
                (default: output/)
"""
import argparse
import csv
import math
import os
import re
from collections import defaultdict

import numpy as np
import pandas as pd
from PIL import Image, ImageDraw, ImageFont

DEFAULT_INPUT_CSV = "data/40_Chem_DESeq_normalized.csv"
DEFAULT_LABEL_XLSX = "data/Table_Labels_A.xlsx"
DEFAULT_OUT_DIR = "output"
FALLBACK_FONTS = [
    "/System/Library/Fonts/Helvetica.ttc",   # macOS
    "/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",  # Linux
    "C:\\Windows\\Fonts\\arial.ttf",         # Windows
]


def resolve_font():
    for path in FALLBACK_FONTS:
        if os.path.exists(path):
            return path
    raise FileNotFoundError(
        "No usable font found; set FONT_PATH env var to a .ttf/.ttc file."
    )


def load_font(font_path, size, bold=False):
    """Bold face is font-collection index 1 (only present in .ttc files like
    macOS's Helvetica.ttc); single-face .ttf fallbacks only have index 0."""
    if bold:
        try:
            return ImageFont.truetype(font_path, size, index=1)
        except OSError:
            pass
    return ImageFont.truetype(font_path, size, index=0)


def parse_sample_column(name):
    match = re.match(r"^(.*?)[ _]([0-9]+(?:\.[0-9]+)?)$", name)
    if not match:
        raise ValueError(f"Could not parse sample column: {name}")
    return match.group(1), float(match.group(2))


def canonical_chemical(name):
    if not isinstance(name, str):
        return ""
    cleaned = re.sub(r"\s*\([^)]*\)", "", name).strip()
    replacements = {
        "5-Acetylsalicylic acid": "Acetylsalicylic acid",
        "deltamethrin": "Deltamethrin",
    }
    cleaned = replacements.get(cleaned, cleaned)
    return re.sub(r"\s+", " ", cleaned).lower()


def load_label_a(label_xlsx):
    labels = pd.read_excel(label_xlsx, sheet_name="Labels A_B", header=2)
    labels = labels.dropna(subset=["Chemical", "Label A"])
    category_by_chemical = {}
    for _, row in labels.iterrows():
        key = canonical_chemical(row["Chemical"])
        if key:
            value = str(row["Label A"]).strip()
            if value in {"1", "1.0"}:
                value = "Category1"
            elif value in {"0", "0.0"}:
                value = "Category0"
            elif value.lower() == "unknown":
                value = "Unknown"
            category_by_chemical[key] = value
    return category_by_chemical


def convex_hull(points):
    pts = sorted(set((float(x), float(y)) for x, y in points))
    if len(pts) <= 1:
        return pts

    def cross(o, a, b):
        return (a[0] - o[0]) * (b[1] - o[1]) - (a[1] - o[1]) * (b[0] - o[0])

    lower = []
    for p in pts:
        while len(lower) >= 2 and cross(lower[-2], lower[-1], p) <= 0:
            lower.pop()
        lower.append(p)
    upper = []
    for p in reversed(pts):
        while len(upper) >= 2 and cross(upper[-2], upper[-1], p) <= 0:
            upper.pop()
        upper.append(p)
    return lower[:-1] + upper[:-1]


def expand_polygon(points, factor=1.12):
    if not points:
        return points
    cx = sum(x for x, _ in points) / len(points)
    cy = sum(y for _, y in points) / len(points)
    return [(cx + (x - cx) * factor, cy + (y - cy) * factor) for x, y in points]


def covariance_ellipse(points, n_std=1.45, steps=160):
    arr = np.asarray(points, dtype=float)
    center = arr.mean(axis=0)
    if len(arr) < 3:
        return []
    cov = np.cov(arr.T)
    vals, vecs = np.linalg.eigh(cov)
    order = vals.argsort()[::-1]
    vals, vecs = vals[order], vecs[:, order]
    vals = np.maximum(vals, 1e-9)
    angles = np.linspace(0, 2 * np.pi, steps, endpoint=False)
    circle = np.stack([np.cos(angles), np.sin(angles)], axis=1)
    ellipse = circle @ np.diag(np.sqrt(vals) * n_std) @ vecs.T + center
    return [tuple(p) for p in ellipse]


def draw_dashed_line(draw, points, fill, width=4, dash=22, gap=16):
    for start, end in zip(points, points[1:]):
        x1, y1 = start
        x2, y2 = end
        length = math.hypot(x2 - x1, y2 - y1)
        if length == 0:
            continue
        dx, dy = (x2 - x1) / length, (y2 - y1) / length
        dist = 0
        while dist < length:
            seg_end = min(dist + dash, length)
            draw.line([x1 + dx * dist, y1 + dy * dist, x1 + dx * seg_end, y1 + dy * seg_end], fill=fill, width=width)
            dist += dash + gap


def draw_translucent_polygon(base, points, fill, alpha):
    overlay = Image.new("RGBA", base.size, (255, 255, 255, 0))
    odraw = ImageDraw.Draw(overlay)
    odraw.polygon(points, fill=fill + (alpha,))
    base.alpha_composite(overlay)


def marker_polygon(cx, cy, radius, sides, rotation=-math.pi / 2):
    return [
        (
            cx + radius * math.cos(rotation + 2 * math.pi * i / sides),
            cy + radius * math.sin(rotation + 2 * math.pi * i / sides),
        )
        for i in range(sides)
    ]


def star_polygon(cx, cy, radius, rotation=-math.pi / 2):
    pts = []
    for i in range(10):
        r = radius if i % 2 == 0 else radius * 0.45
        a = rotation + math.pi * i / 5
        pts.append((cx + r * math.cos(a), cy + r * math.sin(a)))
    return pts


def draw_marker(draw, x, y, shape, color, radius=19, fill=True, width=5):
    fill_color = color if fill else (255, 255, 255)
    if shape == "circle":
        draw.ellipse([x - radius, y - radius, x + radius, y + radius], fill=fill_color, outline=color, width=width)
    elif shape == "square":
        draw.rectangle([x - radius, y - radius, x + radius, y + radius], fill=fill_color, outline=color, width=width)
    elif shape == "diamond":
        draw.polygon([(x, y - radius), (x + radius, y), (x, y + radius), (x - radius, y)], fill=fill_color, outline=color)
        if not fill:
            draw.line([(x, y - radius), (x + radius, y), (x, y + radius), (x - radius, y), (x, y - radius)], fill=color, width=width)
    elif shape == "triangle":
        pts = marker_polygon(x, y, radius * 1.15, 3)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "tri_down":
        pts = marker_polygon(x, y, radius * 1.15, 3, rotation=math.pi / 2)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "tri_left":
        pts = marker_polygon(x, y, radius * 1.15, 3, rotation=math.pi)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "tri_right":
        pts = marker_polygon(x, y, radius * 1.15, 3, rotation=0)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "pentagon":
        pts = marker_polygon(x, y, radius * 1.05, 5)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "hexagon":
        pts = marker_polygon(x, y, radius * 1.05, 6)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "octagon":
        pts = marker_polygon(x, y, radius * 1.05, 8)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "star":
        pts = star_polygon(x, y, radius * 1.15)
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "rect_h":
        draw.rectangle([x - radius * 1.35, y - radius * 0.72, x + radius * 1.35, y + radius * 0.72], fill=fill_color, outline=color, width=width)
    elif shape == "rect_v":
        draw.rectangle([x - radius * 0.72, y - radius * 1.35, x + radius * 0.72, y + radius * 1.35], fill=fill_color, outline=color, width=width)
    elif shape == "bowtie":
        pts = [(x - radius, y - radius), (x, y), (x - radius, y + radius), (x + radius, y + radius), (x, y), (x + radius, y - radius)]
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "hourglass":
        pts = [(x - radius, y - radius), (x + radius, y - radius), (x, y), (x + radius, y + radius), (x - radius, y + radius), (x, y)]
        draw.polygon(pts, fill=fill_color, outline=color)
        draw.line(pts + [pts[0]], fill=color, width=width)
    elif shape == "asterisk":
        for angle in (0, math.pi / 3, 2 * math.pi / 3):
            dx = radius * math.cos(angle)
            dy = radius * math.sin(angle)
            draw.line([x - dx, y - dy, x + dx, y + dy], fill=color, width=width + 1)
    elif shape == "line_h":
        draw.line([x - radius * 1.25, y, x + radius * 1.25, y], fill=color, width=width + 3)
    elif shape == "line_v":
        draw.line([x, y - radius * 1.25, x, y + radius * 1.25], fill=color, width=width + 3)
    elif shape == "plus":
        draw.line([x - radius, y, x + radius, y], fill=color, width=width + 2)
        draw.line([x, y - radius, x, y + radius], fill=color, width=width + 2)
    elif shape == "x":
        draw.line([x - radius, y - radius, x + radius, y + radius], fill=color, width=width + 2)
        draw.line([x - radius, y + radius, x + radius, y - radius], fill=color, width=width + 2)


def nice_ticks(vmin, vmax, n=5):
    if math.isclose(vmin, vmax):
        return [vmin]
    span = vmax - vmin
    step0 = span / max(n - 1, 1)
    mag = 10 ** math.floor(math.log10(abs(step0)))
    residual = step0 / mag
    if residual <= 1:
        step = mag
    elif residual <= 2:
        step = 2 * mag
    elif residual <= 5:
        step = 5 * mag
    else:
        step = 10 * mag
    start = math.ceil(vmin / step) * step
    ticks = []
    val = start
    while val <= vmax + step * 0.5:
        ticks.append(0.0 if abs(val) < step * 1e-6 else val)
        val += step
    return ticks


def fmt_tick(value):
    if abs(value) >= 10:
        return f"{value:.0f}"
    if abs(value) >= 1:
        return f"{value:.1f}".rstrip("0").rstrip(".")
    return f"{value:.2f}".rstrip("0").rstrip(".")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input-csv", default=os.environ.get("INPUT_CSV", DEFAULT_INPUT_CSV))
    p.add_argument("--label-xlsx", default=os.environ.get("LABEL_XLSX", DEFAULT_LABEL_XLSX))
    p.add_argument("--out-dir", default=os.environ.get("OUT_DIR", DEFAULT_OUT_DIR))
    p.add_argument("--font-path", default=os.environ.get("FONT_PATH", ""))
    return p.parse_args()


def main():
    args = parse_args()
    font_path = args.font_path or resolve_font()
    os.makedirs(args.out_dir, exist_ok=True)
    df = pd.read_csv(args.input_csv, encoding="utf-8-sig")
    sample_cols = list(df.columns[1:])

    grouped = defaultdict(list)
    for col in sample_cols:
        chemical, dose = parse_sample_column(col)
        grouped[chemical].append((dose, col))

    selected = []
    for chemical in sorted(grouped):
        if chemical.startswith("DMSO"):
            continue
        dose, col = max(grouped[chemical], key=lambda item: item[0])
        selected.append((chemical, dose, col))

    expression = df[[col for _, _, col in selected]].apply(pd.to_numeric, errors="coerce").fillna(0.0)
    normalized_expression = expression.to_numpy(dtype=float)
    variance = normalized_expression.var(axis=1)
    keep = variance > 0
    genes = df.iloc[:, 0].astype(str).to_numpy()[keep]
    x = normalized_expression[keep, :].T

    x_centered = x - x.mean(axis=0, keepdims=True)
    u, s, vt = np.linalg.svd(x_centered, full_matrices=False)
    scores = u[:, :2] * s[:2]
    explained = (s**2) / (len(selected) - 1)
    explained_ratio = explained / explained.sum()
    loadings = vt[:2, :].T

    category_by_chemical = load_label_a(args.label_xlsx)
    categories = []
    missing_labels = []
    for chemical, _, _ in selected:
        category = category_by_chemical.get(canonical_chemical(chemical))
        if category is None:
            category = "Unknown"
            missing_labels.append(chemical)
        categories.append(category)

    score_csv = os.path.join(args.out_dir, "pca_high_concentration_scores.csv")
    with open(score_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["chemical", "dose", "sample_column", "Label A category", "PC1", "PC2"])
        for (chemical, dose, col), category, (pc1, pc2) in zip(selected, categories, scores):
            writer.writerow([chemical, dose, col, category, pc1, pc2])

    loading_csv = os.path.join(args.out_dir, "pca_high_concentration_gene_loadings.csv")
    with open(loading_csv, "w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f)
        writer.writerow(["gene", "PC1_loading", "PC2_loading"])
        for gene, (pc1, pc2) in zip(genes, loadings):
            writer.writerow([gene, pc1, pc2])

    dpi = 600
    width, height = 7200, 6000
    img = Image.new("RGBA", (width, height), "white")
    draw = ImageDraw.Draw(img)
    font = load_font(font_path, 110)
    bold_font = load_font(font_path, 116, bold=True)
    small = load_font(font_path, 82)
    tiny = load_font(font_path, 72)
    title_font = load_font(font_path, 132)
    group_font = load_font(font_path, 136, bold=True)

    left, top, right, bottom = 860, 980, 6480, 5270
    gray = (165, 170, 174)
    dark = (66, 66, 66)
    axis = (138, 138, 138)
    red_axis = (238, 150, 150)
    grid = (220, 224, 228)
    group_style = {
        "Category1": {"color": (186, 84, 174), "label": "Category1"},
        "Category0": {"color": (239, 131, 58), "label": "Category0"},
        "Unknown": {"color": (28, 176, 185), "label": "Unknown"},
    }

    area_points = [tuple(row) for row in scores]
    for category in ["Category1", "Category0", "Unknown"]:
        points = [tuple(score) for score, cat in zip(scores, categories) if cat == category]
        if len(points) >= 3:
            area_points.extend(expand_polygon(convex_hull(points), 1.08))
            area_points.extend(covariance_ellipse(points))
    area_points = np.asarray(area_points, dtype=float)
    sx_min, sx_max = area_points[:, 0].min(), area_points[:, 0].max()
    sy_min, sy_max = area_points[:, 1].min(), area_points[:, 1].max()
    x_pad = max((sx_max - sx_min) * 0.12, 8)
    y_pad = max((sy_max - sy_min) * 0.16, 12)
    sx_min, sx_max = sx_min - x_pad, sx_max + x_pad
    sy_min, sy_max = sy_min - y_pad, sy_max + y_pad

    lx_min, lx_max = np.percentile(loadings[:, 0], [0.5, 99.5])
    ly_min, ly_max = np.percentile(loadings[:, 1], [0.5, 99.5])
    scale = 0.72 * min((sx_max - sx_min) / (lx_max - lx_min), (sy_max - sy_min) / (ly_max - ly_min))
    loading_plot = loadings * scale

    def px(xval):
        return left + (xval - sx_min) / (sx_max - sx_min) * (right - left)

    def py(yval):
        return bottom - (yval - sy_min) / (sy_max - sy_min) * (bottom - top)

    draw.text((left, 170), "PCA: high concentration samples", fill=dark, font=title_font)
    draw.text((left, 335), "Label A groups; DESeq-normalized values; gene-wise centered PCA; DMSO controls excluded", fill=(105, 105, 105), font=small)

    draw.rectangle([left, top, right, bottom], outline=axis, width=5)
    for tick in nice_ticks(sx_min, sx_max, 7):
        if tick < sx_min or tick > sx_max:
            continue
        xpix = px(tick)
        draw.line([xpix, bottom, xpix, bottom + 52], fill=axis, width=4)
        label = fmt_tick(tick)
        bbox = draw.textbbox((0, 0), label, font=small)
        draw.text((xpix - (bbox[2] - bbox[0]) / 2, bottom + 76), label, fill=dark, font=small)
    for tick in nice_ticks(sy_min, sy_max, 7):
        if tick < sy_min or tick > sy_max:
            continue
        ypix = py(tick)
        draw.line([left - 52, ypix, left, ypix], fill=axis, width=4)
        label = fmt_tick(tick)
        bbox = draw.textbbox((0, 0), label, font=small)
        draw.text((left - 92 - (bbox[2] - bbox[0]), ypix - 36), label, fill=dark, font=small)

    if sx_min < 0 < sx_max:
        draw.line([px(0), top, px(0), bottom], fill=grid, width=4)
    if sy_min < 0 < sy_max:
        draw.line([left, py(0), right, py(0)], fill=grid, width=4)

    for tick in nice_ticks(lx_min, lx_max, 6):
        if tick < lx_min or tick > lx_max:
            continue
        xpix = px(tick * scale)
        draw.line([xpix, top - 52, xpix, top], fill=red_axis, width=4)
        label = fmt_tick(tick)
        bbox = draw.textbbox((0, 0), label, font=small)
        draw.text((xpix - (bbox[2] - bbox[0]) / 2, top - 145), label, fill=dark, font=small)
    draw.line([left, top, right, top], fill=red_axis, width=3)
    for tick in nice_ticks(ly_min, ly_max, 6):
        if tick < ly_min or tick > ly_max:
            continue
        ypix = py(tick * scale)
        draw.line([right, ypix, right + 52, ypix], fill=red_axis, width=4)
        label = fmt_tick(tick)
        draw.text((right + 78, ypix - 36), label, fill=dark, font=small)
    draw.line([right, top, right, bottom], fill=red_axis, width=3)

    rng = np.random.default_rng(4)
    if len(loading_plot) > 12000:
        idx = rng.choice(len(loading_plot), 12000, replace=False)
        loading_to_draw = loading_plot[idx]
    else:
        loading_to_draw = loading_plot
    for xval, yval in loading_to_draw:
        if sx_min <= xval <= sx_max and sy_min <= yval <= sy_max:
            xp, yp = px(xval), py(yval)
            draw.ellipse([xp - 3, yp - 3, xp + 3, yp + 3], fill=gray + (150,))

    for category in ["Category1", "Category0", "Unknown"]:
        points = [tuple(score) for score, cat in zip(scores, categories) if cat == category]
        if len(points) < 3:
            continue
        color = group_style[category]["color"]
        hull = expand_polygon(convex_hull(points), 1.08)
        hull_px = [(px(xval), py(yval)) for xval, yval in hull]
        ellipse = covariance_ellipse(points)
        ellipse_px = [(px(xval), py(yval)) for xval, yval in ellipse]
        draw_translucent_polygon(img, ellipse_px, color, 55)
        draw_translucent_polygon(img, hull_px, color, 22)
        draw = ImageDraw.Draw(img)
        draw_dashed_line(draw, hull_px + [hull_px[0]], fill=color + (210,), width=5, dash=28, gap=18)
        cx, cy = np.asarray(points).mean(axis=0)
        label = group_style[category]["label"]
        bbox = draw.textbbox((0, 0), label, font=group_font, stroke_width=10)
        tx = px(cx) - (bbox[2] - bbox[0]) / 2
        ty = py(cy) - (bbox[3] - bbox[1]) / 2
        draw.text(
            (tx, ty),
            label,
            fill=color + (255,),
            font=group_font,
            stroke_width=10,
            stroke_fill=(255, 255, 255, 235),
        )

    shapes = {"Category1": "circle", "Category0": "square", "Unknown": "triangle"}

    for i, ((chemical, _, _), category, (xval, yval)) in enumerate(zip(selected, categories, scores)):
        color = group_style[category]["color"]
        shape = shapes[category]
        xp, yp = px(xval), py(yval)
        draw_marker(draw, xp, yp, shape, color, radius=34, fill=True, width=7)

    xlabel = f"PC1 ({explained_ratio[0] * 100:.1f}%)"
    ylabel = f"PC2 ({explained_ratio[1] * 100:.1f}%)"
    bbox = draw.textbbox((0, 0), xlabel, font=bold_font)
    draw.text(((left + right) / 2 - (bbox[2] - bbox[0]) / 2, height - 235), xlabel, fill=dark, font=bold_font)
    rotated = Image.new("RGBA", (850, 230), (255, 255, 255, 0))
    rdraw = ImageDraw.Draw(rotated)
    rdraw.text((0, 30), ylabel, fill=dark, font=bold_font)
    rotated = rotated.rotate(90, expand=True)
    img.paste(rotated, (135, int((top + bottom) / 2 - rotated.height / 2)), rotated)

    legend_x = left + 90
    legend_y = 655
    for i, category in enumerate(["Category1", "Category0", "Unknown"]):
        color = group_style[category]["color"]
        x0 = legend_x + i * 930
        draw_marker(draw, x0, legend_y, shapes[category], color, radius=42, fill=True, width=8)
        count = sum(1 for cat in categories if cat == category)
        draw.text((x0 + 78, legend_y - 48), f"{category} (n={count})", fill=dark, font=small)

    png = os.path.join(args.out_dir, "pca_high_concentration_600dpi.png")
    tif = os.path.join(args.out_dir, "pca_high_concentration_600dpi.tiff")
    rgb = img.convert("RGB")
    rgb.save(png, dpi=(dpi, dpi))
    rgb.save(tif, dpi=(dpi, dpi), compression="tiff_lzw")
    print(f"Selected high-concentration samples: {len(selected)}")
    print(f"Genes retained after zero-variance filtering: {keep.sum()}")
    print(f"PC1: {explained_ratio[0] * 100:.2f}%")
    print(f"PC2: {explained_ratio[1] * 100:.2f}%")
    if missing_labels:
        print("Missing labels:", ", ".join(missing_labels))
    print(png)
    print(tif)
    print(score_csv)
    print(loading_csv)


if __name__ == "__main__":
    main()
