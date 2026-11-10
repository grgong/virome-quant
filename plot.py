"""Heatmap of viral_rpm.tsv: vOTUs in rows, samples in columns, reads per million in five bins.

    python3 plot.py results/viral_rpm.tsv [-s samples.txt] [-g column] [-o prefix]

-s  text file with the samples to draw, one name per line, drawn in that order (default: all)
-g  column of the table (for example condition) that colours and groups the samples; samples with
    the same value next to each other form a group
-o  output prefix, writes <prefix>.pdf and <prefix>.png (default: next to the table)

vOTUs without reads in the drawn samples are left out. Rows are ordered by their highest RPM.
"""
import argparse
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import BoundaryNorm, ListedColormap
from matplotlib.patches import Patch, Rectangle

BIN_EDGES = [1, 10, 100, 1e3, 1e4]
BIN_COLS = ["#E6E6E6", "#CAC5F3", "#A79FE1", "#8175CB", "#594CA0", "#312271"]
BIN_LABELS = ["< 1", "1–10", "10–100", "100–1,000", "1,000–10,000", "≥ 10,000"]
GROUP_COLS = ["#0072B2", "#D55E00", "#009E73", "#E69F00", "#56B4E9", "#CC79A7", "#F0E442"]  # Okabe-Ito
FONT = 6                                  # pt
MM = 1 / 25.4                             # inch
CELL, GAP, STRIP = 3.2 * MM, 1.5 * MM, 1.2 * MM

plt.rcParams.update({"font.family": "sans-serif", "font.size": FONT, "pdf.fonttype": 42,
                     "font.sans-serif": ["Arial", "Liberation Sans", "Helvetica", "DejaVu Sans"]})

ap = argparse.ArgumentParser()
ap.add_argument("table")
ap.add_argument("-s", "--samples")
ap.add_argument("-g", "--group")
ap.add_argument("-o", "--out")
args = ap.parse_args()
out = args.out or args.table.rsplit(".", 1)[0] + "_heatmap"

rows = list(csv.reader(open(args.table), delimiter="\t"))
head, body = rows[0], rows[1:]
first = next((i for i, h in enumerate(head) if h.endswith("]")), None)
if first is None:
    sys.exit("no vOTU column in " + args.table)
by_name = {r[0]: r for r in body}

if args.samples:
    names = [l.strip() for l in open(args.samples) if l.strip() and not l.startswith("#")]
    names = list(dict.fromkeys(names))
    missing = [n for n in names if n not in by_name]
    if missing:
        sys.exit("not in the table: " + ", ".join(missing))
else:
    names = [r[0] for r in body]
if not names:
    sys.exit("no samples to draw")

rpm = np.array([[float(x) for x in by_name[n][first:]] for n in names]).T      # vOTU x sample
labels = np.array([h.rsplit(" [", 1)[0] for h in head[first:]])
if len(set(labels)) < len(labels):
    labels = np.array(head[first:])
keep = np.flatnonzero(rpm.max(axis=1) > 0)
if not len(keep):
    sys.exit("no vOTU has reads in these samples")
keep = keep[np.argsort(-rpm[keep].max(axis=1))]
rpm, labels = rpm[keep], labels[keep]
bins = np.digitize(rpm, BIN_EDGES)
n_row, n_col = bins.shape

# runs of samples with the same group value; each run is drawn on its own, with a gap before it
if args.group:
    if args.group not in head[:first]:
        sys.exit(f"no column {args.group} in the table; columns: {', '.join(head[1:first])}")
    value = [by_name[n][head.index(args.group)] for n in names]
else:
    value = [""] * n_col
runs, start = [], 0
for i in range(1, n_col + 1):
    if i == n_col or value[i] != value[start]:
        runs.append((start, i, value[start]))
        start = i
# a group keeps its colour whichever samples are drawn
if args.group:
    every = dict.fromkeys(r[head.index(args.group)] for r in body)
    colour = {v: GROUP_COLS[i % len(GROUP_COLS)] for i, v in enumerate(every)}
gap = GAP / CELL                          # in cell widths
width = n_col + gap * (len(runs) - 1)

# the axes is exactly as large as its cells, so a data unit is one square cell
left, right, top, bottom = 4.0, 3.0, 2.5, 0.5
W, H = left + width * CELL + right, bottom + n_row * CELL + top
fig = plt.figure(figsize=(W, H))
ax = fig.add_axes([left / W, bottom / H, width * CELL / W, n_row * CELL / H])
cmap, norm = ListedColormap(BIN_COLS), BoundaryNorm(np.arange(-0.5, 6.5), 6)
centres, starts, x = [], [], 0.0
for lo, hi, v in runs:
    starts.append(x)
    ax.pcolormesh(x + np.arange(hi - lo + 1), np.arange(n_row + 1), bins[:, lo:hi], cmap=cmap, norm=norm,
                  edgecolors="white", linewidth=0.8)
    centres += list(x + np.arange(hi - lo) + 0.5)
    if args.group:
        ax.add_patch(Rectangle((x, -0.15 - STRIP / CELL), hi - lo, STRIP / CELL, color=colour[v], lw=0, clip_on=False))
    x += hi - lo + gap
ax.set_xlim(0, width)
ax.set_ylim(n_row, 0)
for side in ax.spines.values():
    side.set_visible(False)
ax.set_yticks(np.arange(n_row) + 0.5, labels)
ax.xaxis.tick_top()
pad = 2 + (0.15 + STRIP / CELL) * CELL * 72 if args.group else 2          # clear of the strip, in pt
ax.set_xticks(centres, names, rotation=90)
ax.tick_params(length=0)
ax.tick_params(axis="x", pad=pad)

if args.group:
    fig.canvas.draw()
    label_h = max(t.get_window_extent().height for t in ax.get_xticklabels()) * 72 / fig.dpi
    for (lo, hi, v), x in zip(runs, starts):
        ax.annotate(v, (x + (hi - lo) / 2, 0), xytext=(0, pad + label_h + 3), textcoords="offset points",
                    ha="center", va="bottom")

handles = [Patch(facecolor=c, edgecolor="white", linewidth=0.8, label=l) for c, l in zip(BIN_COLS, BIN_LABELS)]
box = CELL * 72 / FONT                    # legend square, in font sizes
ax.legend(handles=handles, title="Reads per million", loc="upper left", bbox_to_anchor=(1.0, 1.0),
          bbox_transform=ax.transAxes, borderaxespad=1.2, frameon=False, alignment="left", borderpad=0,
          handlelength=box, handleheight=box, labelspacing=0, handletextpad=0.5,
          title_fontproperties={"weight": "bold", "size": FONT})

for ext in ("pdf", "png"):
    fig.savefig(f"{out}.{ext}", dpi=300, bbox_inches="tight", pad_inches=0.05)
size = fig.get_tightbbox().width * 25.4
print(f"{n_row} vOTUs x {n_col} samples -> {out}.pdf/.png ({size:.0f} mm wide)")
if size > 183:
    print("wider than 183 mm (Nature's limit): draw fewer samples")
