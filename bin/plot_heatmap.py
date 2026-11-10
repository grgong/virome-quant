"""Heatmap of results/viral_rpm.tsv: samples (rows, in samplesheet order) x high-quality vOTUs
(columns), reads per million on a log colour scale; grey = not detected.

Usage: python3 -I plot_heatmap.py <run_dir>   -> results/viral_rpm.{pdf,png}
"""
import csv
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm

RUN = sys.argv[1]
rows = list(csv.reader(open(f"{RUN}/results/viral_rpm.tsv"), delimiter="\t"))
head, body = rows[0], rows[1:]
first_votu = next((i for i, h in enumerate(head) if h.endswith("]")), len(head))
votus = head[first_votu:]
if not votus or not body:
    print("no high-quality vOTU detected: no heatmap")
    sys.exit(0)
meta = [i for i in range(1, first_votu) if head[i] != "flowcell_lane"]
labels = [r[0] + "  " + " | ".join(r[i] for i in meta) for r in body]   # sample, then its metadata
m = np.array([[float(x) for x in r[first_votu:]] for r in body])

# the axes are sized by the cells; bbox_inches="tight" below adds the room the labels take
fig, ax = plt.subplots(figsize=(1.0 + 0.5 * len(votus), 0.6 + 0.25 * len(body)))
ax.set_facecolor("0.9")
pos = m[m > 0]
im = ax.imshow(np.ma.masked_less_equal(m, 0), aspect="auto", cmap="viridis",
               norm=LogNorm(vmin=max(pos.min(), 0.1), vmax=pos.max()), interpolation="nearest")
ax.set_xticks(range(len(votus)), votus, rotation=60, ha="right", fontsize=7)
ax.set_yticks(range(len(body)), labels, fontsize=7)
ax.set_xticks(np.arange(-0.5, len(votus)), minor=True)
ax.set_yticks(np.arange(-0.5, len(body)), minor=True)
ax.grid(which="minor", color="white", linewidth=0.8)
ax.tick_params(which="minor", length=0)
fig.colorbar(im, ax=ax, shrink=0.6).set_label("reads per million (grey: not detected)", fontsize=7)
for ext in ("pdf", "png"):
    fig.savefig(f"{RUN}/results/viral_rpm.{ext}", dpi=200, bbox_inches="tight")
print(f"{RUN}/results/viral_rpm.pdf: {len(body)} samples x {len(votus)} vOTUs")
