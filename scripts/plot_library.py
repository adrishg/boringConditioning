#!/usr/bin/env python3
"""Plot the clustered library: representatives per topology and length.

Two panels:
  1. a topology x length grid of representative counts -- the library as it
     will actually be handed to RFdiffusion, with empty cells called out
  2. how far each length collapsed under clustering, which is where the
     redundancy hiding in the raw counts becomes visible

Usage:
    python3 plot_library.py --clusters inventory/clusters.tsv \
        --output inventory/library.png --band 58 85
"""
import argparse
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

TEAL = mcolors.LinearSegmentedColormap.from_list(
    "teal", ["#eef7f4", "#9fd8c4", "#3faf96", "#1f8a7a", "#14425c"])
EMPTY = "#f6dcd6"
GRID_KW = dict(linestyle="--", alpha=0.4, zorder=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--band", nargs=2, type=int, default=[58, 85],
                    metavar=("LO", "HI"))
    ap.add_argument("--topology", type=int, default=4,
                    help="topology the gap audit is about (default 4)")
    args = ap.parse_args()
    lo, hi = args.band

    with open(args.clusters) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    reps = [r for r in rows if r["is_representative"] == "1"]

    grid = defaultdict(int)
    total = defaultdict(int)
    for r in rows:
        key = (int(r["n_helix"]), int(r["n_res"]))
        total[key] += 1
        if r["is_representative"] == "1":
            grid[key] += 1

    topos = sorted({int(r["n_helix"]) for r in rows})
    lengths = list(range(min(k[1] for k in total), max(k[1] for k in total) + 1))

    fig, axes = plt.subplots(2, 1, figsize=(14, 8),
                             gridspec_kw={"height_ratios": [1.15, 1]})

    # --- panel 1: topology x length grid --------------------------------
    ax = axes[0]
    mat = np.array([[grid.get((t, L), 0) for L in lengths] for t in topos],
                   dtype=float)
    vmax = max(mat.max(), 1)
    ax.imshow(np.ma.masked_where(mat == 0, mat), cmap=TEAL, aspect="auto",
              vmin=0, vmax=vmax,
              extent=(lengths[0] - .5, lengths[-1] + .5, len(topos) - .5, -.5))
    # empty cells
    for i, t in enumerate(topos):
        for j, L in enumerate(lengths):
            if mat[i, j] == 0:
                ax.add_patch(plt.Rectangle((L - .5, i - .5), 1, 1,
                                           color=EMPTY, zorder=2))
            else:
                ax.text(L, i, int(mat[i, j]), ha="center", va="center",
                        fontsize=7, zorder=3,
                        color="white" if mat[i, j] > vmax * .55 else "#14425c")
    ax.axvline(lo - .5, color="#b2432f", lw=1.2, ls="--", zorder=4)
    ax.axvline(hi + .5, color="#b2432f", lw=1.2, ls="--", zorder=4)
    ax.set_yticks(range(len(topos)))
    ax.set_yticklabels([f"{t}-helix" for t in topos])
    ax.set_xticks([L for L in lengths if L % 2 == 0])
    ax.tick_params(axis="x", labelsize=7, rotation=90)
    ax.set_xlabel("Residues")
    gaps = [L for L in range(lo, hi + 1) if not grid.get((args.topology, L))]
    ax.set_title(
        f"Library after clustering: {len(reps)} representatives "
        f"(from {len(rows)} structures)\n"
        f"{args.topology}-helix gaps in {lo}-{hi}: "
        f"{', '.join(map(str, gaps)) if gaps else 'none'}   "
        f"(dashed lines mark the band; pink = nothing)")

    # --- panel 2: collapse ratio ----------------------------------------
    ax = axes[1]
    t = args.topology
    designs = [total.get((t, L), 0) for L in lengths]
    folds = [grid.get((t, L), 0) for L in lengths]
    ax.bar(lengths, designs, color="#cfe6de", edgecolor="white", linewidth=.4,
           label=f"{t}-helix designs", zorder=3)
    ax.bar(lengths, folds, color="#1f8a7a", edgecolor="white", linewidth=.4,
           label="distinct folds kept", zorder=4)
    for L, d, f in zip(lengths, designs, folds):
        if d >= 20:
            ax.text(L, d, f"{d}→{f}", ha="center", va="bottom", fontsize=7,
                    color="#b2432f", rotation=90)
    ax.axvspan(lo - .5, hi + .5, color="#1f8a7a", alpha=.06, zorder=1)
    ax.set_yscale("symlog")
    ax.set_xlabel("Residues")
    ax.set_ylabel("Count (symlog)")
    ax.set_xticks([L for L in lengths if L % 2 == 0])
    ax.tick_params(axis="x", labelsize=7, rotation=90)
    ax.legend(frameon=False, fontsize=9)
    ax.grid(axis="y", **GRID_KW)
    ax.set_title(f"How much of the {t}-helix pool was redundant "
                 f"(labels on groups of 20+)")

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
