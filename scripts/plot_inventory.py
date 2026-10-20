#!/usr/bin/env python3
"""Plot the scaffold inventory: what lengths and topologies the pool covers.

Three panels:
  1. stacked length histogram, colored by helix count
  2. coverage of the target band -- one bar per length, so gaps are obvious
  3. distinct backbones per length (true diversity, not design count)

Usage:
    python3 plot_inventory.py --tsv inventory/inventory.tsv \
        --output inventory/coverage.png --band 60 85
"""
import argparse
import csv
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# Teal-and-green, CVD-checked; light -> dark with increasing helix count.
TOPOLOGY_COLORS = {3: "#7ccfb8", 4: "#1f8a7a", 5: "#14425c"}
GAP_COLOR = "#d94f4f"
GRID_KW = dict(axis="y", linestyle="--", alpha=0.45, zorder=0)


def load(tsv):
    with open(tsv) as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--band", nargs=2, type=int, default=[60, 85],
                    metavar=("LO", "HI"), help="target length band to audit")
    ap.add_argument("--topology", type=int, default=4,
                    help="helix count the band audit is about (default 4)")
    ap.add_argument("--bundles-only", action="store_true",
                    help="count only rows that passed the bundle geometry check")
    args = ap.parse_args()

    rows = load(args.tsv)
    if args.bundles_only:
        rows = [r for r in rows if r["is_bundle"] == "1"]
    if not rows:
        raise SystemExit("no rows to plot")

    lo, hi = args.band
    by_topology = defaultdict(lambda: defaultdict(int))
    for r in rows:
        by_topology[int(r["n_helix"])][int(r["n_res"])] += 1

    lengths = sorted({int(r["n_res"]) for r in rows})
    topologies = sorted(by_topology)

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))

    # --- panel 1: stacked histogram over every length present ------------
    ax = axes[0]
    bottom = defaultdict(int)
    for topo in topologies:
        heights = [by_topology[topo].get(L, 0) for L in lengths]
        ax.bar(lengths, heights,
               bottom=[bottom[L] for L in lengths],
               color=TOPOLOGY_COLORS.get(topo, "#999999"),
               edgecolor="white", linewidth=0.4,
               label=f"{topo}-helix", zorder=3)
        for L, h in zip(lengths, heights):
            bottom[L] += h
    ax.set_title("Scaffold pool by length and topology")
    ax.set_xlabel("Residues")
    ax.set_ylabel("Structures")
    ax.legend(frameon=False)
    ax.grid(**GRID_KW)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # --- panel 2: the target band, one slot per length -------------------
    ax = axes[1]
    band = list(range(lo, hi + 1))
    counts = [by_topology[args.topology].get(L, 0) for L in band]
    colors = [TOPOLOGY_COLORS.get(args.topology, "#1f8a7a") if c else GAP_COLOR
              for c in counts]
    ax.bar(band, [c if c else 0.35 for c in counts], color=colors,
           edgecolor="white", linewidth=0.4, zorder=3)
    n_gaps = sum(1 for c in counts if c == 0)
    ax.set_title(f"{args.topology}-helix coverage of {lo}-{hi} res "
                 f"- {n_gaps} length(s) with nothing (red)")
    ax.set_xlabel("Residues")
    ax.set_ylabel("Structures")
    ax.set_xticks(band)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(**GRID_KW)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    # --- panel 3: distinct backbones, not design count -------------------
    ax = axes[2]
    distinct = defaultdict(set)
    for r in rows:
        if int(r["n_helix"]) != args.topology:
            continue
        # designs with no backbone id are each their own structure
        key = (f'{r["campaign"]}__{r["backbone_id"]}' if r["backbone_id"]
               else r["source_file"])
        distinct[int(r["n_res"])].add(key)
    uniq = [len(distinct.get(L, ())) for L in band]
    ax.bar(band, [u if u else 0.35 for u in uniq],
           color=["#1f8a7a" if u else GAP_COLOR for u in uniq],
           edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_title(f"Distinct RFdiffusion backbones per length "
                 f"({args.topology}-helix) - the number that actually matters")
    ax.set_xlabel("Residues")
    ax.set_ylabel("Distinct backbones")
    ax.set_xticks(band)
    ax.tick_params(axis="x", labelsize=8)
    ax.grid(**GRID_KW)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
