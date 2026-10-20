#!/usr/bin/env python3
"""Plot the structural quality and redundancy of the pooled scaffolds.

Companion to plot_inventory.py, which covers length coverage. This one asks
whether the structures are actually good bundles and how many of them are
really independent.

  1. helix-length distribution, split by topology
  2. loop-length distribution (fold conditioning cares about tight loops)
  3. bundle geometry: worst helix-pair angle vs worst pair separation,
     with the pass thresholds drawn in
  4. redundancy: how many designs came off each RFdiffusion backbone

Usage:
    python3 plot_quality.py --tsv inventory/inventory.tsv \
        --output inventory/quality.png
"""
import argparse
import csv
from collections import Counter, defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

TOPOLOGY_COLORS = {3: "#7ccfb8", 4: "#1f8a7a", 5: "#14425c"}
FAIL_COLOR = "#d94f4f"
GRID_KW = dict(linestyle="--", alpha=0.45, zorder=0)

# must match inventory.py
MAX_BUNDLE_ANGLE = 45.0
MAX_BUNDLE_DIST = 18.0


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    with open(args.tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not rows:
        raise SystemExit("no rows")

    topologies = sorted({int(r["n_helix"]) for r in rows})
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    # --- 1. helix lengths ------------------------------------------------
    ax = axes[0][0]
    for topo in topologies:
        vals = [int(x) for r in rows if int(r["n_helix"]) == topo
                for x in r["helix_lens"].split(",") if x]
        if vals:
            ax.hist(vals, bins=range(min(vals), max(vals) + 2),
                    alpha=0.75, color=TOPOLOGY_COLORS.get(topo, "#999"),
                    label=f"{topo}-helix", zorder=3)
    ax.set_title("Helix lengths")
    ax.set_xlabel("Residues per helix")
    ax.set_ylabel("Helices")
    ax.legend(frameon=False)
    ax.grid(axis="y", **GRID_KW)

    # --- 2. loop lengths -------------------------------------------------
    ax = axes[0][1]
    loops = [int(x) for r in rows for x in r["loop_lens"].split(",")
             if x and x != "-"]
    counts = Counter(loops)
    ks = sorted(counts)
    ax.bar(ks, [counts[k] for k in ks], color="#1f8a7a",
           edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_title(f"Loop lengths (median {sorted(loops)[len(loops) // 2]} res)")
    ax.set_xlabel("Residues between helices")
    ax.set_ylabel("Loops")
    ax.set_xticks(ks)
    ax.grid(axis="y", **GRID_KW)

    # --- 3. bundle geometry ----------------------------------------------
    ax = axes[1][0]
    for topo in topologies:
        sub = [r for r in rows if int(r["n_helix"]) == topo]
        ax.scatter([float(r["worst_pair_angle"]) for r in sub],
                   [float(r["worst_pair_dist"]) for r in sub],
                   s=26, alpha=0.7, edgecolor="white", linewidth=0.4,
                   color=TOPOLOGY_COLORS.get(topo, "#999"),
                   label=f"{topo}-helix", zorder=3)
    ax.axvline(MAX_BUNDLE_ANGLE, color=FAIL_COLOR, linestyle="--", linewidth=1.2)
    ax.axhline(MAX_BUNDLE_DIST, color=FAIL_COLOR, linestyle="--", linewidth=1.2)
    n_fail = sum(1 for r in rows if r["is_bundle"] == "0")
    ax.set_title(f"Bundle geometry - {n_fail} outside the dashed limits")
    ax.set_xlabel("Worst consecutive helix-pair angle (deg)")
    ax.set_ylabel("Worst pair separation (A)")
    ax.legend(frameon=False)
    ax.grid(**GRID_KW)

    # --- 4. redundancy ---------------------------------------------------
    ax = axes[1][1]
    per_backbone = defaultdict(int)
    for r in rows:
        key = (f'{r["campaign"]}__{r["backbone_id"]}' if r["backbone_id"]
               else r["source_file"])
        per_backbone[key] += 1
    spread = Counter(per_backbone.values())
    ks = sorted(spread)
    ax.bar(ks, [spread[k] for k in ks], color="#1f8a7a",
           edgecolor="white", linewidth=0.4, zorder=3)
    for k in ks:
        ax.text(k, spread[k], str(spread[k]), ha="center", va="bottom",
                fontsize=9)
    ax.set_title(f"Redundancy - {len(rows)} designs from "
                 f"{len(per_backbone)} distinct backbones")
    ax.set_xlabel("Designs sharing one backbone")
    ax.set_ylabel("Backbones")
    ax.set_xticks(ks)
    ax.grid(axis="y", **GRID_KW)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
