#!/usr/bin/env python3
"""Plot the staged complexes waiting for binder extraction.

Reads every complexes/*/manifest.csv written by harvest_complexes.py and
shows what the harvest actually bought: which binder lengths arrived, how
good the multimer predictions are, and how much of it is redundant.

Usage:
    python3 plot_complexes.py --complexes-dir complexes \
        --output inventory/complexes.png
"""
import argparse
import csv
import glob
import os
from collections import defaultdict

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker

# teal-and-green ramp, sampled to however many batches there are so that
# no two batches ever land on the same hue
import matplotlib.colors as mcolors
TEAL_RAMP = mcolors.LinearSegmentedColormap.from_list(
    "teal", ["#bfe8d9", "#7ccfb8", "#3faf96", "#1f8a7a", "#2e6f72", "#14425c"])
GRID_KW = dict(linestyle="--", alpha=0.45, zorder=0)


def as_float(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--complexes-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--band", nargs=2, type=int, default=[60, 85],
                    metavar=("LO", "HI"),
                    help="target band, drawn for reference")
    args = ap.parse_args()

    batches = {}
    for path in sorted(glob.glob(os.path.join(args.complexes_dir,
                                              "*", "manifest.csv"))):
        name = os.path.basename(os.path.dirname(path))
        with open(path) as fh:
            batches[name] = list(csv.DictReader(fh))
    if not batches:
        raise SystemExit(f"no manifest.csv under {args.complexes_dir}")

    names = list(batches)
    denom = max(len(names) - 1, 1)
    colors = {n: TEAL_RAMP(i / denom) for i, n in enumerate(names)}
    total = sum(len(v) for v in batches.values())

    fig, axes = plt.subplots(3, 1, figsize=(11, 11))

    # --- 1. binder lengths, stacked by batch -----------------------------
    ax = axes[0]
    per_batch = {n: defaultdict(int) for n in names}
    for n, rows in batches.items():
        for r in rows:
            per_batch[n][int(r["binder_len"])] += 1
    lengths = sorted({L for d in per_batch.values() for L in d})
    per_len_total = defaultdict(int)
    for d in per_batch.values():
        for L, v in d.items():
            per_len_total[L] += v
    bottom = defaultdict(int)
    for n in names:
        h = [per_batch[n].get(L, 0) for L in lengths]
        ax.bar(lengths, h, bottom=[bottom[L] for L in lengths],
               color=colors[n], edgecolor="white", linewidth=0.4,
               label=f"{n} (n={len(batches[n])})", zorder=3)
        for L, v in zip(lengths, h):
            bottom[L] += v
    lo, hi = args.band
    ax.axvspan(lo - 0.5, hi + 0.5, color="#1f8a7a", alpha=0.08, zorder=1)
    in_band = sorted(L for L in lengths if lo <= L <= hi)
    empty = [L for L in range(lo, hi + 1) if L not in per_len_total]
    note = (f"target band {lo}-{hi}: "
            f"{sum(per_len_total[L] for L in in_band)} structures\n"
            f"still empty: {', '.join(map(str, empty)) if empty else 'none'}")
    ax.text(0.015, 0.95, note, transform=ax.transAxes, fontsize=9,
            color="#1f6f62", va="top")
    ax.set_title(f"Staged complexes: binder length ({total} structures)")
    ax.set_xlabel("Binder residues")
    ax.set_ylabel("Structures")
    ax.legend(frameon=False, fontsize=8, loc="upper left",
              bbox_to_anchor=(1.01, 1.0))
    ax.grid(axis="y", **GRID_KW)

    # --- 2. prediction quality -------------------------------------------
    ax = axes[1]
    for n in names:
        xs, ys = [], []
        for r in batches[n]:
            a, b = as_float(r.get("iptm")), as_float(r.get("plddt"))
            if a is not None and b is not None:
                xs.append(a)
                ys.append(b)
        if xs:
            ax.scatter(xs, ys, s=26, alpha=0.75, color=colors[n],
                       edgecolor="white", linewidth=0.4, label=n, zorder=3)
    ax.axvline(0.8, color="#d94f4f", linestyle="--", linewidth=1.2)
    ax.text(0.805, ax.get_ylim()[0], " ipTM 0.8", fontsize=9, color="#a33",
            va="bottom")
    ax.set_title("Multimer prediction quality (all shown already passed "
                 "the campaign filter)")
    ax.set_xlabel("ipTM")
    ax.set_ylabel("pLDDT")
    ax.legend(frameon=False, fontsize=9)
    ax.grid(**GRID_KW)

    # --- 3. distinct backbones per length --------------------------------
    ax = axes[2]
    distinct = defaultdict(set)
    for n, rows in batches.items():
        for r in rows:
            bid = r.get("backbone_id_seq", "")
            if bid:
                # several designs off one backbone: 10022_6 -> 10022
                key = bid.rsplit("_", 1)[0]
            else:
                # no backbone column: curated sets (visually inspected,
                # cluster representatives) where each file is already
                # one independent fold
                key = r["file"]
            distinct[int(r["binder_len"])].add((n, key))
    ks = sorted(distinct)
    ax.bar(ks, [len(distinct[k]) for k in ks], color="#1f8a7a",
           edgecolor="white", linewidth=0.4, zorder=3)
    ax.set_title(f"Distinct backbones per length - "
                 f"{sum(len(v) for v in distinct.values())} independent "
                 f"folds among {total} designs")
    ax.set_xlabel("Binder residues")
    ax.set_ylabel("Distinct backbones")
    ax.set_xticks(ks)
    ax.tick_params(axis="x", labelsize=7, rotation=90)
    ax.grid(axis="y", **GRID_KW)
    ax.yaxis.set_major_locator(ticker.MaxNLocator(integer=True))

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
