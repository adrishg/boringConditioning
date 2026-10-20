#!/usr/bin/env python3
"""Plot what binder lengths a surveyed source tree actually contains.

Used to answer "is it worth harvesting?" before copying anything. Takes the
TSV produced by a length sweep (columns binder_len, file) and draws the
length distribution against the target band, marking which lengths in the
band the source covers and which it leaves empty.

Usage:
    python3 plot_source_survey.py --tsv inventory/pore_survey_lengths.tsv \
        --output inventory/pore_survey.png --band 58 85 \
        --title "a source tree"
"""
import argparse
import csv
from collections import Counter

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

PRESENT = "#1f8a7a"
EMPTY = "#d94f4f"
OUT_OF_BAND = "#b8ccd4"
GRID_KW = dict(axis="y", linestyle="--", alpha=0.45, zorder=0)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--band", nargs=2, type=int, default=[58, 85],
                    metavar=("LO", "HI"))
    ap.add_argument("--title", default="Source survey")
    ap.add_argument("--min-len", type=int, default=40,
                    help="ignore shorter chains; they are parsing artifacts, "
                         "not binders (default 40)")
    ap.add_argument("--max-len", type=int, default=115,
                    help="ignore longer chains; a single mislabeled file with "
                         "a whole target in chain A otherwise stretches the "
                         "axis until nothing is readable (default 115)")
    args = ap.parse_args()

    with open(args.tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    lengths = [int(r["binder_len"]) for r in rows]
    short = sum(1 for L in lengths if L < args.min_len)
    longs = sorted(L for L in lengths if L > args.max_len)
    counts = Counter(L for L in lengths
                     if args.min_len <= L <= args.max_len)
    if not counts:
        raise SystemExit("nothing left after --min-len")

    lo, hi = args.band
    span = list(range(min(min(counts), lo), max(max(counts), hi) + 1))
    heights = [counts.get(L, 0) for L in span]
    colors = []
    for L, h in zip(span, heights):
        if not (lo <= L <= hi):
            colors.append(OUT_OF_BAND)
        else:
            colors.append(PRESENT if h else EMPTY)

    fig, ax = plt.subplots(figsize=(13, 5.5))
    ax.axvspan(lo - 0.5, hi + 0.5, color="#1f8a7a", alpha=0.07, zorder=1)
    ax.bar(span, [h if h else 0.9 for h in heights], color=colors,
           edgecolor="white", linewidth=0.4, zorder=3)
    for L, h in zip(span, heights):
        if h:
            ax.text(L, h * 1.15, str(h), ha="center", va="bottom",
                    fontsize=8)

    empty = [L for L in range(lo, hi + 1) if not counts.get(L)]
    covered = [L for L in range(lo, hi + 1) if counts.get(L)]
    ax.set_yscale("log")
    ax.set_title(
        f"{args.title}: {sum(counts.values())} binders sampled, "
        f"only {len(counts)} distinct lengths exist\n"
        f"in band {lo}-{hi} -> covers {covered or 'nothing'};  "
        f"empty at {', '.join(map(str, empty)) if empty else 'nothing'}")
    ax.set_xlabel("Binder chain length (residues)")
    ax.set_ylabel("Files sampled (log)")
    ax.set_xticks([L for L in span if L % 5 == 0] +
                  [L for L in span if counts.get(L)])
    ax.tick_params(axis="x", labelsize=8, rotation=90)
    ax.grid(**GRID_KW)

    handles = [plt.Rectangle((0, 0), 1, 1, color=c) for c in
               (PRESENT, EMPTY, OUT_OF_BAND)]
    ax.legend(handles, ["in band, present", "in band, EMPTY",
                        "outside band"], frameon=False, loc="upper left")

    notes = []
    if short:
        notes.append(f"{short} file(s) with no chain A")
    if longs:
        notes.append(f"{len(longs)} over {args.max_len} res "
                     f"({', '.join(map(str, longs))}) - mislabeled chains")
    if notes:
        ax.text(0.99, 0.94, "off-plot: " + "; ".join(notes),
                transform=ax.transAxes, ha="right", fontsize=8, color="#777")

    fig.tight_layout()
    fig.savefig(args.output, dpi=200)
    print(f"wrote {args.output}")


if __name__ == "__main__":
    main()
