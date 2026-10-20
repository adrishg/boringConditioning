#!/usr/bin/env python3
"""Cluster scaffolds structurally within each length and keep representatives.

Fold conditioning takes the design length from the scaffold's ss/adj tensors,
so a scaffold is only ever interchangeable with another of the *same* length.
Clustering therefore runs per (topology, length) group, never across lengths.

Similarity is foldseek TM-score rather than index-wise RMSD: two 67-residue
bundles can place their helix boundaries differently, and position-by-position
RMSD would call that a big difference when the fold is the same. Foldseek
aligns first, so it compares folds instead of registers.

Clustering is greedy by centrality. Within a group each structure is scored by
its summed similarity to the rest; the most central unassigned structure
becomes a seed, everything within --tm-cutoff of it joins that cluster, and
the loop repeats. The seed is the representative, so representatives are
typical members rather than outliers.

Usage:
    python3 cluster_scaffolds.py --tsv inventory/inventory.tsv --pool pool \
        --out-dir library --clusters inventory/clusters.tsv --tm-cutoff 0.5
"""
import argparse
import csv
import os
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict

USABLE = ("canonical", "splayed")
SENSITIVITY = (0.4, 0.5, 0.6, 0.7, 0.8, 0.9)


def find_foldseek(explicit=None):
    if explicit:
        return explicit
    hit = shutil.which("foldseek")
    if hit:
        return hit
    guess = os.path.expanduser("~/miniforge3/envs/foldseek/bin/foldseek")
    return guess if os.path.exists(guess) else None


def all_vs_all(foldseek, pdbs, workdir):
    """Return {(a, b): tm} for every pair, keyed by file stem."""
    indir = os.path.join(workdir, "in")
    os.makedirs(indir, exist_ok=True)
    for p in pdbs:
        os.symlink(os.path.abspath(p), os.path.join(indir, os.path.basename(p)))
    aln = os.path.join(workdir, "aln.tsv")
    cmd = [foldseek, "easy-search", indir, indir, aln,
           os.path.join(workdir, "tmp"),
           "--format-output", "query,target,alntmscore",
           "-e", "1e5", "--exhaustive-search", "1", "-v", "1"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0 or not os.path.exists(aln):
        raise RuntimeError(f"foldseek failed: {res.stderr[-300:]}")
    sim = {}
    with open(aln) as fh:
        for line in fh:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 3:
                continue
            q, t, tm = parts[0], parts[1], parts[2]
            try:
                # foldseek's alntmscore can read slightly above 1
                v = min(float(tm), 1.0)
            except ValueError:
                continue
            key = (q, t)
            sim[key] = max(sim.get(key, 0.0), v)
    return sim


def tm(sim, a, b):
    if a == b:
        return 1.0
    # the matrix is not perfectly symmetric; take the stronger direction
    return max(sim.get((a, b), 0.0), sim.get((b, a), 0.0))


def cluster(names, sim, cutoff):
    """Greedy by centrality. Returns [(representative, [members])]."""
    unassigned = set(names)
    out = []
    while unassigned:
        # most central remaining structure seeds the next cluster
        seed = max(unassigned,
                   key=lambda n: (sum(tm(sim, n, m) for m in unassigned), n))
        members = [m for m in unassigned if tm(sim, seed, m) >= cutoff]
        if seed not in members:
            members.append(seed)
        out.append((seed, sorted(members)))
        unassigned -= set(members)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tsv", required=True)
    ap.add_argument("--pool", required=True)
    ap.add_argument("--out-dir", required=True, help="representatives land here")
    ap.add_argument("--clusters", required=True, help="cluster table to write")
    ap.add_argument("--tm-cutoff", type=float, default=0.6,
                    help="TM-score at or above which two scaffolds are the "
                         "same cluster. Default 0.6, not the usual 0.5 "
                         "same-fold line: at 0.5 the ARD_alt 67-mers (TM "
                         "0.55-0.59 to the dominant fold) get swallowed, and the "
                         "greedy seeding is unstable -- adding one structure "
                         "at 85 res collapsed 5 clusters into 1. At 0.6 the "
                         "singleton set matches 0.7, so it is the point where "
                         "the answer stops moving.")
    ap.add_argument("--topology", type=int, action="append",
                    help="restrict to these helix counts (repeatable)")
    ap.add_argument("--include-all", action="store_true",
                    help="cluster every structure, not just usable bundles")
    ap.add_argument("--foldseek", help="path to the foldseek binary")
    ap.add_argument("--sensitivity", action="store_true",
                    help="also report cluster counts at other cutoffs")
    args = ap.parse_args()

    foldseek = find_foldseek(args.foldseek)
    if not foldseek:
        sys.exit("foldseek not found; pass --foldseek /path/to/foldseek")

    with open(args.tsv) as fh:
        rows = list(csv.DictReader(fh, delimiter="\t"))
    if not args.include_all:
        rows = [r for r in rows if r["bundle_class"] in USABLE]
    if args.topology:
        keep = set(args.topology)
        rows = [r for r in rows if int(r["n_helix"]) in keep]
    if not rows:
        sys.exit("no structures left after filtering")

    groups = defaultdict(list)
    for r in rows:
        groups[(int(r["n_helix"]), int(r["n_res"]))].append(r)

    by_stem = {os.path.splitext(os.path.basename(r["source_file"]))[0]: r
               for r in rows}
    os.makedirs(args.out_dir, exist_ok=True)
    out_rows, sens = [], defaultdict(int)

    print(f"{len(rows)} structures in {len(groups)} (topology, length) groups; "
          f"TM cutoff {args.tm_cutoff}")
    for (nh, nres), members in sorted(groups.items()):
        stems = [os.path.splitext(os.path.basename(m["source_file"]))[0]
                 for m in members]
        if len(members) == 1:
            clusters = [(stems[0], stems)]
            sim = {}
        else:
            with tempfile.TemporaryDirectory() as wd:
                paths = [os.path.join(args.pool, m["source_file"])
                         for m in members]
                sim = all_vs_all(foldseek, paths, wd)
            clusters = cluster(stems, sim, args.tm_cutoff)
            if args.sensitivity:
                for c in SENSITIVITY:
                    sens[(nh, nres, c)] = len(cluster(stems, sim, c))

        dest = os.path.join(args.out_dir, f"{nh}h_{nres:03d}res")
        os.makedirs(dest, exist_ok=True)
        for idx, (rep, mem) in enumerate(clusters, 1):
            src = os.path.join(args.pool, by_stem[rep]["source_file"])
            shutil.copy2(src, os.path.join(dest, f"{nh}h_{nres:03d}res_"
                                                 f"c{idx:02d}.pdb"))
            for m in mem:
                out_rows.append({
                    "clean_id": by_stem[m]["clean_id"],
                    "n_helix": nh, "n_res": nres,
                    "cluster": f"{nh}h_{nres:03d}res_c{idx:02d}",
                    "is_representative": int(m == rep),
                    "tm_to_rep": round(tm(sim, rep, m), 3),
                    "bundle_class": by_stem[m]["bundle_class"],
                    "campaign": by_stem[m]["campaign"],
                    "source_file": by_stem[m]["source_file"],
                })
        print(f"  {nh}h {nres:3d} res: {len(members):4d} -> "
              f"{len(clusters):3d} representatives")

    cols = ["clean_id", "n_helix", "n_res", "cluster", "is_representative",
            "tm_to_rep", "bundle_class", "campaign", "source_file"]
    os.makedirs(os.path.dirname(os.path.abspath(args.clusters)), exist_ok=True)
    with open(args.clusters, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(out_rows)

    n_rep = sum(1 for r in out_rows if r["is_representative"])
    print(f"\n{len(rows)} structures -> {n_rep} representatives "
          f"({args.out_dir})")
    print(f"cluster table: {args.clusters}")

    if args.sensitivity:
        print("\ncutoff sensitivity (total representatives):")
        for c in SENSITIVITY:
            total = sum(v for (_, _, cc), v in sens.items() if cc == c)
            singles = sum(1 for (nh, nres), m in groups.items() if len(m) == 1)
            print(f"  TM {c:.1f}: {total + singles}")


if __name__ == "__main__":
    main()
