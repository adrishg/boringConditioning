#!/usr/bin/env python3
"""Write the library's manifest and README so provenance travels with it.

Every scaffold in library/ carries a different pedigree: some passed AF
multimer filtering, some were curated by eye, and some are raw design decoys
accepted on purpose to close a length gap. That difference matters to whoever
uses the library next, and it must not live only in a report -- it goes in a
manifest sitting beside the structures.

Reads clusters.tsv (which structure represents which cluster), inventory.tsv
(topology and geometry) and batches.tsv (where each batch came from and how
it was validated), and writes library/MANIFEST.tsv plus library/README.md.

Usage:
    python3 make_library_manifest.py --clusters inventory/clusters.tsv \
        --inventory inventory/inventory.tsv --batches batches.tsv \
        --library library
"""
import argparse
import csv
import os
from collections import defaultdict

BAND = (58, 85)


def load(path, delim="\t"):
    with open(path) as fh:
        return list(csv.DictReader(fh, delimiter=delim))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--clusters", required=True)
    ap.add_argument("--inventory", required=True)
    ap.add_argument("--batches", required=True)
    ap.add_argument("--library", required=True)
    ap.add_argument("--band", nargs=2, type=int, default=list(BAND))
    ap.add_argument("--topology", type=int, default=4)
    args = ap.parse_args()
    lo, hi = args.band

    clusters = load(args.clusters)
    inv = {r["clean_id"]: r for r in load(args.inventory)}
    batches = {r["batch"]: r for r in load(args.batches)}

    size = defaultdict(int)
    for r in clusters:
        size[r["cluster"]] += 1

    rows = []
    for r in clusters:
        if r["is_representative"] != "1":
            continue
        batch = r["source_file"].split(os.sep)[0]
        meta = batches.get(batch, {})
        iv = inv.get(r["clean_id"], {})
        nh, nres = int(r["n_helix"]), int(r["n_res"])
        rows.append({
            "file": f"{nh}h_{nres:03d}res/{r['cluster']}.pdb",
            "n_helix": nh,
            "n_res": nres,
            "bundle_class": r["bundle_class"],
            "cluster_size": size[r["cluster"]],
            "helix_lens": iv.get("helix_lens", ""),
            "loop_lens": iv.get("loop_lens", ""),
            "validation": meta.get("validation", "unknown"),
            "batch": batch,
            "source": meta.get("source", ""),
            "original_file": r["source_file"],
        })
    rows.sort(key=lambda r: (r["n_helix"], r["n_res"], r["file"]))

    cols = ["file", "n_helix", "n_res", "bundle_class", "cluster_size",
            "helix_lens", "loop_lens", "validation", "batch", "source",
            "original_file"]
    man = os.path.join(args.library, "MANIFEST.tsv")
    with open(man, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    # ---- README ------------------------------------------------------
    t = args.topology
    per_len = defaultdict(list)
    for r in rows:
        if r["n_helix"] == t:
            per_len[r["n_res"]].append(r)
    floor = min(per_len) if per_len else 0
    gaps = [L for L in range(lo, hi + 1) if L not in per_len]
    unreachable = [L for L in gaps if L < floor]
    reachable = [L for L in gaps if L > floor]
    singles = sorted(L for L in per_len if len(per_len[L]) == 1)

    by_val = defaultdict(int)
    for r in rows:
        by_val[r["validation"]] += 1
    by_topo = defaultdict(int)
    for r in rows:
        by_topo[r["n_helix"]] += 1

    def row_line(L):
        rs = per_len[L]
        val = sorted({x["validation"] for x in rs})
        mark = " **" if len(rs) == 1 else ""
        return (f"| {L} | {len(rs)}{mark} | {', '.join(val)} | "
                f"{', '.join(sorted({x['batch'] for x in rs}))} |")

    readme = f"""# Scaffold library

{len(rows)} representative scaffolds for **RFdiffusion v1 fold conditioning**,
one per structural cluster, grouped as `<topology>h_<length>res/`.

Built by clustering {len(clusters)} usable bundles at foldseek TM 0.5, within
each (topology, length) group. `MANIFEST.tsv` carries the full record for each
file; everything below is derived from it.

| Topology | Representatives |
|---|---|
""" + "".join(f"| {k}-helix | {v} |\n" for k, v in sorted(by_topo.items())) + f"""

## Read this before using it

Scaffolds here do **not** all have the same pedigree:

| Validation | Scaffolds | Meaning |
|---|---|---|
""" + "".join(
        f"| {k} | {v} | "
        + {"multimer": "passed AF multimer filtering (ipTM/pLDDT)",
           "curated": "selected by hand or by eye, no score table",
           "unknown": "provenance not recorded",
           "UNVALIDATED": "**raw design decoys, never validated**"}.get(
               k, "") + " |\n"
        for k, v in sorted(by_val.items())) + f"""
Multimer scores say whether a design *bound its target*, which is not what a
scaffold is for -- fold conditioning reads only the secondary structure and
the block adjacency. A passing score is still a useful signal that the
backbone is a clean, foldable bundle. Unvalidated scaffolds were accepted
deliberately, to close length gaps that nothing else reaches.

## {t}-helix coverage, {lo}-{hi} res

Floor is **{floor} res** -- `sampled_insertion` only ever adds residues, so
nothing below this is reachable from the library.

| Residues | Representatives | Validation | Batch |
|---|---|---|---|
""" + "".join(row_line(L) + "\n" for L in sorted(per_len) if lo <= L <= hi) + f"""

`**` marks a length held up by a single structure.

- **Gaps (nothing at all):** {', '.join(map(str, gaps)) if gaps else 'none'}
  - reachable by insertion from the {floor}-res scaffold: {', '.join(map(str, reachable)) or 'none'}
  - below the floor, need new structures: {', '.join(map(str, unreachable)) or 'none'}
- **Single point of failure:** {', '.join(map(str, singles))}

## Using it

Convert on the cluster, inside the RFdiffusion checkout:

```bash
./helper_scripts/make_secstruc_adj.py --pdb_dir library/ --out_dir scaffolds/
```

That writes one `_ss.pt` and one `_adj.pt` per scaffold. Then:

```bash
./scripts/run_inference.py \\
  scaffoldguided.scaffoldguided=True \\
  scaffoldguided.scaffold_dir=scaffolds/ \\
  scaffoldguided.target_pdb=True \\
  scaffoldguided.target_path=input_pdbs/target.pdb \\
  scaffoldguided.target_ss=target_folds/target_ss.pt \\
  scaffoldguided.target_adj=target_folds/target_adj.pt \\
  'ppi.hotspot_res=[A59,A83,A91]' \\
  scaffoldguided.mask_loops=False \\
  denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
```

The target needs its own `_ss.pt` / `_adj.pt` from the same script.
`scaffoldguided.scaffold_list=path/to/list.txt` restricts a run to a subset --
use it to exclude the unvalidated scaffolds, or to pin a single length.

To reach the lengths between scaffolds, let the model insert residues:

```bash
scaffoldguided.mask_loops=True scaffoldguided.sampled_insertion=15 \\
  scaffoldguided.sampled_N=5 scaffoldguided.sampled_C=5
```

**RFdiffusion v1 only.** v2 does not document fold conditioning, and v3's
input specification answers "Can I guide on secondary structure?" with
"Currently no".
"""
    with open(os.path.join(args.library, "README.md"), "w") as fh:
        fh.write(readme)

    print(f"wrote {man} ({len(rows)} scaffolds)")
    print(f"wrote {os.path.join(args.library, 'README.md')}")
    for k, v in sorted(by_val.items()):
        print(f"  {k}: {v}")


if __name__ == "__main__":
    main()
