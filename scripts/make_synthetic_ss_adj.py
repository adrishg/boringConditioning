#!/usr/bin/env python3
"""Build _ss.pt / _adj.pt for lengths no scaffold in the library reaches.

Unconditional generation plus a compactness potential does not produce
bundles: it controls radius of gyration, not helix count, and at 58 residues
it happily returns a two-helix hairpin. Fold conditioning does control
topology -- but it needs a scaffold at the target length, which is exactly
what is missing.

The way out is that fold conditioning never reads the scaffold's coordinates.
It reads two tensors:

    <name>_ss.pt   (L,)   per-residue secondary structure, 0=H 1=E 2=L 3=mask
    <name>_adj.pt  (L,L)  block adjacency, 1 where two secondary-structure
                          blocks touch

and RFdiffusion's own helper builds the adjacency *blockwise*: loops are
skipped entirely, and any helix pair with a Cb contact under 6 A has its
whole submatrix set to 1. Nothing else about the geometry survives. So the
tensors can be written directly, and a length with no PDB is not a barrier.

Rather than invent a topology, this takes a real four-helix bundle from the
library, reads its actual helix/loop segmentation and its measured adjacency
pattern, and rescales the segments to the target length. Every output is an
observed topology at a new size, not a guess. Using several source bundles
per length gives several distinct topologies.

Two steps, because the two halves need different environments: reading a
bundle needs mdtraj (workstation) and writing a .pt needs torch (SE3nv on
Hive), and neither environment has both.

    # 1. workstation -- measure the real topologies into a small JSON
    python3 make_synthetic_ss_adj.py --measure \
        --sources library/4h_060res/*.pdb library/4h_063res/*.pdb \
        --json topologies.json

    # 2. Hive -- turn those into tensors at the target lengths
    micromamba run -p .../SE3nv python make_synthetic_ss_adj.py \
        --json topologies.json --lengths 58,59 --out-dir synthetic_scaffolds
"""
import argparse
import json
import os
import sys

import numpy as np

# RFdiffusion's encoding (helper_scripts/make_secstruc_adj.py)
SS_HELIX, SS_STRAND, SS_LOOP, SS_MASK = 0, 1, 2, 3
CB_CUTOFF = 6.0
MIN_HELIX = 8


def cbeta(N, Ca, C):
    """Rosetta's virtual Cb, matching the helper script exactly."""
    b, c = Ca - N, C - Ca
    a = np.cross(b, c)
    return -0.57910144 * a + 0.5689693 * b - 0.5441217 * c + Ca


def segments_from_ss(ss, min_helix=MIN_HELIX):
    """Return [(kind, start, end)] with kind 'H' or 'L', merging short runs."""
    helices, start = [], None
    for i, ch in enumerate(ss):
        if ch == "H" and start is None:
            start = i
        elif ch != "H" and start is not None:
            if i - start >= min_helix:
                helices.append((start, i))
            start = None
    if start is not None and len(ss) - start >= min_helix:
        helices.append((start, len(ss)))

    segs, cursor = [], 0
    for s, e in helices:
        if s > cursor:
            segs.append(("L", cursor, s))
        segs.append(("H", s, e))
        cursor = e
    if cursor < len(ss):
        segs.append(("L", cursor, len(ss)))
    return segs


def measure_adjacency(pdb_path):
    """Read a bundle: its segment layout and which helix pairs touch."""
    import mdtraj as md
    traj = md.load(pdb_path)
    top = traj.topology
    xyz = traj.xyz[0] * 10.0
    idx_n = top.select("name N")
    idx_ca = top.select("name CA")
    idx_c = top.select("name C")
    if not (len(idx_n) == len(idx_ca) == len(idx_c)):
        raise ValueError(f"{pdb_path}: incomplete backbone")
    Cb = cbeta(xyz[idx_n], xyz[idx_ca], xyz[idx_c])

    ss = "".join(md.compute_dssp(traj, simplified=False)[0])
    segs = segments_from_ss(ss)
    helices = [(s, e) for kind, s, e in segs if kind == "H"]
    if len(helices) < 3:
        raise ValueError(f"{pdb_path}: only {len(helices)} helices")

    d = np.linalg.norm(Cb[:, None, :] - Cb[None, :, :], axis=-1)
    n = len(helices)
    touch = np.zeros((n, n), dtype=bool)
    for i in range(n):
        for j in range(i + 1, n):
            a, b = helices[i], helices[j]
            if (d[a[0]:a[1], b[0]:b[1]] < CB_CUTOFF).any():
                touch[i, j] = touch[j, i] = True
    return segs, touch


def rescale(segs, target):
    """Rescale segment lengths to sum to `target`, keeping loops short.

    Loops are held at their original length where possible -- they set the
    turn geometry and squeezing them distorts the fold. The difference is
    taken from the helices in proportion to their size, so a bundle shrinks
    evenly rather than losing one helix.
    """
    loops = [(k, e - s) for k, s, e in segs if k == "L"]
    helix_lens = [e - s for k, s, e in segs if k == "H"]
    loop_total = sum(n for _, n in loops)
    helix_total = sum(helix_lens)

    avail = target - loop_total
    if avail < MIN_HELIX * len(helix_lens):
        # loops must give way too, but never below 2 residues each
        need = MIN_HELIX * len(helix_lens)
        shrink = min(loop_total - 2 * len(loops), need - avail) if loops else 0
        loop_total -= max(shrink, 0)
        avail = target - loop_total
        if avail < MIN_HELIX * len(helix_lens):
            return None

    scaled = [max(MIN_HELIX, round(h * avail / helix_total)) for h in helix_lens]
    # fix rounding drift on the longest helix so the total lands exactly
    drift = avail - sum(scaled)
    while drift != 0:
        k = int(np.argmax(scaled)) if drift < 0 else int(np.argmin(scaled))
        if drift < 0 and scaled[k] <= MIN_HELIX:
            return None
        scaled[k] += 1 if drift > 0 else -1
        drift = avail - sum(scaled)

    # rebuild the segment list with the new sizes
    out, hi, li = [], 0, 0
    loop_sizes = [n for _, n in loops]
    if loop_total != sum(loop_sizes) and loop_sizes:
        per = loop_total // len(loop_sizes)
        loop_sizes = [max(2, per)] * len(loop_sizes)
        loop_sizes[-1] += loop_total - sum(loop_sizes)
    for kind, s, e in segs:
        if kind == "H":
            out.append(("H", scaled[hi])); hi += 1
        else:
            out.append(("L", loop_sizes[li])); li += 1
    total = sum(n for _, n in out)
    return out if total == target else None


def build_tensors(layout, touch):
    """Turn a rescaled layout plus a touch matrix into (ss, adj) arrays."""
    ss, spans, pos, hi = [], [], 0, 0
    for kind, n in layout:
        if kind == "H":
            ss.extend([SS_HELIX] * n)
            spans.append((pos, pos + n))
            hi += 1
        else:
            ss.extend([SS_LOOP] * n)
        pos += n
    L = len(ss)
    adj = np.zeros((L, L), dtype=np.float32)
    for i in range(len(spans)):
        for j in range(i + 1, len(spans)):
            if touch[i, j]:
                a, b = spans[i], spans[j]
                adj[a[0]:a[1], b[0]:b[1]] = 1.0
                adj[b[0]:b[1], a[0]:a[1]] = 1.0
    return np.array(ss, dtype=np.float32), adj


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--measure", action="store_true",
                    help="step 1: read PDBs with mdtraj and write --json")
    ap.add_argument("--sources", nargs="+", help="real bundle PDBs (step 1)")
    ap.add_argument("--json", required=True,
                    help="topology file written by step 1, read by step 2")
    ap.add_argument("--lengths", help="e.g. 58,59 (step 2)")
    ap.add_argument("--out-dir", help="where the .pt files go (step 2)")
    args = ap.parse_args()

    if args.measure:
        if not args.sources:
            sys.exit("--measure needs --sources")
        out = []
        for src in args.sources:
            tag = os.path.splitext(os.path.basename(src))[0]
            try:
                segs, touch = measure_adjacency(src)
            except Exception as exc:
                print(f"  skip {tag}: {exc}", file=sys.stderr)
                continue
            out.append({"tag": tag,
                        "segs": [[k, int(s), int(e)] for k, s, e in segs],
                        "touch": touch.astype(int).tolist()})
            print(f"  measured {tag}: "
                  f"{sum(1 for k, _, _ in segs if k == 'H')} helices")
        with open(args.json, "w") as fh:
            json.dump(out, fh, indent=1)
        print(f"wrote {len(out)} topologies to {args.json}")
        return

    if not (args.lengths and args.out_dir):
        sys.exit("step 2 needs --lengths and --out-dir")
    try:
        import torch
    except ImportError:
        sys.exit("torch is required; run inside the RFdiffusion env (SE3nv)")

    with open(args.json) as fh:
        topologies = json.load(fh)

    lengths = [int(x) for x in args.lengths.split(",")]
    os.makedirs(args.out_dir, exist_ok=True)
    made = 0
    for rec in topologies:
        tag = rec["tag"]
        segs = [(k, s, e) for k, s, e in rec["segs"]]
        touch = np.array(rec["touch"], dtype=bool)
        n_helix = int(sum(1 for k, _, _ in segs if k == "H"))
        for L in lengths:
            layout = rescale(segs, L)
            if layout is None:
                print(f"  skip {tag} @ {L}: cannot fit {n_helix} helices",
                      file=sys.stderr)
                continue
            ss, adj = build_tensors(layout, touch)
            name = f"syn_{n_helix}h_{L:03d}res_{tag}"
            torch.save(torch.from_numpy(ss), os.path.join(args.out_dir, f"{name}_ss.pt"))
            torch.save(torch.from_numpy(adj), os.path.join(args.out_dir, f"{name}_adj.pt"))
            helices = [n for k, n in layout if k == "H"]
            print(f"  {name}: helices {helices} loops "
                  f"{[n for k, n in layout if k == 'L']} total {int(ss.shape[0])}")
            made += 1
    print(f"wrote {made} scaffold tensor pairs to {args.out_dir}")


if __name__ == "__main__":
    main()
