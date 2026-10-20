#!/usr/bin/env python3
"""Pull the binder chain out of design complexes into a flat pool folder.

Successor to the original getAchain.py. Same job -- across every campaign
checked in the lab notebook the binder is chain A and the target is chain B,
so a plain chain-A pull handles all of them -- with four fixes:

  * No length tag in the output name by default. The original always appended
    `_NNres`, which turned an already-single-chain folder into 418 files that
    were really 209 structures, each written twice under two names.
  * Keeps TER and END, so the output is a valid PDB rather than a bare pile
    of ATOM records.
  * Skips a file instead of silently writing nothing when the chain is absent,
    and says so.
  * Refuses to overwrite across batches, so pooling two sources that share a
    filename cannot lose a structure.

Usage:
    python3 get_chain.py --input-dir <source> --output-dir pool/<batch>
    python3 get_chain.py --input-dir <source> --output-dir pool/<batch> \
        --chain shortest --prefix vsd4_
"""
import argparse
import csv
import os
import sys
import warnings


def chains_from_cif(path):
    """mmCIF via mdtraj. Latent-X and AF3 write .cif, not .pdb, so the pool
    has to accept both. Returns {chain_id: (atom_indices, n_residues)}."""
    warnings.filterwarnings("ignore")
    import mdtraj as md
    traj = md.load(path)
    out = {}
    for ch in traj.topology.chains:
        cid = getattr(ch, "chain_id", None) or chr(ord("A") + ch.index)
        atoms = [a.index for a in ch.atoms]
        if atoms:
            out[cid] = (atoms, ch.n_residues)
    return traj, out


def write_pdb(traj, atom_indices, path, chain_id="A"):
    """Write selected atoms as a PDB.

    mdtraj's own save_pdb raises on these files (its atom-serial formatting
    assumes a serial the mmCIF does not carry), so the records are written
    here where the numbering is ours and renumbering from 1 is guaranteed.
    """
    xyz = traj.xyz[0] * 10.0  # nm -> angstrom
    top = traj.topology
    with open(path, "w") as fh:
        for serial, ai in enumerate(atom_indices, 1):
            a = top.atom(ai)
            x, y, z = xyz[ai]
            name = a.name
            # PDB column rule: short element names start one column in
            pad = f" {name:<3}" if len(name) < 4 else name
            el = (a.element.symbol if a.element is not None else "").upper()
            fh.write(
                f"ATOM  {serial % 100000:5d} {pad}{'':1}"
                f"{a.residue.name:>3} {chain_id}"
                f"{a.residue.resSeq % 10000:4d}{'':4}"
                f"{x:8.3f}{y:8.3f}{z:8.3f}"
                f"{1.00:6.2f}{0.00:6.2f}{'':10}{el:>2}\n")
        fh.write("TER\nEND\n")


def read_chains(pdb_path):
    """Map chain id -> (atom lines, ordered residue keys)."""
    chains, order = {}, {}
    with open(pdb_path) as fh:
        for line in fh:
            if not line.startswith(("ATOM", "HETATM")):
                continue
            chain = line[21]
            res_key = (line[22:27], line[17:20])  # resseq+icode, resname
            chains.setdefault(chain, []).append(line)
            seen = order.setdefault(chain, [])
            if not seen or seen[-1] != res_key:
                if res_key not in seen:
                    seen.append(res_key)
    return chains, order


def pick_chain(order, requested):
    if not order:
        return None
    if requested == "shortest":
        return min(order, key=lambda c: len(order[c]))
    if requested == "longest":
        return max(order, key=lambda c: len(order[c]))
    return requested if requested in order else None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--output-dir", required=True)
    ap.add_argument("--chain", default="A",
                    help="chain id to extract, or 'shortest' / 'longest' "
                         "(default A -- the binder in every campaign checked)")
    ap.add_argument("--prefix", default="",
                    help="prepended to every output name; use it to keep "
                         "batches apart when filenames collide")
    ap.add_argument("--tag-length", action="store_true",
                    help="append _NNres to output names (off by default -- "
                         "this is what duplicated the original set)")
    ap.add_argument("--min-res", type=int, default=0,
                    help="skip chains shorter than this")
    ap.add_argument("--max-res", type=int, default=0,
                    help="skip chains longer than this (0 = no limit)")
    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)
    rows, skipped = [], []

    names = sorted(f for f in os.listdir(args.input_dir)
                   if f.lower().endswith((".pdb", ".cif")))
    if not names:
        sys.exit(f"no PDB or mmCIF files in {args.input_dir}")

    for fname in names:
        in_path = os.path.join(args.input_dir, fname)
        is_cif = fname.lower().endswith(".cif")
        if is_cif:
            traj, cif_chains = chains_from_cif(in_path)
            order = {c: list(range(n)) for c, (_, n) in cif_chains.items()}
        else:
            chains, order = read_chains(in_path)
        chain = pick_chain(order, args.chain)
        if chain is None:
            skipped.append((fname, f"chain {args.chain} not found "
                                   f"(has: {','.join(sorted(order)) or 'none'})"))
            continue

        n_res = len(order[chain])
        if n_res < args.min_res:
            skipped.append((fname, f"{n_res} res < min {args.min_res}"))
            continue
        if args.max_res and n_res > args.max_res:
            skipped.append((fname, f"{n_res} res > max {args.max_res}"))
            continue

        base = os.path.splitext(fname)[0]
        out_name = f"{args.prefix}{base}"
        if args.tag_length:
            out_name += f"_{n_res}res"
        out_path = os.path.join(args.output_dir, out_name + ".pdb")
        if os.path.exists(out_path):
            skipped.append((fname, f"{os.path.basename(out_path)} already "
                                   f"exists -- pass a --prefix"))
            continue

        if is_cif:
            write_pdb(traj, cif_chains[chain][0], out_path, chain)
        else:
            with open(out_path, "w") as out:
                out.writelines(chains[chain])
                out.write("TER\nEND\n")

        rows.append({"output": os.path.basename(out_path), "source": fname,
                     "chain": chain, "n_res": n_res})

    csv_path = os.path.join(args.output_dir, "extraction_log.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=["output", "source", "chain", "n_res"])
        w.writeheader()
        w.writerows(rows)

    print(f"extracted {len(rows)} chains -> {args.output_dir}")
    if rows:
        lengths = sorted(r["n_res"] for r in rows)
        print(f"  lengths {lengths[0]}-{lengths[-1]} res")
    print(f"  log: {csv_path}")
    if skipped:
        print(f"skipped {len(skipped)}:", file=sys.stderr)
        for fname, why in skipped[:10]:
            print(f"  {fname}: {why}", file=sys.stderr)
        if len(skipped) > 10:
            print(f"  ... and {len(skipped) - 10} more", file=sys.stderr)


if __name__ == "__main__":
    main()
