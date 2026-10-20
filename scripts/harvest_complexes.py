#!/usr/bin/env python3
"""Copy design complexes that passed multimer testing into the staging area.

Copies, never moves -- the source campaigns stay untouched. Binders are
extracted later with get_chain.py, so complexes are kept whole here.

The awkward part is that a campaign's score CSV rarely names files the way
they sit on disk; the pipeline renames them after scoring. So the filename is
rebuilt from CSV columns with --name-template, e.g.

    case A csv: design__10022_6_unrelaxed_rank_001_alphafold2_multimer_v3_...
           disk: design__10022_6__multimer-af2__seed-001__model-2__rank-001.pdb
           --name-template 'design__{backbone_id_seq}__multimer-af2__seed-{seed}__model-{model_num}__rank-{rank}.pdb'

    case B --name-template '{design_folder}_af3_model_{model}.pdb'

When the CSV's own column already matches disk, use --file-column instead.
Every row that cannot be resolved is reported rather than silently dropped.

Usage:
    python3 harvest_complexes.py \
        --source-dir <campaign>/models --score-csv <campaign>/filtered.csv \
        --name-template '...' --output-dir complexes/<batch> \
        --binder-chain A --require-passed
"""
import argparse
import csv
import os
import shutil
import sys


def ca_count(pdb_path, chain):
    n = 0
    with open(pdb_path) as fh:
        for line in fh:
            if line.startswith("ATOM") and line[21] == chain and line[12:16] == " CA ":
                n += 1
    return n


def as_float(value):
    try:
        return float(value)
    except (TypeError, ValueError):
        return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--source-dir", required=True, help="folder holding the PDBs")
    ap.add_argument("--score-csv",
                    help="score table; omit for an already-curated folder "
                         "(visually inspected sets, cluster representatives) "
                         "where every PDB is taken and only length filters apply")
    ap.add_argument("--output-dir", required=True)
    group = ap.add_mutually_exclusive_group()
    group.add_argument("--name-template",
                       help="disk filename built from CSV columns, e.g. "
                            "'{design_folder}_af3_model_{model}.pdb'")
    group.add_argument("--file-column",
                       help="CSV column that already holds the disk filename")
    ap.add_argument("--binder-chain", default="A",
                    help="chain that is the binder (default A); recorded in "
                         "the manifest so extraction later is unambiguous")
    ap.add_argument("--prefix", default="", help="prepended to copied names")
    ap.add_argument("--require-passed", action="store_true",
                    help="keep only rows whose 'passed' column is true")
    ap.add_argument("--min-iptm", type=float)
    ap.add_argument("--min-plddt", type=float)
    ap.add_argument("--max-rmsd", type=float, help="upper bound on rmsd_A")
    ap.add_argument("--min-len", type=int, help="skip binders shorter than this")
    ap.add_argument("--max-len", type=int,
                    help="skip binders longer than this; also catches files "
                         "whose chains are mislabeled and hand back a whole target")
    ap.add_argument("--dry-run", action="store_true",
                    help="report what would be copied, write nothing")
    args = ap.parse_args()

    if args.score_csv:
        if not (args.name_template or args.file_column):
            sys.exit("--score-csv needs --name-template or --file-column")
        with open(args.score_csv) as fh:
            rows = list(csv.DictReader(fh))
        if not rows:
            sys.exit(f"{args.score_csv} has no rows")
    else:
        if args.require_passed or args.min_iptm or args.min_plddt or args.max_rmsd:
            sys.exit("score filters need --score-csv")
        rows = [{"__file__": f} for f in
                sorted(os.listdir(args.source_dir)) if f.endswith(".pdb")]
        if not rows:
            sys.exit(f"no PDB files in {args.source_dir}")
        args.file_column = "__file__"

    kept, dropped, missing = [], [], []
    for row in rows:
        if args.require_passed:
            if str(row.get("passed", "")).strip().lower() not in ("true", "1", "yes"):
                dropped.append("not passed")
                continue
        iptm, plddt = as_float(row.get("iptm")), as_float(row.get("plddt"))
        rmsd = as_float(row.get("rmsd_A"))
        if args.min_iptm is not None and (iptm is None or iptm < args.min_iptm):
            dropped.append("iptm")
            continue
        if args.min_plddt is not None and (plddt is None or plddt < args.min_plddt):
            dropped.append("plddt")
            continue
        if args.max_rmsd is not None and (rmsd is None or rmsd > args.max_rmsd):
            dropped.append("rmsd_A")
            continue

        if args.file_column:
            name = os.path.basename(row.get(args.file_column, ""))
        else:
            try:
                name = args.name_template.format(**row)
            except KeyError as exc:
                sys.exit(f"--name-template refers to column {exc} which the "
                         f"CSV does not have. Columns: {', '.join(rows[0])}")
        src = os.path.join(args.source_dir, name)
        if not os.path.isfile(src):
            missing.append(name)
            continue
        if args.min_len or args.max_len:
            blen = ca_count(src, args.binder_chain)
            if args.min_len and blen < args.min_len:
                dropped.append("binder too short")
                continue
            if args.max_len and blen > args.max_len:
                dropped.append("binder too long")
                continue
        kept.append((src, name, row, iptm, plddt, rmsd))

    if missing:
        print(f"WARNING: {len(missing)} scored rows have no file on disk, "
              f"e.g. {missing[0]}", file=sys.stderr)

    if args.dry_run:
        print(f"would copy {len(kept)} of {len(rows)} rows "
              f"({len(dropped)} filtered out, {len(missing)} missing)")
        return

    os.makedirs(args.output_dir, exist_ok=True)
    manifest, collisions = [], 0
    for src, name, row, iptm, plddt, rmsd in kept:
        out_name = args.prefix + name
        dst = os.path.join(args.output_dir, out_name)
        if os.path.exists(dst):
            collisions += 1
            continue
        shutil.copy2(src, dst)
        manifest.append({
            "file": out_name,
            "binder_chain": args.binder_chain,
            "binder_len": ca_count(dst, args.binder_chain),
            "iptm": row.get("iptm", ""),
            "plddt": row.get("plddt", ""),
            "rmsd_A": row.get("rmsd_A", ""),
            "rmsd_B": row.get("rmsd_B", ""),
            "backbone_id_seq": row.get("backbone_id_seq", ""),
            "source_dir": os.path.abspath(args.source_dir),
            "source_file": name,
        })

    csv_path = os.path.join(args.output_dir, "manifest.csv")
    with open(csv_path, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(manifest[0]) if manifest else
                           ["file", "binder_chain", "binder_len"])
        w.writeheader()
        w.writerows(manifest)

    print(f"copied {len(manifest)} complexes -> {args.output_dir}")
    if manifest:
        lens = sorted(m["binder_len"] for m in manifest)
        print(f"  binder chain {args.binder_chain}: {lens[0]}-{lens[-1]} res")
    if collisions:
        print(f"  {collisions} skipped, name already present (use --prefix)")
    print(f"  manifest: {csv_path}")


if __name__ == "__main__":
    main()
