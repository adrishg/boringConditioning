#!/usr/bin/env python3
"""Inventory helical-bundle scaffold candidates for RFdiffusion fold conditioning.

For every PDB under --input-dir (recursively), this script:
  * drops byte-identical duplicates (keeps the first name alphabetically),
  * runs DSSP to assign secondary structure,
  * counts the helices and measures helix / loop lengths,
  * checks how bundle-like the packing is (helix-helix angles and distances),
  * assigns a clean, project-neutral ID of the form 4h_071res_03.

DSSP backend: mdtraj's Kabsch-Sander implementation by default, because it runs
in-process. The real mkdssp (--dssp mkdssp) reloads a 477 MB CCD dictionary on
every single call, which costs ~1.2 s per structure and does not scale as the
pool grows. The two agree residue-for-residue on these poly-glycine backbones
-- DSSP is pure backbone geometry and H-bonding, and the files carry N, CA, C
and O -- so the fast path loses nothing. Cross-check any time with:

    python3 inventory.py --input-dir pool --output /tmp/check.tsv --dssp mkdssp

Usage:
    python3 inventory.py --input-dir pool --output inventory/inventory.tsv
"""
import argparse
import csv
import glob
import hashlib
import os
import re
import shutil
import subprocess
import sys
import tempfile
import warnings

import numpy as np

warnings.filterwarnings("ignore")

# A run of H shorter than this is a turn, not a bundle helix.
MIN_HELIX = 8

# Helix pairs further apart than this (centroid-to-centroid) are not packed
# against each other. Measured: 97% of four-helix scaffolds here sit at
# 8-13 A, which is canonical helix-helix packing.
MAX_BUNDLE_DIST = 15.0

# Beyond this the helices are not forming a bundle at all. Deliberately loose:
# the crossing angle is bimodal across these campaigns -- a main mode at 10-39
# deg and a second, equally real cluster at 50-59 deg from well-packed but
# splayed bundles. An earlier 45 deg cutoff rejected hundreds of
# well-packed structures for having the wrong flavor of bundle. Fold
# conditioning only needs the secondary structure and which helices are
# adjacent, so splayed packing is still perfectly usable.
MAX_BUNDLE_ANGLE = 70.0

# Below this the bundle is the classic near-(anti)parallel up-down kind.
CANONICAL_ANGLE = 40.0


# --------------------------------------------------------------------------
# DSSP
# --------------------------------------------------------------------------

def find_cifpp_data():
    """mkdssp 4.x needs the CCD dictionary; conda installs it off PATH."""
    if os.environ.get("LIBCIFPP_DATA_DIR"):
        return os.environ["LIBCIFPP_DATA_DIR"]
    exe = shutil.which("mkdssp")
    if not exe:
        return None
    prefix = os.path.dirname(os.path.dirname(os.path.realpath(exe)))
    for cand in (os.path.join(prefix, "share", "libcifpp"),
                 "/usr/share/libcifpp", "/usr/local/share/libcifpp"):
        if os.path.exists(os.path.join(cand, "components.cif")):
            return cand
    return None


def dssp_mkdssp(pdb_path, data_dir):
    """Return the one-letter DSSP string, or None if mkdssp could not run."""
    env = dict(os.environ)
    if data_dir:
        env["LIBCIFPP_DATA_DIR"] = data_dir
    with tempfile.NamedTemporaryFile(suffix=".dssp", delete=False) as tmp:
        out_path = tmp.name
    try:
        # mkdssp sniffs mmCIF first and prints a harmless parse warning for
        # PDB input, so the return code and stderr are both unreliable here;
        # the only trustworthy signal is whether the output file got written.
        subprocess.run(
            ["mkdssp", "--output-format", "dssp", pdb_path, out_path],
            env=env, capture_output=True, timeout=120,
        )
        if not os.path.getsize(out_path):
            return None
        ss, started = [], False
        with open(out_path) as fh:
            for line in fh:
                if started:
                    if len(line) > 16 and line[13] != "!":
                        code = line[16]
                        ss.append(code if code != " " else "-")
                elif line.startswith("  #  RESIDUE"):
                    started = True
        return "".join(ss) if ss else None
    except (subprocess.SubprocessError, OSError):
        return None
    finally:
        if os.path.exists(out_path):
            os.unlink(out_path)


def dssp_mdtraj(traj):
    codes = md.compute_dssp(traj, simplified=False)[0]
    return "".join("-" if c in (" ", "NA", "") else c for c in codes)


# --------------------------------------------------------------------------
# geometry
# --------------------------------------------------------------------------

def segments(ss, codes="H"):
    """Contiguous runs of any character in `codes`, as (start, end_exclusive)."""
    out, start = [], None
    for i, c in enumerate(ss):
        if c in codes and start is None:
            start = i
        elif c not in codes and start is not None:
            out.append((start, i))
            start = None
    if start is not None:
        out.append((start, len(ss)))
    return out


def helix_axis(ca):
    """Unit vector along a helix: first principal component of its CA trace."""
    centered = ca - ca.mean(axis=0)
    axis = np.linalg.svd(centered, full_matrices=False)[2][0]
    if np.dot(ca[-1] - ca[0], axis) < 0:  # orient N -> C
        axis = -axis
    return axis


def md5(path):
    h = hashlib.md5()
    with open(path, "rb") as fh:
        for chunk in iter(lambda: fh.read(1 << 16), b""):
            h.update(chunk)
    return h.hexdigest()


def parse_name(fname):
    """Split a pooled filename into (campaign, backbone_id).

    Names look like campaign_tag_98__13206_22_unrelaxed_rank_001_alphafold2_..,
    project_name__796_13_81res, or 4hb_82res_2. The campaign sits left of the
    double underscore; the first number to its right is the RFdiffusion
    backbone the design came from. That id matters: several designs off one
    backbone are near-duplicate folds, so diversity filtering keys on it
    rather than on the filename.
    """
    stem = os.path.splitext(os.path.basename(fname))[0]
    stem = re.sub(r"_unrelaxed_rank_.*$", "", stem)  # ColabFold/AF2 model tag
    stem = re.sub(r"_\d+res(?=$|_)", "", stem)       # length tags from get_chain
    if "__" in stem:
        campaign, run = stem.split("__", 1)
        nums = re.findall(r"\d+", run)
        return campaign or "unnamed", (nums[0] if nums else "")
    return re.sub(r"_\d+$", "", stem) or "unnamed", ""


# --------------------------------------------------------------------------

def analyze(path, data_dir, min_helix=MIN_HELIX, backend="mdtraj",
            max_angle=MAX_BUNDLE_ANGLE, max_dist=MAX_BUNDLE_DIST):
    traj = md.load(path)
    ca_idx = traj.topology.select("name CA")
    xyz = traj.xyz[0][ca_idx] * 10.0  # nm -> angstrom

    ss = None
    if backend == "mkdssp":
        ss = dssp_mkdssp(path, data_dir)
    if ss is None or len(ss) != len(ca_idx):
        if backend == "mkdssp":
            backend = "mdtraj-fallback"
        ss = dssp_mdtraj(traj)

    helices = [(s, e) for s, e in segments(ss, "H") if e - s >= min_helix]
    hlens = [e - s for s, e in helices]
    loops = [helices[i + 1][0] - helices[i][1] for i in range(len(helices) - 1)]

    axes = [helix_axis(xyz[s:e]) for s, e in helices]
    centers = [xyz[s:e].mean(axis=0) for s, e in helices]

    angles, dists = [], []
    for i in range(len(axes) - 1):
        cosang = float(np.clip(np.dot(axes[i], axes[i + 1]), -1.0, 1.0))
        # abs(): a bundle is fine whether the pair runs up-down or up-up.
        angles.append(float(np.degrees(np.arccos(abs(cosang)))))
        dists.append(float(np.linalg.norm(centers[i] - centers[i + 1])))

    n = len(ss)
    worst_angle = max(angles) if angles else 0.0
    worst_dist = max(dists) if dists else 0.0
    packed = len(helices) >= 3 and worst_dist <= max_dist and ss.count("E") == 0
    is_bundle = packed and worst_angle <= max_angle
    if not packed:
        bundle_class = "not_bundle"
    elif worst_angle <= CANONICAL_ANGLE:
        bundle_class = "canonical"
    elif worst_angle <= max_angle:
        bundle_class = "splayed"
    else:
        bundle_class = "not_bundle"

    return dict(
        n_res=n,
        n_helix=len(helices),
        helix_lens=",".join(map(str, hlens)),
        helix_min=min(hlens) if hlens else 0,
        helix_max=max(hlens) if hlens else 0,
        loop_lens=",".join(map(str, loops)) or "-",
        loop_max=max(loops) if loops else 0,
        helix_frac=round(ss.count("H") / n, 3),
        strand_frac=round(ss.count("E") / n, 3),
        worst_pair_angle=round(worst_angle, 1),
        worst_pair_dist=round(worst_dist, 1),
        rg=round(float(md.compute_rg(traj.atom_slice(ca_idx))[0] * 10.0), 2),
        is_bundle=int(is_bundle),
        bundle_class=bundle_class,
        dssp_backend=backend,
        ss_string=ss,
    )


COLUMNS = [
    "clean_id", "n_res", "n_helix", "is_bundle", "bundle_class", "helix_lens", "helix_min",
    "helix_max", "loop_lens", "loop_max", "helix_frac", "strand_frac",
    "worst_pair_angle", "worst_pair_dist", "rg", "campaign", "backbone_id",
    "source_file", "md5", "dssp_backend", "ss_string",
]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--min-helix", type=int, default=MIN_HELIX,
                    help=f"shortest H run counted as a helix (default {MIN_HELIX})")
    ap.add_argument("--dssp", choices=("mdtraj", "mkdssp"), default="mdtraj",
                    help="DSSP backend; mkdssp is the reference but ~1.2 s/structure")
    ap.add_argument("--max-angle", type=float, default=MAX_BUNDLE_ANGLE,
                    help=f"worst helix-pair angle still called a bundle "
                         f"(default {MAX_BUNDLE_ANGLE:.0f} deg)")
    ap.add_argument("--max-dist", type=float, default=MAX_BUNDLE_DIST,
                    help=f"worst helix-pair separation still called packed "
                         f"(default {MAX_BUNDLE_DIST:.0f} A)")
    args = ap.parse_args()

    data_dir = find_cifpp_data()
    if args.dssp == "mkdssp":
        if not shutil.which("mkdssp"):
            sys.exit("--dssp mkdssp requested but mkdssp is not on PATH")
        if not data_dir:
            print("note: mkdssp found but its CCD dictionary was not; "
                  "falling back to mdtraj", file=sys.stderr)

    paths = sorted(glob.glob(os.path.join(args.input_dir, "**", "*.pdb"),
                             recursive=True))
    if not paths:
        sys.exit(f"no PDB files under {args.input_dir}")

    seen, rows, failed = {}, [], []
    for path in paths:
        digest = md5(path)
        if digest in seen:
            continue
        seen[digest] = path
        try:
            rec = analyze(path, data_dir, args.min_helix, args.dssp,
                          args.max_angle, args.max_dist)
        except Exception as exc:
            failed.append((os.path.basename(path), str(exc)[:70]))
            continue
        rec["source_file"] = os.path.relpath(path, args.input_dir)
        rec["campaign"], rec["backbone_id"] = parse_name(path)
        rec["md5"] = digest[:10]
        rows.append(rec)

    # Clean IDs: group by topology and length, number within the group.
    rows.sort(key=lambda r: (r["n_helix"], r["n_res"], r["source_file"]))
    counter = {}
    for r in rows:
        key = (r["n_helix"], r["n_res"])
        counter[key] = counter.get(key, 0) + 1
        r["clean_id"] = f"{r['n_helix']}h_{r['n_res']:03d}res_{counter[key]:02d}"

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=COLUMNS, delimiter="\t",
                           extrasaction="ignore")
        w.writeheader()
        w.writerows(rows)

    dupes = len(paths) - len(seen)
    print(f"{len(paths)} files scanned, {dupes} exact duplicates dropped, "
          f"{len(rows)} unique structures -> {args.output}")
    if failed:
        print(f"{len(failed)} failed:", file=sys.stderr)
        for name, err in failed[:10]:
            print(f"  {name}: {err}", file=sys.stderr)


if __name__ == "__main__":
    import mdtraj as md  # imported late so --help works without it
    globals()["md"] = md
    main()
