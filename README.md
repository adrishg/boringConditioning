# boringConditioning

![Helices are all you need](figures/banner_misscav_bell_curve.png)

> Banner adapted from a meme posted by **Martin Pacesa** on LinkedIn.
> <!-- TODO: add the post URL -->

Build a library of small helical-bundle scaffolds for **RFdiffusion fold
conditioning**, from design structures you already have.

Instead of letting RFdiffusion pick topology freely and discarding everything
that is not a 4-helix bundle, condition it up front so every design is one.
"Boring" on purpose: predictable, well-behaved folds.

## What this repository is

**Tools, not structures.** No PDBs and no `.pt` tensors are distributed here —
they are unpublished design data. What ships is the pipeline that turns a pile
of design outputs into a deduplicated, length-indexed scaffold library.

Point it at your own structures and it builds your own library.

## What our library looks like

For scale, the library this pipeline produced from our own design campaigns —
**364 distinct folds**, one per structural cluster, every length from 54 to 102
residues:

| Topology | Folds |
|---|---|
| 3-helix | 68 |
| **4-helix** | **179** |
| 5-helix | 103 |
| 6-helix | 14 |

4-helix bundles per length, 54–102 res — no gaps from 58 up:

```
54:1   57:1   58:5   59:4   60:5   61:6   62:6   63:4   64:6   65:5
66:6   67:14  68:3   69:9   70:6   71:4   72:2   73:2   74:5   75:1
76:3   77:1   78:3   79:4   80:2   81:3   82:5   83:2   84:2   85:3
86:2   87:2   88:1   89:3   90:3   91:2   92:5   93:2   94:1   95:3
96:6   97:7   98:12  99:3   100:4
```

Those 364 folds came from 2,177 usable structures — the rest were redundant.

## Requirements

```bash
conda create -n boring python=3.11 -c conda-forge mdtraj numpy matplotlib
conda activate boring
```

Plus, for the clustering and conversion steps:

- [foldseek](https://github.com/steineggerlab/foldseek) — structural clustering
- [RFdiffusion](https://github.com/RosettaCommons/RFdiffusion) — its
  `helper_scripts/make_secstruc_adj.py` writes the tensors, and its `SE3nv`
  environment runs them

`mkdssp` is optional; secondary structure defaults to mdtraj's Kabsch–Sander,
which agrees with mkdssp 4.6 residue-for-residue on these backbones.

## Quick start

```bash
git clone <this repo> boringConditioning
cd boringConditioning
conda activate boring
```

### 1. Harvest the designs that passed your filters

Most design pipelines rename files after scoring, so the score table rarely
matches what is on disk. `--name-template` rebuilds the disk name from CSV
columns; use `--dry-run` first to check it resolves.

```bash
python3 scripts/harvest_complexes.py \
    --source-dir  <campaign>/multimer_results/models \
    --score-csv   <campaign>/multimer_results/multimerTest_filtered.csv \
    --name-template 'design__{backbone_id_seq}__multimer-af2__seed-{seed}__model-{model_num}__rank-{rank}.pdb' \
    --output-dir  complexes/<batch> --binder-chain A --require-passed \
    --prefix '<batch>__' --dry-run
```

For an already-curated folder with no score table, drop `--score-csv` and
filter on length instead.

### 2. Extract the binder chain

```bash
python3 scripts/get_chain.py --input-dir complexes/<batch> \
    --output-dir pool/<batch> --chain A --min-res 40 --max-res 130
```

Reads PDB and mmCIF. Across every campaign we checked, the binder is chain A
and the target is chain B.

### 3. Inventory — topology, geometry, redundancy

```bash
python3 scripts/inventory.py --input-dir pool --output inventory/inventory.tsv
```

Deduplicates by md5, runs DSSP, counts helices, measures helix/loop lengths
and packing, and classifies each structure `canonical` / `splayed` /
`not_bundle`. Splayed bundles are usable: fold conditioning reads only the
secondary structure and which helices touch.

### 4. Cluster and keep one representative per fold

```bash
python3 scripts/cluster_scaffolds.py --tsv inventory/inventory.tsv \
    --pool pool --out-dir library --clusters inventory/clusters.tsv
```

Clusters **within** each (topology, length) group — a scaffold is only
interchangeable with another of the same length, because fold conditioning
takes the design length from the tensors. Similarity is foldseek TM-score, not
index-wise RMSD, since two same-length bundles can place their helix
boundaries differently.

This is where the real size of a library shows up. Ours went from 2,177 usable
structures to **374 distinct folds**; one 180-design set collapsed to a single
fold.

### 5. Record provenance

```bash
python3 scripts/make_library_manifest.py --clusters inventory/clusters.tsv \
    --inventory inventory/inventory.tsv --batches batches.tsv --library library
```

Writes `library/MANIFEST.tsv` with each scaffold's topology, cluster size and
**validation status**, so a scaffold that never passed a filter is labelled as
such rather than silently mixed in. Describe each source in `batches.tsv`.

### 6. Convert to the tensors RFdiffusion reads

On the cluster, inside the RFdiffusion environment:

```bash
sbatch cluster/1-make_secstruc_adj.sbatch [LIBRARY_DIR] [OUT_DIR]
```

One `_ss.pt` and one `_adj.pt` per scaffold. The script flattens the nested
library, refuses to run on a filename collision, and fails if the output count
does not match the input count — a partial conversion would otherwise leave
you sampling from a smaller set than you think.

### 7. Run fold conditioning

```bash
$RFDIFFUSION/scripts/run_inference.py \
    inference.output_prefix=out/design \
    scaffoldguided.scaffoldguided=True \
    scaffoldguided.scaffold_dir=scaffolds/ \
    scaffoldguided.target_pdb=True \
    scaffoldguided.target_path=input_pdbs/target.pdb \
    scaffoldguided.target_ss=target_folds/target_ss.pt \
    scaffoldguided.target_adj=target_folds/target_adj.pt \
    'ppi.hotspot_res=[A59,A83,A91]' \
    scaffoldguided.mask_loops=False \
    denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
```

The **target needs its own** `_ss.pt` / `_adj.pt` from the same helper.
`scaffoldguided.scaffold_list=list.txt` restricts a run to a subset — use it to
exclude unvalidated scaffolds, or to pin one length.

## Looking at what you built

Optional, all reading the tables from the steps above:

```bash
python3 scripts/plot_inventory.py   --tsv inventory/inventory.tsv \
    --output coverage.png --band 58 85 --topology 4 --bundles-only
python3 scripts/plot_library.py     --clusters inventory/clusters.tsv \
    --output library.png --band 58 85 --topology 4
python3 scripts/plot_quality.py     --tsv inventory/inventory.tsv --output quality.png
python3 scripts/plot_complexes.py   --complexes-dir complexes --output complexes.png
```

`plot_library.py` is the useful one: a topology-by-length grid of what you
have, and how much of each group was redundant. `plot_source_survey.py` answers
a different question — *should I harvest this at all?* — from a length sweep of
a candidate source, before copying anything.

## Checking whether unvalidated scaffolds are real

Scaffolds that never passed a filter are still worth keeping if they are the
only fold at their length, but it is worth knowing whether they are
*designable*: does a sequence designed onto the backbone fold back to it?

```bash
bash cluster/3-revise-unvalidated.sh [SCAFFOLD_DIR] [OUT_DIR] [NUM_SEQS]
```

A submit wrapper, not a job — it chains sequence design and a monomer
foldability test with an `afterok` dependency. It expects a binder-design
pipeline providing those two steps; set `BINDER_REPO_DIR` to it. Run it on
cluster representatives only, so nothing is spent re-checking folds that
clustering already called redundant.

## Filling lengths you do not have

Fold conditioning sets the design length from the scaffold, so a length with no
scaffold is a length you cannot design at. Two ways out.

**Insertion sampling**, when the gap is *above* your shortest scaffold:

```bash
scaffoldguided.mask_loops=True scaffoldguided.sampled_insertion=15 \
    scaffoldguided.sampled_N=5 scaffoldguided.sampled_C=5
```

**Synthetic conditioning tensors**, when it is below the floor and nothing
exists. Fold conditioning never reads a scaffold's coordinates, only its
secondary-structure string and block-adjacency matrix — so a real bundle's
topology can be rescaled to a length that has no PDB:

```bash
# workstation (needs mdtraj)
python3 scripts/make_synthetic_ss_adj.py --measure \
    --sources library/4h_060res/*.pdb --json topologies.json
# cluster, inside SE3nv (needs torch)
python3 scripts/make_synthetic_ss_adj.py --json topologies.json \
    --lengths 58,59 --out-dir synthetic_scaffolds
sbatch cluster/2-design_missing_lengths.sbatch
```

The two steps are split because reading a structure needs mdtraj and writing a
`.pt` needs torch, and the two environments rarely have both.

Expect roughly a third of designs to be usable at the shortest lengths, rising
with length — at 58 residues there is barely room for four helices and many
collapse into a two-helix hairpin. Unconditional generation with a compactness
potential does **not** work here: it controls radius of gyration, not helix
count.

## Version support

| Version | Fold conditioning |
|---|---|
| **RFdiffusion v1** | **Yes** — the only version this library works with |
| RFdiffusion2 | No; its README documents motif scaffolding only |
| RFdiffusion3 | No. Its input spec answers "Can I guide on secondary structure?" with "Currently no" |

## Repository layout

```
scripts/    the pipeline, step by step (10 tools)
cluster/    SLURM jobs (conversion, gap-filling design, revision)
docs/       NOTES.md — what cost time when building one of these
batches.tsv template for recording where each source batch came from
figures/    banner and coverage plots
```

`docs/NOTES.md` records the failure modes — mismatched score tables, corrupt
coordinates that parse cleanly, clustering thresholds, why unconditional
generation gives hairpins — so they do not have to be rediscovered.
