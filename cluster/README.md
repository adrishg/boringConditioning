# Hive cluster jobs

Two SLURM scripts. They follow the `sacredBinderPipeline` conventions:
partition, account and mail-user are **passed at submit time**, not baked in,
and the RFdiffusion environment is activated through
`${BINDER_REPO_DIR}/bash_tools/load_env.sh` when it is available, falling back
to activating `SE3nv` by name and then by prefix
(`$SE3NV_PREFIX`).

Verified on Hive 2026-10-12: RFdiffusion is at
`$RFDIFFUSION_ROOT` and has `make_secstruc_adj.py`.
The `/share/yarovlab/...` path in the pipeline repo's `env.defaults.sh` does
**not** exist on this cluster.

Defaults assumed, override by exporting before submitting:

```bash
export BINDER_RFDIFFUSION_ROOT=$RFDIFFUSION_ROOT
export BINDER_ENV_BACKBONE=SE3nv
```

## 1. Convert the library to ss/adj tensors

This is the step that has never been run — **no `.pt` file exists anywhere
yet**, so nothing in the library is usable by RFdiffusion until this finishes.

```bash
sbatch --partition="$SLURM_CPU_PARTITION" --mail-user="$SLURM_MAIL_USER" \
       1-make_secstruc_adj.sbatch /path/to/library /path/to/scaffolds
```

CPU only, ~4 h wall limit is generous for 267 scaffolds. The script flattens
the nested `library/<topology>h_<len>res/` tree into symlinks first, because
`make_secstruc_adj.py` takes a single `--pdb_dir`. It refuses to run if two
scaffolds share a basename, and fails if the count of `_ss.pt`/`_adj.pt` does
not match the number of inputs — a partial conversion would otherwise leave a
scaffold set quietly smaller than intended.

The **target** needs its own `_ss.pt`/`_adj.pt` from the same helper:

```bash
python "$BINDER_RFDIFFUSION_ROOT/helper_scripts/make_secstruc_adj.py" \
    --input_pdb target.pdb --out_dir target_folds/
```

## 2. Design the lengths that do not exist

A sweep of all 62,806 structures in the lab notebook found **zero** binders at
58, 59 or 64 residues, and 66 only as 3-helix. Those cannot be recovered from
anything on disk. 60, 61 and 65 exist but each rests on a single unvalidated
structure, so they are included to give them independent folds.

```bash
sbatch --partition="$SLURM_GPU_PARTITION" --mail-user="$SLURM_MAIL_USER" \
       --account="$SLURM_GPU_ACCOUNT" \
       2-design_missing_lengths.sbatch /path/to/gapDesigns 200
```

A 7-task array, one length each: 58, 59, 60, 61, 64, 65, 66.

This is **unconditional monomer generation**, not binder design — a
fold-conditioning scaffold only supplies topology, so there is no target and
no hotspot. `monomer_ROG` pushes toward compact folds, which is what yields
bundles instead of extended chains. Topology is not guaranteed; the output is
filtered afterwards. At 200 designs per length, expect to keep roughly a
quarter to a half as 4-helix bundles, which is ample.

Noise is left at its default rather than zeroed as in the binder-design
script: this run wants fold diversity, not one best pose.

## 3. Bring the results back

```bash
# per length, on the workstation
python3 scripts/get_chain.py --input-dir gapDesigns/len_58 \
    --output-dir pool/gapDesigns_58 --chain A --min-res 40 --max-res 130
python3 scripts/inventory.py --input-dir pool --output inventory/inventory.tsv
python3 scripts/cluster_scaffolds.py --tsv inventory/inventory.tsv --pool pool \
    --out-dir library --clusters inventory/clusters.tsv
python3 scripts/make_library_manifest.py --clusters inventory/clusters.tsv \
    --inventory inventory/inventory.tsv --batches batches.tsv --library library
```

Add a `batches.tsv` row for the new batch so its provenance is recorded —
these will be `UNVALIDATED` unless they are put through multimer testing.

Then re-run step 1 to regenerate the tensors for the enlarged library.

## Running fold conditioning

```bash
"$BINDER_RFDIFFUSION_ROOT/scripts/run_inference.py" \
    inference.output_prefix=out/boring \
    scaffoldguided.scaffoldguided=True \
    scaffoldguided.scaffold_dir=/path/to/scaffolds \
    scaffoldguided.target_pdb=True \
    scaffoldguided.target_path=input_pdbs/target.pdb \
    scaffoldguided.target_ss=target_folds/target_ss.pt \
    scaffoldguided.target_adj=target_folds/target_adj.pt \
    'ppi.hotspot_res=[A59,A83,A91]' \
    scaffoldguided.mask_loops=False \
    denoiser.noise_scale_ca=0.5 denoiser.noise_scale_frame=0.5
```

`scaffoldguided.scaffold_list=list.txt` restricts a run to a subset — use it
to exclude the `UNVALIDATED` scaffolds (see `library/MANIFEST.tsv`) or to pin
a single length.

To cover lengths between scaffolds without new structures:

```bash
scaffoldguided.mask_loops=True scaffoldguided.sampled_insertion=15 \
    scaffoldguided.sampled_N=5 scaffoldguided.sampled_C=5
```

**RFdiffusion v1 only.** v2 does not document fold conditioning, and v3's
input spec answers "Can I guide on secondary structure?" with "Currently no".
