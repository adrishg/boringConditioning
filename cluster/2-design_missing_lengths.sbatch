#!/bin/bash
# Generate four-helix bundles at the lengths the library cannot reach.
#
# Lengths 58, 59, 64 and 66 have no 4-helix scaffold anywhere in the lab
# notebook (checked against all 62,806 structures); 60, 61 and 65 exist but
# rest on a single unvalidated structure each.
#
# This uses FOLD CONDITIONING, not unconditional generation. The first
# attempt used monomer_ROG and returned a two-helix hairpin at 58 res: a
# compactness potential controls radius of gyration, not helix count, and at
# these lengths there is barely room for four helices unless the model is
# told to make four.
#
# The conditioning tensors come from scripts/make_synthetic_ss_adj.py. Fold
# conditioning never reads a scaffold's coordinates -- only its
# secondary-structure string and its block-adjacency matrix -- so a real
# bundle's topology can be rescaled to a length that has no PDB. Each target
# length gets 25 distinct topologies copied from real library bundles.
#
# Submit as an array, one task per length:
#   sbatch 2-design_missing_lengths.sbatch [OUT_DIR] [NUM_DESIGNS]
#
# OUT_DIR defaults to <submit dir>/gapDesigns, NUM_DESIGNS to 200.
# Set the partition and account below to ones your user is associated
# with (`sacctmgr show assoc user=$USER`).
#
#SBATCH --partition=CHANGEME_GPU_PARTITION
#SBATCH --account=CHANGEME_GPU_ACCOUNT
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH -t 1-00:00:00
#SBATCH --array=0-6
#SBATCH --job-name=gapScaffolds
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-%x-%A_%a.out
#SBATCH --error=slurm-%x-%A_%a.err

set -euo pipefail

# 58, 59, 64, 66 are missing outright; 60, 61, 65 are single unvalidated folds.
LENGTHS=(58 59 60 61 64 65 66)

if [[ $# -gt 2 ]]; then
    echo "Usage: sbatch $0 [OUT_DIR] [NUM_DESIGNS]" >&2
    exit 2
fi

submit_dir="${SLURM_SUBMIT_DIR:-$PWD}"
out_root="${1:-${submit_dir}/gapDesigns}"
out_root="${out_root%/}"
num_designs="${2:-200}"
[[ "$num_designs" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: NUM_DESIGNS must be a positive integer" >&2; exit 2; }

idx="${SLURM_ARRAY_TASK_ID:-0}"
if [[ "$idx" -ge "${#LENGTHS[@]}" ]]; then
    echo "ERROR: array index $idx exceeds the ${#LENGTHS[@]} configured lengths" >&2
    exit 2
fi
length="${LENGTHS[$idx]}"

: "${BINDER_RFDIFFUSION_ROOT:?set this to your RFdiffusion checkout}"
runner="${BINDER_RFDIFFUSION_ROOT}/scripts/run_inference.py"
[[ -x "$runner" ]] || { echo "ERROR: runner not found: $runner" >&2; exit 2; }

scaffold_src="${SYNTHETIC_SCAFFOLDS:-${submit_dir}/synthetic_scaffolds}"
[[ -d "$scaffold_src" ]] || {
    echo "ERROR: no synthetic scaffolds at $scaffold_src" >&2
    echo "Generate them first -- see cluster/README.md step 0." >&2
    exit 2
}

if [[ -n "${BINDER_REPO_DIR:-}" && -f "${BINDER_REPO_DIR}/bash_tools/load_env.sh" ]]; then
    source "${BINDER_REPO_DIR}/bash_tools/load_env.sh"
    binder_activate_role backbone
else
    # Activate by PREFIX, never by name: micromamba cannot resolve the bare
    # name "SE3nv" without MAMBA_ROOT_PREFIX pointing at this lab's env store;
    # it exits 1 before printing anything. Nothing here redirects stderr -- a
    # silent activation failure is far more expensive than a noisy one.
    env_prefix="${BINDER_ENV_PREFIX:?set this to your SE3nv env prefix}"
    [[ -d "$env_prefix" ]] || { echo "ERROR: env prefix not found: $env_prefix" >&2; exit 2; }
    # set -u must come off around activation: this env ships
    # activate.d/libblas_mkl_activate.sh, which reads MKL_INTERFACE_LAYER
    # while it is unset, and under `set -u` that is fatal.
    set +u
    if command -v micromamba >/dev/null 2>&1; then
        eval "$(micromamba shell hook -s bash)"
        micromamba activate "$env_prefix"
    elif command -v conda >/dev/null 2>&1; then
        source "$(conda info --base)/etc/profile.d/conda.sh"
        conda activate "$env_prefix"
    else
        echo "ERROR: neither micromamba nor conda on PATH" >&2
        exit 2
    fi
    set -u
    [[ -n "${CONDA_PREFIX:-}" ]] || { echo "ERROR: activation produced no CONDA_PREFIX" >&2; exit 2; }
    echo "Environment: ${CONDA_PREFIX}"
fi

out_dir="${out_root}/len_${length}"
mkdir -p "$out_dir"

# RFdiffusion samples a scaffold at random from scaffold_dir, so that
# directory must hold only this length -- otherwise one task emits a mix of
# lengths and the per-length accounting is wrong. Build a per-length view.
len_dir="${out_dir}/scaffolds"
rm -rf "$len_dir"; mkdir -p "$len_dir"
n_scaf=0
for ss in "$scaffold_src"/syn_*_"$(printf '%03d' "$length")"res_*_ss.pt; do
    [[ -e "$ss" ]] || continue
    adj="${ss%_ss.pt}_adj.pt"
    [[ -e "$adj" ]] || { echo "ERROR: $ss has no matching _adj.pt" >&2; exit 3; }
    ln -s "$(readlink -f "$ss")" "$len_dir/"
    ln -s "$(readlink -f "$adj")" "$len_dir/"
    n_scaf=$((n_scaf + 1))
done
[[ $n_scaf -gt 0 ]] || { echo "ERROR: no scaffolds for length ${length} in $scaffold_src" >&2; exit 3; }

prefix="${out_dir}/bundle_${length}res"
echo "Length ${length}: ${num_designs} designs over ${n_scaf} conditioning topologies"

"$runner" \
    inference.output_prefix="$prefix" \
    inference.num_designs="$num_designs" \
    scaffoldguided.scaffoldguided=True \
    scaffoldguided.target_pdb=False \
    scaffoldguided.scaffold_dir="$len_dir" \
    scaffoldguided.mask_loops=False \
    denoiser.noise_scale_ca=0.5 \
    denoiser.noise_scale_frame=0.5

n_out=$(compgen -G "${prefix}*.pdb" | wc -l || true)
if [[ "$n_out" -eq 0 ]]; then
    echo "ERROR: RFdiffusion produced no PDB files for length ${length}" >&2
    exit 3
fi

echo "Length ${length}: ${n_out} backbones in ${out_dir}"
echo
echo "Next, on the workstation:"
echo "  python3 scripts/get_chain.py --input-dir <copied ${out_dir}> \\"
echo "      --output-dir pool/gapDesigns_${length} --chain A"
echo "  python3 scripts/inventory.py --input-dir pool --output inventory/inventory.tsv"
echo "  python3 scripts/cluster_scaffolds.py --tsv inventory/inventory.tsv --pool pool \\"
echo "      --out-dir library --clusters inventory/clusters.tsv"
