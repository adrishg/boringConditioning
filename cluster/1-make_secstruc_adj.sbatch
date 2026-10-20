#!/bin/bash
# Convert the scaffold library into the _ss.pt / _adj.pt tensors that
# RFdiffusion v1 fold conditioning reads.
#
# Partition, account and mail-user are passed at submit time, matching
# sacredBinderPipeline:
#   sbatch 1-make_secstruc_adj.sbatch [LIBRARY_DIR] [OUT_DIR]
# LIBRARY_DIR defaults to <submit dir>/library, OUT_DIR to its sibling
# scaffolds/.
#
# CPU only -- the helper script is DSSP plus an adjacency matrix, no GPU.
# CPU partition: your GPU account may not be associated with a CPU
# partition. Check with `sacctmgr show assoc user=$USER` and set both
# the partition and account below accordingly.
#SBATCH --partition=CHANGEME_CPU_PARTITION
#SBATCH --account=CHANGEME_CPU_ACCOUNT
#SBATCH -t 0-04:00:00
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --job-name=secstruc_adj
#SBATCH --mail-type=END,FAIL
#SBATCH --output=slurm-%x-%j.out
#SBATCH --error=slurm-%x-%j.err

set -euo pipefail

if [[ $# -gt 2 ]]; then
    echo "Usage: sbatch $0 [LIBRARY_DIR] [OUT_DIR]" >&2
    exit 2
fi

# Defaults so a bare sbatch from the project root does the obvious thing.
library_dir="${1:-${SLURM_SUBMIT_DIR:-$PWD}/library}"
library_dir="${library_dir%/}"
out_dir="${2:-$(dirname "$library_dir")/scaffolds}"

: "${BINDER_RFDIFFUSION_ROOT:?set this to your RFdiffusion checkout}"

helper="${BINDER_RFDIFFUSION_ROOT}/helper_scripts/make_secstruc_adj.py"
[[ -f "$helper" ]] || { echo "ERROR: helper not found: $helper" >&2; exit 2; }
[[ -d "$library_dir" ]] || { echo "ERROR: no library dir: $library_dir" >&2; exit 2; }

# Activate the RFdiffusion env the same way the pipeline does when it is
# available; otherwise fall back to a plain conda activate of SE3nv.
if [[ -n "${BINDER_REPO_DIR:-}" && -f "${BINDER_REPO_DIR}/bash_tools/load_env.sh" ]]; then
    source "${BINDER_REPO_DIR}/bash_tools/load_env.sh"
    binder_activate_role backbone
else
    # Activate by PREFIX, never by name. micromamba cannot resolve the bare
    # name "SE3nv" without MAMBA_ROOT_PREFIX pointing at this lab's env store;
    # it exits 1 before printing anything, which is what killed array 24765197.
    # Nothing here redirects stderr -- a silent activation failure is far more
    # expensive than a noisy one.
    env_prefix="${BINDER_ENV_PREFIX:?set this to your SE3nv env prefix}"
    [[ -d "$env_prefix" ]] || { echo "ERROR: env prefix not found: $env_prefix" >&2; exit 2; }
    # set -u must come off around activation: this env ships
    # activate.d/libblas_mkl_activate.sh, which reads MKL_INTERFACE_LAYER
    # while it is unset. Under `set -u` that is fatal and kills the job.
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

# make_secstruc_adj.py takes ONE --pdb_dir and writes outputs named after each
# input's basename. Our library is nested by topology and length, so build a
# flat directory of symlinks first. Representative names are already unique
# across subdirectories (4h_060res_c01.pdb), so nothing collides -- but check,
# because a silent overwrite here would lose a scaffold.
flat_dir="${out_dir}/_flat_inputs"
mkdir -p "$flat_dir" "$out_dir"
find "$flat_dir" -type l -delete

n_in=0
while IFS= read -r -d '' pdb; do
    base="$(basename "$pdb")"
    if [[ -e "${flat_dir}/${base}" ]]; then
        echo "ERROR: duplicate scaffold name ${base}; refusing to overwrite" >&2
        exit 3
    fi
    ln -s "$(readlink -f "$pdb")" "${flat_dir}/${base}"
    n_in=$((n_in + 1))
done < <(find "$library_dir" -name '*.pdb' -print0)

echo "Converting ${n_in} scaffolds from ${library_dir}"
[[ $n_in -gt 0 ]] || { echo "ERROR: no PDB files found" >&2; exit 2; }

python "$helper" --pdb_dir "$flat_dir" --out_dir "$out_dir"

# One _ss.pt and one _adj.pt per input, or the library is incomplete and a
# scaffoldguided run would silently sample from a smaller set than intended.
n_ss=$(find "$out_dir" -name '*_ss.pt' | wc -l)
n_adj=$(find "$out_dir" -name '*_adj.pt' | wc -l)
echo "Produced ${n_ss} _ss.pt and ${n_adj} _adj.pt for ${n_in} inputs"

if [[ "$n_ss" -ne "$n_in" || "$n_adj" -ne "$n_in" ]]; then
    echo "ERROR: expected ${n_in} of each; some conversions failed" >&2
    exit 3
fi

rm -rf "$flat_dir"
echo "Scaffold tensors written to ${out_dir}"
echo
echo "Use with:"
echo "  scaffoldguided.scaffoldguided=True scaffoldguided.scaffold_dir=${out_dir}"
echo "The target needs its own _ss.pt/_adj.pt from this same helper."
