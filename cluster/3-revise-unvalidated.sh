#!/bin/bash
# Data revision for the library's unvalidated scaffolds.
#
# 72 of the 364 cluster representatives have no validation: the RFdiffusion
# gap-designs, the Latent-X designs and the raw CaS decoys. They are in the
# library because each is a distinct fold -- clustering already discarded the
# redundant ones -- so these are exactly the folds worth checking and nothing
# is spent re-checking duplicates.
#
# A multimer test does not apply: these are de novo scaffolds with no target.
# The right question for a fold-conditioning scaffold is whether it is
# *designable* -- does a sequence designed onto this backbone fold back to it?
# That is the pipeline's existing monomer foldability test:
#
#   2-sequenceDesign_solubleMPNN_de_novo.sh   backbone -> sequences
#   3-foldabilityTest_runAF2_monomer.sh       sequences -> AF2 -> RMSD to design
#
# A scaffold whose sequences fold back at low RMSD is real. One that does not
# is a backbone the model cannot realise, and conditioning on it wastes GPU
# time on every future run.
#
# This is a submit wrapper, not a job: it chains the two pipeline steps with
# a dependency so step 3 only starts if step 2 succeeds.
#
# Usage, from the project root on Hive:
#   bash cluster/3-revise-unvalidated.sh [SCAFFOLD_DIR] [OUT_DIR] [NUM_SEQS]

set -euo pipefail

scaffold_dir="${1:-$PWD/revision/scaffolds}"
out_root="${2:-$PWD/revision}"
num_seqs="${3:-8}"

: "${BINDER_REPO_DIR:?set this to your binder-design pipeline checkout}"
: "${SLURM_GPU_PARTITION:=CHANGEME_GPU_PARTITION}"
: "${SLURM_GPU_ACCOUNT:=CHANGEME_GPU_ACCOUNT}"
: "${SLURM_MAIL_USER:=you@example.edu}"

seqdesign="${BINDER_REPO_DIR}/2-SequenceDesign/2-sequenceDesign_solubleMPNN_de_novo.sh"
foldability="${BINDER_REPO_DIR}/3-FoldabilityTest/3-foldabilityTest_runAF2_monomer.sh"

for f in "$seqdesign" "$foldability"; do
    [[ -f "$f" ]] || { echo "ERROR: pipeline step not found: $f" >&2; exit 2; }
done
[[ -d "$scaffold_dir" ]] || { echo "ERROR: no scaffolds at $scaffold_dir" >&2; exit 2; }

n_scaf=$(find "$scaffold_dir" -name '*.pdb' | wc -l)
[[ "$n_scaf" -gt 0 ]] || { echo "ERROR: no PDBs in $scaffold_dir" >&2; exit 2; }

seqs_dir="${out_root}/seqs"
fold_dir="${out_root}/foldability"
mkdir -p "$seqs_dir" "$fold_dir"

export BINDER_REPO_DIR
echo "Revising ${n_scaf} unvalidated scaffolds"
echo "  sequences  -> ${seqs_dir}"
echo "  foldability-> ${fold_dir}"
echo "  repo       -> ${BINDER_REPO_DIR}"

# Step 2: design sequences onto each backbone. Chain A is the whole scaffold
# here -- these are monomers, there is no target chain.
jid_seq=$(sbatch --parsable \
    --partition="$SLURM_GPU_PARTITION" --account="$SLURM_GPU_ACCOUNT" \
    --mail-user="$SLURM_MAIL_USER" \
    "$seqdesign" \
      --input-dir "$scaffold_dir" \
      --output-dir "$seqs_dir" \
      --chains-to-design A \
      --num-seqs "$num_seqs")
echo "  submitted sequence design: ${jid_seq}"

# Step 3: fold each designed sequence with AF2 and compare back to the
# backbone. afterok, so a failed design step does not burn GPU hours here.
jid_fold=$(sbatch --parsable --dependency=afterok:"$jid_seq" \
    --partition="$SLURM_GPU_PARTITION" --account="$SLURM_GPU_ACCOUNT" \
    --mail-user="$SLURM_MAIL_USER" \
    "$foldability" \
      --seqs-dir "$seqs_dir" \
      --output-dir "$fold_dir" \
      --project-name boringRevision)
echo "  submitted foldability test: ${jid_fold}  (afterok:${jid_seq})"

echo
echo "Watch with:  squeue -u \$USER -j ${jid_seq},${jid_fold}"
echo "When done, bring ${fold_dir} back and judge each scaffold on the"
echo "RMSD of its best designed sequence to the original backbone."
