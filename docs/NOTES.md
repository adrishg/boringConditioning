# Notes

Things that cost time when building a scaffold library. Nothing here is
specific to one project.

## Harvesting

**Score tables rarely name files the way they sit on disk.** Pipelines rename
outputs after scoring, so joining on the CSV's `file` column often matches
nothing at all. Rebuild the disk name from the other columns instead, and
check the join resolves before copying anything. A join that silently matches
zero rows looks identical to a campaign with no passing designs.

**The same designs often appear in several directories.** `input_pdbs` and a
sibling `models` directory are frequently byte-identical; a `lowest_dG_poses`
directory is usually a genuinely different pose of the same design. Compare
before treating them as independent sources.

**Check the binder chain, not the shortest chain.** The convention is usually
binder = chain A, target = chain B, but when the target is a short peptide the
"shortest chain" heuristic inverts and silently gives you the target.

## Validate coordinates before trusting anything

Files can be malformed in ways that parse cleanly. If coordinates are large
enough, the three columns run together with no separator and a fixed-column
parser mis-splits them — every value is a number, so nothing errors.

The cheap test is geometry: consecutive CA atoms are ~3.8 Å apart. Anything
much larger inside one chain means missing residues or corrupt coordinates.
One set we found had a quarter of its residues scattered by tens of thousands
of ångströms; DSSP reported zero helices for all of them, which is how it
surfaced.

Also worth checking: multi-MODEL files with differing atom counts, and
filenames whose stated length disagrees with the actual chain length.

## Clustering

**Cluster within each (topology, length) group, never across.** Fold
conditioning takes the design length from the tensors, so a scaffold is only
interchangeable with another of the same length.

**Use structural alignment, not index-wise RMSD.** Two same-length bundles can
place their helix boundaries differently; position-by-position RMSD scores
that as a large difference when the fold is identical.

**TM 0.5 is too loose for this.** It is the usual same-fold line, but it
swallows distinguishable sub-folds, and greedy clustering is seed-order
sensitive — adding one structure collapsed five clusters into one at a length
we were not even touching. 0.6 gave the same answer as 0.7, which is the sign
the threshold has stopped mattering.

**Raw design counts mean nothing.** Large sets collapse hard; one set of 180
designs turned out to be a single fold. Count distinct folds, not files.

## Judging bundles

Helix crossing angle is bimodal: a main mode near 10–40° and a second,
equally real cluster near 50–60° from well-packed but splayed bundles. A 45°
cutoff rejects hundreds of usable structures. Fold conditioning reads only the
secondary structure and which helices touch, so splayed packing is fine —
judge on packing distance and absence of strand, and treat the angle as a
loose sanity bound.

## Filling lengths you do not have

**Unconditional generation with a compactness potential does not work.** It
controls radius of gyration, not helix count. At the short end there is barely
room for four helices, and a two-helix hairpin satisfies the potential
perfectly.

**Fold conditioning does work, and does not need a real scaffold.** It never
reads coordinates — only the per-residue secondary structure and a block
adjacency matrix, and that matrix is purely block-structured (loops excluded,
any helix pair in contact gets its whole submatrix set). Both can be written
directly, so a length with no PDB is not a barrier. Rescaling a real bundle's
measured topology to the target length beats inventing one.

Expect roughly a third of designs to be usable at the shortest lengths, rising
to about two thirds as length increases. Budget accordingly and filter after.

## Environments

Reading structures and writing tensors usually need different environments —
one has the structure library, the other has torch, and rarely both. Split the
work into a measure step and an emit step with a small JSON between them
rather than trying to unify them.

Secondary structure from an in-process DSSP implementation and from the
reference `mkdssp` agree residue-for-residue on backbone-only models, and the
in-process one is orders of magnitude faster. Cross-check once on a sample,
then use the fast path.

## RFdiffusion versions

Fold conditioning exists only in **RFdiffusion v1**. v2 documents motif
scaffolding and never mentions it. v3's input specification answers "Can I
guide on secondary structure?" with "Currently no", offering only a
fewer-loops flag that gives no control over helix count.

## Cluster hygiene

Generate tensors and run analysis in a job, not on a login node, even when the
work looks small. And treat a scheduler's "estimated start" from a dry-run
submission as meaningless — it comes from a single scheduling pass. One job
predicted to wait two days started immediately and finished in under a minute.
