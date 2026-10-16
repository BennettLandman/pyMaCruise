# P41z LesionTOADS thinning-limit sensitivity — 2026-10-15

## Source-call correction and scope

**Correction, 2026-10-15:** The original source-call interpretation in this
report was wrong. The supplied
`VandyCRUISE/src/edu/jhmi/rad/medic/algorithms/AlgorithmLesionToads.java`
(SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`) forwards
its `firstLimit` field to `LesionToadSegmentation` at line 513. The selected
plugin passes `3.0f` into that field at line 493 of
`MedicAlgorithmLesionToads.java`; the segmentation constructor assigns the
received `firstLimit_` to `limit`. Therefore, the supplied source call path
uses effective thinning limit **3.0**, not zero. This same algorithm source
hash is cited by the earlier iteration-6 Java comparison. The compiled V3
runtime's exact version still needs to be recorded, but the supplied source
does not support the earlier claim that revision 1.3 hard-codes zero.

The limit-zero runs below remain valid **sensitivity experiments**, not source
wrapper probes. The limit-3 selected Java comparison is the branch that
matches the supplied plugin setting. Earlier P41z iterations 1–99 used limit
1000 as an exploratory setting; they are not the supplied wrapper setting.

## Fresh source-aligned synthetic continuation

To remove the exploratory limit-1000 ancestry, a new trajectory was started
from the saved folded synthetic registered recrop. The registered branch
recomputed initial memberships and thinning with `first_thinning_limit=3.0`,
then each continuing loop chunk used `thinning_limit=3.0`. This is source-
aligned LesionTOADS loop evidence on synthetic inputs, not a complete original
Java/MATLAB PLUS workflow.

At completed iteration **27**, the energy ratio is **0.1333700567**; the run
remains nonterminal and has not reached its iteration cap. The latest chunk
(five passes ending at iteration 27) peaked at **1.049 GB** sampled RSS under
a 4 GB watchdog. Its nine manifest outputs all pass their SHA-256 checks; all
eight NIfTIs are finite and RAS with grid `(71, 89, 67)` (channel counts 13,
6, and 17 for the applicable multi-channel products).

A selected Java replay of the fresh limit-3 transition from completed
iteration 6 to 7 matches the Python state exactly: registered shape values
(5,503,849), all **164,674** thinning events and serialized state, all seven
refreshed statistics, membership-energy fields, optimizer transform/rotation/
bounds, **175,225** growth events, post-growth relations, and the saved
checkpoint arrays. The thinning queue began with 97,884 items and drained to
zero. This is selected-method comparison with the registered input held
fixed, not whole-loop or whole-workflow equivalence. The focused comparator
regression tests passed **4/4 in 0.56 s**.

Start or reproduce this sequence from the repository root:

```sh
PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-branch \
  --registered-recrop-manifest build/plus-vandy-atlas-aligned-fine-2iter/registered_recrop/lesion-registered-recrops.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-source-limit3/registered-branch \
  --iterations 99 --first-thinning-limit 3 --first-thinning-volume-ratio 0.5

PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-source-limit3/registered-branch/lesion-registered-branch.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-source-limit3/iterations-1-5 \
  --iterations 99 --thinning-limit 3 --max-passes 5

PYTHONPATH=src .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-source-limit3/iterations-1-5/lesion-registered-loop.json \
  --expected-next-checkpoint build/p41z-source-limit3/iteration-7-from-6/lesion-registered-loop.json \
  --thinning-limit 3 \
  --output-dir build/p41z-source-limit3/iteration-7-java-comparison
```

Continue from each completed loop manifest using the same loop command and
`--thinning-limit 3 --max-passes 5`; this preserves a restartable chain and
keeps each execution bounded. Detailed per-iteration comparison output and
checkpoints are in ignored `build/p41z-source-limit3/`.

To examine sensitivity without relabeling prior work, iteration 89 was
continued with limit zero through iteration 99. This is a **hybrid synthetic
continuation**: iterations before 90 remain the saved limit-1000 state; only
iterations 90–99 use limit zero. It is not the supplied source wrapper
profile, a fresh full Java workflow, or an end-to-end PLUS equivalence run.

## Loop result and selected Java evidence

The limit-0 sensitivity branch reached the selected cap at iteration 99,
with energy
ratio **0.1519980282**, still above the `1e-4` convergence threshold. Ratios
for iterations 91–99 were 0.152720, 0.152066, 0.152070, 0.152005, 0.152007,
0.151989, 0.151995, 0.151988, and 0.151998. The five-pass iteration-91–95
segment peaked at **0.786 GB** sampled RSS; iteration 96 peaked at **0.577 GB**;
iterations 97–99 peaked at **0.736 GB**. A selected Java replay of the
iteration-95→96 transition, explicitly run with limit zero, matched shape
sampling, thinning serialization, seven refreshed statistics, memberships,
registration, growth/checkpoint arrays, and relations exactly. Java and Python
both processed **0 thinning events** and left **97,507 queue items**. The
replay peaked at **0.899 GB** sampled RSS. This remains selected-method
evidence, not a full Java loop execution.

As a sensitivity check, iteration 95 was also continued once with limit **3.0**,
the plugin argument, from the same exact checkpoint. The selected Java replay
at limit 3.0 also matched its Python transition exactly. The Java/Python
comparison records **174,376 growth events**, versus **0** at limit zero. The
two Python iteration-96 segmentations differ at **419 voxels**; skeletons at
424 voxels, availability at 65 voxels, and ordering at 13,427 voxels. Energy
ratios are 0.152386740 (limit 3) and 0.151989207 (limit 0). The small ratio
change therefore does not imply equivalent state evolution. The limit-3 branch
is the one aligned with the supplied plugin-to-segmentation source call.

The same limit-3 branch continued from iteration 96 through 99 and also hit
the cap (ratio **0.131115973**, `stop_reason=iteration_limit`,
`convergence_stop=false`). The 98→99 Java transition matched the selected
methods exactly; terminal-only omitted arrays are listed as unavailable by
the comparator rather than treated as mismatch. This confirms that neither
explicit setting naturally converges on this hybrid synthetic state.

The limit-0 profile changes behavior materially from the exploratory
limit-1000 profile: the 91–99 energy ratios are about **0.152**, compared with
about **0.130** in the earlier continuation. Both synthetic continuations hit their configured cap and neither naturally
converges. The supplied source path passes 3.0; a source-versioned original
run is still needed to establish complete PLUS workflow equivalence and
record the exact compiled runtime.

## Source-profile surface and ribbon audit

The terminal state was composed through the selected Dura, JIST transfer,
MaSeg separation, and source-profile TGDM/surface stages. The composition is
again downstream of the hybrid loop state described above. On the 81×101×77
RAS synthetic fixture, the explicit ribbon audit found:

- GM-between-surfaces support **92.02%**; WM-inside-inner **99.56%**; CSF
  outside outer **99.26%**.
- All three meshes are watertight with no boundary/nonmanifold edges, but
  each has **2 connected components** and Euler characteristic **−16**.
- Level-set mean thickness is **4.765 mm**. Closest-surface median/p05/p95 is
  **6.339/0.231/12.106 mm**.
- The independent diagnostic audit reports `passed=false` for tissue support,
  topology, and production surface-thickness bounds. These anatomy diagnostics
  are separate from source-equivalence gates.

Selected Java TGDM was then run on the exact same Python-generated upstream
arrays and saved initial-level-set seed. All three delivered surfaces pass
the agreed gates (Dice ≥0.97, symmetric mean distance ≤0.01 mm, Hausdorff
≤1 mm): inner and central have zero sign changes; outer has 12, Dice
**0.999949**, mean distance **0.001636 mm**, and Hausdorff **0.354332 mm**.
Java TGDM RSS peaked at **0.658 GB**. This localizes these surface differences
away from the selected TGDM implementation on these fixed inputs; it does not
validate upstream TOADS/LesionTOADS or complete PLUS equivalence.

The limit-3 branch was independently composed and audited as well. GM support
is **93.34%**, WM support **99.44%**, and CSF-outside support **99.25%**. Its
meshes are watertight and connected, but each has Euler characteristic **−16**;
level-set mean thickness is **4.850 mm**, and closest-surface median/p05/p95
is **6.509/0.237/12.341 mm**. The anatomy diagnostic still reports
`passed=false` for support, topology, and thickness. Thus the limit change
alters components and some support values, but neither branch passes the
diagnostic ribbon checks. Selected Java TGDM on the limit-3 branch's same
Python upstream arrays also passes all three delivered-surface gates: the
outer has 9 sign changes, Dice **0.999962**, mean distance **0.001917 mm**, and
Hausdorff **0.465420 mm**; inner and central have zero sign changes. Java TGDM
RSS peaked at **0.662 GB**.

## Reproduction

From the repository root, the limit-zero sensitivity run, selected Java
transition comparison, surface composition, and ribbon audit can be
recreated with:

```sh
PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iteration89-complete/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration89-source-wrapper-limit0 \
  --iterations 99 --thinning-limit 0 --max-passes 1

PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iteration89-source-wrapper-limit0/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iterations90-94-source-wrapper-limit0 \
  --iterations 99 --thinning-limit 0 --max-passes 5

PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iterations90-94-source-wrapper-limit0/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration96-source-wrapper-limit0 \
  --iterations 99 --thinning-limit 0 --max-passes 1

PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iteration96-source-wrapper-limit0/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration99-source-wrapper-limit0 \
  --iterations 99 --thinning-limit 0 --max-passes 3

PYTHONPATH=src .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-case-normalization/iterations90-94-source-wrapper-limit0/lesion-registered-loop.json \
  --expected-next-checkpoint build/p41z-case-normalization/iteration96-source-wrapper-limit0/lesion-registered-loop.json \
  --thinning-limit 0 \
  --output-dir build/p41z-case-normalization/iteration89-source-wrapper-limit0-java-comparison

# From the same iteration-95 checkpoint, test the plugin's configured value.
PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iterations90-94-source-wrapper-limit0/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration96-plugin-limit3 \
  --iterations 99 --thinning-limit 3 --max-passes 1

PYTHONPATH=src .venv/bin/python scripts/compare_p41z_registered_thinning_java.py \
  build/p41z-case-normalization/iterations90-94-source-wrapper-limit0/lesion-registered-loop.json \
  --expected-next-checkpoint build/p41z-case-normalization/iteration96-plugin-limit3/lesion-registered-loop.json \
  --thinning-limit 3 \
  --output-dir build/p41z-case-normalization/iteration96-plugin-limit3-java-comparison

PYTHONPATH=src .venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/p41z-case-normalization/iteration99-source-wrapper-limit0/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --output-dir build/p41z-case-normalization/iteration99-source-wrapper-limit0-surface-audit

PYTHONPATH=src .venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/p41z-case-normalization/iteration99-source-wrapper-limit0-surface-audit/synthetic-plus-terminal-surface-composition.json \
  --report build/p41z-case-normalization/iteration99-source-wrapper-limit0-surface-audit/ribbon-quality.json \
  --maseg-profile synthetic

# Continue the independent limit-3 branch and audit its terminal surfaces.
PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iteration96-plugin-limit3/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iterations97-98-plugin-limit3 \
  --iterations 99 --thinning-limit 3 --max-passes 2

PYTHONPATH=src .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iterations97-98-plugin-limit3/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iteration99-plugin-limit3 \
  --iterations 99 --thinning-limit 3 --max-passes 1

PYTHONPATH=src .venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/p41z-case-normalization/iteration99-plugin-limit3/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --output-dir build/p41z-case-normalization/iteration99-plugin-limit3-surface-audit

PYTHONPATH=src .venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/p41z-case-normalization/iteration99-plugin-limit3-surface-audit/synthetic-plus-terminal-surface-composition.json \
  --report build/p41z-case-normalization/iteration99-plugin-limit3-surface-audit/ribbon-quality.json \
  --maseg-profile synthetic
```

The selected Java TGDM commands and limits are the same as in the [P41z
iteration-99 selected Java report](p41z-iterations90-99-selected-java-2026-10-15.md#same-input-selected-java-tgdm-comparison), with
`iteration99-source-wrapper-limit0-surface-audit` substituted for
`iteration99-terminal-surface-audit`. All generated imaging data remains in
ignored `build/`; no human data or generated images were added to the commit.

## Code-check evidence

The Java comparison utility now accepts a nonnegative `--thinning-limit`
override so it can test a source-passed value without editing or falsifying
the parent checkpoint settings. Its focused regression file passes **4
tests**. This documentation/evidence update did not run the full pytest suite;
the most recent full-suite result remains the separate dated report at the top
of `plan.md`.
