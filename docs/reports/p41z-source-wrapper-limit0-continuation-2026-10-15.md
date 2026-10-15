# P41z source-wrapper-limit-0 continuation — 2026-10-15

## Source setting and scope

In the supplied VandyCRUISE source, the plugin passes `firstLimit=3.0` to
`AlgorithmLesionToads` near line 488 of `MedicAlgorithmLesionToads.java`.
However, `AlgorithmLesionToads.java` constructs `LesionToadSegmentation` with
`firstLimit=0` near line 512, rather than forwarding its stored field. The
segmentation class assigns that argument to `limit` near line 290, and
`propagateThinning()` processes its queue only while `val > -limit` (near line
2113). With zero, that pass has no thinning events. The algorithm identifies
itself as CVS revision 1.3 (January 2008); this source does not establish that
this exact call path is the implementation used by the owner's bug-fixed V3
runtime. Therefore, zero is a **source-call-path probe**, not a confirmed V3
production setting. The prior P41z iteration-1–99 continuation used Python
`--thinning-limit 1000` to explore loop continuation behavior. Its exact Java
replays are valid for that explicit method setting, but neither setting
establishes the compiled V3 workflow's behavior.

To examine this call-path probe without relabeling prior work, iteration 89
was continued with limit zero through iteration 99. This is a **hybrid
synthetic continuation**: iterations before 90 remain the saved limit-1000
state; only iterations 90–99 use limit zero. It is not a fresh full Java
workflow, source-equivalent end-to-end PLUS run, or validation of the V3
compiled workflow.

## Loop result and selected Java evidence

The limit-0 continuation reached the selected cap at iteration 99, with energy
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

The limit-0 profile changes behavior materially from the exploratory
limit-1000 profile: the 91–99 energy ratios are about **0.152**, compared with
about **0.130** in the earlier continuation. Both synthetic continuations hit
their configured cap and neither naturally converges. The plugin's 3.0
argument, the algorithm's internal zero, and the compiled V3 behavior need to
be resolved against an executable original reference or version-matched
runtime log.

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

## Reproduction

From the repository root, the zero-limit call-path probe, Java transition
comparison, surface composition, and ribbon audit can be recreated with:

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
