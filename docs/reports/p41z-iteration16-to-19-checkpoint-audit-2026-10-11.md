# P41z checkpoint continuation and surface audit, iterations 16–19 — 2026-10-11

## Runnable workflow and results

The source-selected synthetic LesionTOADS settings were kept unchanged:
`n_iter_max=99`, four iterations per restartable segment, thinning limit 1000,
source zero-variance NaN semantics, and convergence threshold `1e-4`. The
continuation started from the committed iteration-15 workflow artifact in the
ignored build tree and completed iterations 16–19 in **719.66 seconds**. Peak
sampled RSS was **1.091 GB** under a **4 GB** watchdog. All nine registered
outputs listed in its manifest passed SHA-256 verification.

| Iteration | Energy ratio | Stop |
| ---: | ---: | --- |
| 16 | 0.13298921 | No |
| 17 | 0.13240865 | No |
| 18 | 0.13197617 | No |
| 19 | 0.13316543 | No |

The output remains a nonterminal checkpoint. Across iterations 8–19 the ratio
oscillates between **0.131976 and 0.133960** without a sustained downward
trend.

The same source-selected Dura, physical-affine transfer, MaSeg separation, and
plus surface tail were composed from iteration 19. The explicit ribbon-quality
audit reports:

| Measure | Iteration 15 | Iteration 19 | Acceptance |
| --- | ---: | ---: | ---: |
| WM label support inside inner | 99.435% | 99.440% | ≥97% |
| GM label support between surfaces | 92.890% | 92.906% | ≥97% |
| CSF label support outside outer | 99.258% | 99.255% | ≥97% |
| Mesh Euler characteristic | −18 | −16 | 2 |
| Surface thickness median | 6.496 mm | 6.505 mm | 4–6 mm |
| Surface thickness p05 / p95 | 0.236 / 12.298 mm | 0.236 / 12.321 mm | 2–8 mm |
| Level-set mean thickness | 4.842 mm | 4.845 mm | 4–6 mm |
| Nested voxel regions | Pass | Pass | Pass |

At iteration 19 all surfaces are closed and connected, but their Euler
characteristic is −16 rather than 2. The slight GM-support increase and Euler
change do not bring the surfaces within acceptance. The level-set mean and
nesting pass; the production surface thickness and tissue support do not.

## Selected Java TGDM comparison

The selected Java `MaGenericTGDM` sequence was replayed against this
iteration-19 saved Python state. It passed all five surface gates: four
surfaces have zero sign changes; the outer surface has six and Dice
**0.999975**, symmetric mean distance **0.001899 mm**, and Hausdorff
**0.421216 mm**. The Java watchdog peak was **0.622 GB** under an 8 GB limit.
The six saved products used to prepare the Java inputs replayed exactly in
Python.

This is a selected-stage synthetic comparison. It reuses the saved Python
signed-distance seed, reconstructs tissue/ACE/GVF inputs on the same crop, and
holds TOADS/LesionTOADS fixed. It is not full registered-loop or whole-workflow
Java equivalence. The TGDM surface match does not explain or remove the
iteration-19 anatomy-quality failures.

## Reproduction

Continue the registered loop:

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 4 \
  --interval-seconds 1 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
  --registered-checkpoint-manifest build/p41z-case-normalization/iterations-0012-to-0015/lesion-registered-loop.json \
  --atlas-descriptor tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/p41z-case-normalization/iterations-0016-to-0019 \
  --iterations 99 --thinning-limit 1000 --max-passes 4
```

Compose the iteration-19 surface tail and perform the explicit quality audit:

```sh
MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/p41z-case-normalization/iterations-0016-to-0019/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --output-dir build/p41z-case-normalization/iterations-0016-to-0019-surfaces \
  --allow-checkpoint

MPLCONFIGDIR="$PWD/build/matplotlib-config" PYTHONPATH=src \
  .venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/synthetic-plus-checkpoint-surface-composition.json \
  --allow-checkpoint \
  --report build/p41z-case-normalization/iterations-0016-to-0019-surfaces/ribbon-quality.json
```

The iteration-19 Java comparison is reproducible with:

```sh
scripts/compare_plus_tgdm_saved_java.sh \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/plus-surfaces \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/tgdm-java \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/tgdm-java-comparison.json \
  build/p41z-case-normalization/iterations-0016-to-0019-surfaces/plus-surfaces/run.json
```

## Validation state and next milestone

The latest full `pytest` suite was **1,559 passed, 4 skipped, 0 failed** on
implementation commit `3bf4aeb2` in 827.07 seconds under a 16 GB watchdog
(11.211 GB peak); this iteration continuation, surface composition, and Java
replay are subsequent integration/source-comparison runs, not a new full test
suite. No Python implementation changed for these runs.

The current blocker is no longer memory or a demonstrated TGDM surface
mismatch. It is the nonterminal registered loop and failed GM support, genus,
and production thickness checks. The next integration milestone is to compare
a later registered-loop cycle against selected Java methods at the same saved
state, identify any state divergence before spending more time extending the
loop, and retain the failed ribbon audit as the acceptance check. This is
synthetic evidence, not real-human-MRI acceptance or a no-go conclusion.
