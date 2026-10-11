# P41z corrected-type synthetic rerun — 2026-10-11

## Result

I resumed the committed iteration-5 registered checkpoint with the uppercase
standard-13 object-type normalization fix from commit `0929c81f`. The one-pass
iteration-6 run completed and wrote a restartable registered-loop manifest.
It remains nonterminal, with energy ratio **0.13159732520580292** and objective
**3267.361083984375**. It completed in **181.59 s** with **0.782 GB** peak
sampled RSS under an **8 GB** watchdog.

The run was compared with the earlier saved iteration-6 output created before
the type-normalization fix. The corrected run changes **1,306** segmentation
and skeleton voxels, **105,653** membership values (maximum absolute
difference `2.3841858e-7`), **3,462** ordering voxels, and **3,818** relation
map voxels. Availability is unchanged. The growth trace has 254,394 events,
three fewer than the earlier 254,397. The objective, previous objective,
registration transform, and energy ratio remain bitwise unchanged. This
confirms the production fix affects discrete growth state, while the current
energy convergence blocker lies elsewhere or is not sensitive to this change.

The corrected checkpoint was composed through Dura, native transfer, MaSeg,
and the surface workflow. The explicit synthetic ribbon audit still fails:

- GM between the inner and outer surfaces: **0%** (97% required).
- WM inside the inner surface: **59.7655%** (97% required).
- CSF outside the outer surface: **99.9827%**.
- Each mesh: **2 components**, watertight, Euler characteristic **−16**.
- Closest-point triangle thickness: median **0.3725 mm**, p05 **0.1644 mm**,
  p95 **0.6415 mm**; nominal ribbon is 5 mm.
- Level-set subtraction thickness: mean **0.8106 mm** across 291 crossing
  voxels.

Thus the upstream growth correction changes the surface slightly (previous
Euler characteristic −20; now −16) but does not resolve the tissue-support,
topology, or thickness failures. The machine-readable metrics are in
[`p41z-case-normalization-ribbon-quality-2026-10-11.json`](p41z-case-normalization-ribbon-quality-2026-10-11.json).

## Reproduction

From the repository root, resume the saved iteration-5 checkpoint:

```sh
env PYTHONPATH=src PYTHONNOUSERSITE=1 \
  MPLCONFIGDIR=build/p41z-case-normalization/mpl-cache \
  OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
    --limit-gb 8 --interval-seconds 1 --report-interval-seconds 30 -- \
    .venv/bin/python -m pymacruise.cli plus-lesion-registered-loop \
      --registered-checkpoint-manifest \
        build/plus-vandy-atlas-aligned-fine/loop-chunks/iteration-0005/lesion-registered-loop.json \
      --atlas-descriptor \
        tests/data/synthetic/plus_vandy_atlas_aligned_fine/cruise-atlas-12obj-lesiontoads2012.txt \
      --output-dir build/p41z-case-normalization/iteration-0006 \
      --iterations 99 --thinning-limit 1000 --max-passes 1
```

Compose and audit the nonterminal checkpoint with:

```sh
env PYTHONPATH=src .venv/bin/python scripts/compose_terminal_synthetic_plus.py \
  --loop-manifest build/p41z-case-normalization/iteration-0006/lesion-registered-loop.json \
  --t1 tests/data/synthetic/plus_vandy_atlas_aligned_fine/t1.nii.gz \
  --brain-mask tests/data/synthetic/plus_vandy_atlas_aligned_fine/brain-mask.nii.gz \
  --maseg tests/data/synthetic/plus_vandy_atlas_aligned_fine/maseg.nii.gz \
  --output-dir build/p41z-case-normalization/iteration-0006-surfaces \
  --allow-checkpoint

env PYTHONPATH=src .venv/bin/python scripts/audit_plus_ribbon_quality.py \
  build/p41z-case-normalization/iteration-0006-surfaces/synthetic-plus-checkpoint-surface-composition.json \
  --allow-checkpoint \
  --report build/p41z-case-normalization/iteration-0006-surfaces/ribbon-quality.json
```

The focused lesion-growth/thinning regression group passed **17 tests, 0
failed, 0 skipped in 0.86 s**. The selected Java handoff comparison passed
**1 test in 200.71 s** under a 4 GB watchdog; it covers iteration-6 growth,
post-growth relations, 5,503,849 registered shape values, and the full
iteration-7 thinning serialization. These are synthetic selected-method
comparisons and workflow checks, not complete Java-loop equivalence or
real-MRI validation. No full suite was run after commit `0929c81f`.

## Next checkpoint

Obtain a Java-authenticated previous objective and compare the next complete
registered iteration state before changing any ribbon or convergence logic.
Then trace where MaSeg cortical GM leaves the ribbon using the matched
iteration outputs. The 99-iteration workflow still has not naturally
converged, and the corrected iteration-6 ratio remains above `1e-4`.
