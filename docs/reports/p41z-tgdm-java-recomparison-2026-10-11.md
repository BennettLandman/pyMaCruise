# P41z TGDM Java re-comparison after topology-boundary correction

## Result

The earlier corrected-case TGDM comparison that reported two sign changes and
approximately 32.4-mm sampled Hausdorff outliers is superseded. That comparison
used Python level sets written before the source-bounds fix in commit
`6d1d0162`. The fix excludes padded zero-status voxels when MaGenericTGDM's
`SimplePointTest18` fills interior edge and corner neighbors. It matches the
Java bounds checks at crop faces.

On the same held-fixed P41z Java inputs, the current Python nested TGDM
replay now has **zero sign-changing voxels in all six compared fields**. Every
surface has Dice **1.0** and passes the agreed symmetric mean distance
(`≤0.01 mm`) and sampled Hausdorff (`≤1 mm`) limits:

| Field | Sign changes | Dice | Mean distance (mm) | Sampled Hausdorff (mm) | Result |
|---|---:|---:|---:|---:|---|
| Smoothed initial | 0 | 1.0 | 1.33e-7 | 8.35e-5 | Pass |
| Initial inner | 0 | 1.0 | 1.01e-8 | 1.19e-5 | Pass |
| Central | 0 | 1.0 | 5.70e-9 | 7.65e-6 | Pass |
| Inner | 0 | 1.0 | 1.16e-8 | 1.53e-5 | Pass |
| Outer | 0 | 1.0 | 5.39e-7 | 0.00670 | Pass |

`initial_level_set` also has zero sign changes, mean absolute voxel error
`2.92e-10`, and maximum absolute error `4.77e-7`. Java used **0.535 GB** peak
sampled RSS under an 8 GB watchdog. The implementation was tested at tree
`17a3e807`; the source-boundary correction itself is in `6d1d0162`. Full-suite
validation at `cabb52aa` includes that correction: 1,554 passed, 4 skipped, 0
failed.

## Scope and limitations

This is a synthetic, held-fixed selected-Java comparison of the nested TGDM
surface stage. It reuses the hash-checked initial signed-distance input and
reconstructs tissue/ACE/GVF fields on the same 73×91×63 RAS grid. It does not
validate upstream TOADS, LesionTOADS, atlas registration, the live plus loop,
or real human data. The current Python recomputation differs from the older
saved Python TGDM outputs in 8,504 `initial_inner`, 8,493 `inner`, 8,375
`central`, and 8,335 `outer` voxel values. These stale products predate the
correction; the newly regenerated Python outputs pass the held-fixed Java
comparison. The separate P41z workflow checkpoint still fails its ribbon
quality and loop-convergence checks, so a TGDM source match does not make that
case's workflow acceptable.

The previous diagnostic is preserved at
[`p41z-mask-corrected-tgdm-java-comparison-2026-10-11.json`](p41z-mask-corrected-tgdm-java-comparison-2026-10-11.json)
as historical evidence; it must not be treated as the current result. Current
machine-readable metrics are in
[`p41z-tgdm-java-recomparison-2026-10-11.json`](p41z-tgdm-java-recomparison-2026-10-11.json).

## Reproduction

From `pymacruise/`:

```sh
PYTHONPATH=src PYTHONNOUSERSITE=1 MPLCONFIGDIR=build/p41z-current-mpl \
OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  ./scripts/compare_plus_tgdm_saved_java.sh \
  build/p41z-case-normalization/mask-corrected-iteration-0006/plus-surfaces \
  build/p41z-case-normalization/mask-corrected-tgdm-java-current-2026-10-11 \
  docs/reports/p41z-tgdm-java-recomparison-2026-10-11.json \
  build/p41z-case-normalization/mask-corrected-iteration-0006/synthetic-plus-checkpoint-surface-composition.json
```

The command writes Java harness intermediates under ignored `build/` and
returns nonzero if any agreed surface gate fails. The test is a selected-method
source comparison, not a complete PLUS workflow equivalence claim.
