# Fine-grid TGDM smoothing pass localization — 2026-10-10

## Scope

This diagnostic uses the held-fixed 73×91×63 RAS crop from the labeled
81×101×77 synthetic MaCRUISEplus case. It compares pre- and post-reinitialization
fields after each of the six source-selected `SmoothCortexProfile` passes. The
Java harness captures each pass after `GenericTGDM.endNarrowBand`; Python records
the corresponding checkpoints. This is synthetic source-comparison evidence,
not a production algorithm change or real-MRI validation.

The initial Java and Python distance fields agree in sign at all 418,509 voxels
(maximum absolute error `9.54e-7`). Their mean error is `2.16e-8` without
float32 update rounding and `1.26e-10` with rounding. At the end of pass 1,
before fast-marching reinitialization, signs still agree and the maximum field
error is `5.53e-5`. Reinitialization retains the signs but produces a maximum
field difference of `0.0251` level-set units. The first sign disagreement appears in the
curvature update on pass 2. It persists at pass 4 and pass 6; intervening pass
snapshots can temporarily return to equal signs.

| Pass | Sign differences before reinitialization | Sign differences after reinitialization | Maximum field error after reinitialization (level-set units) |
|---:|---:|---:|---:|
| 1 | 0 | 0 | 0.025085 |
| 2 | 1 | 1 | 0.236575 |
| 3 | 0 | 0 | 0.023448 |
| 4 | 2 | 2 | 0.259711 |
| 5 | 0 | 0 | 0.032883 |
| 6 | 1 | 1 | 0.260182 |

Rounding each fast-marching value to float32, as Java does when writing the
image and heap keys, improves the initial field agreement and changes the
number of differing voxels. It does not remove the pass-2 sign difference or
change the final JIST surface gate: Dice `0.999986`, symmetric mean distance
`0.000184 mm`, Hausdorff `0.313 mm`. This experiment therefore does not justify
changing the current production numerical mode. The remaining source discrepancy
spans repeated evolution and distance reinitialization; its root operation is
still open.

## Reproduction

Generate the Java snapshots with the reference harness command in the [TGDM
comparison report](plus-aligned-fine-grid-tgdm-java-2026-10-10.md), then run:

```sh
./.venv/bin/python scripts/compare_plus_tgdm_smoothing_passes.py \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --report docs/reports/plus-aligned-fine-grid-tgdm-smoothing-passes-2026-10-10.json

./.venv/bin/python scripts/compare_plus_tgdm_smoothing_passes.py \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --source-float-rounding \
  --report docs/reports/plus-aligned-fine-grid-tgdm-smoothing-passes-float32-diagnostic-2026-10-10.json
```

See the [default-mode JSON](plus-aligned-fine-grid-tgdm-smoothing-passes-2026-10-10.json)
and [float32-rounding JSON](plus-aligned-fine-grid-tgdm-smoothing-passes-float32-diagnostic-2026-10-10.json)
for all per-pass voxel metrics.
