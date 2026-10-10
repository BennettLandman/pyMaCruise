# Fine-grid TGDM smoothing pass localization — 2026-10-10

## Scope

This diagnostic uses the held-fixed 73×91×63 RAS crop from the labeled
81×101×77 synthetic MaCRUISEplus case. It compares pre- and post-reinitialization
fields after each of the six source-selected `SmoothCortexProfile` passes. The
Java harness captures each pass after `GenericTGDM.endNarrowBand`; Python records
the corresponding checkpoints. This is synthetic selected-source comparison,
not real-MRI validation or complete MaCRUISEplus parity.

The earlier TGDM replay used the ACE signed-distance routine's binary heapq
path. Its initial distance field was close to Java, but small field differences
grew during reinitialization and the first sign difference appeared in the
curvature update on pass 2. Rounding those updates to float32 did not clear the
surface gate. See the preserved [heapq baseline](plus-aligned-fine-grid-tgdm-smoothing-passes-heapq-legacy-2026-10-10.json)
and [float32 heapq diagnostic](plus-aligned-fine-grid-tgdm-smoothing-passes-float32-heapq-diagnostic-2026-10-10.json).

TGDM now reuses the repository's existing source-ordered indexed-heap FMM
implementation at unit speed and applies the Java distance cap. This uses the
same six-neighbor scan order, heap value precision, and image update semantics
as the selected Java `DistanceField`, without adding a second FMM port. The
initial distance fields now agree with Java at all 418,509 voxels by sign, with
mean absolute error `1.26e-10` and maximum `4.77e-7`.

| Pass | Sign differences before reinitialization | Sign differences after reinitialization | Maximum post-reinitialization field error (level-set units) |
|---:|---:|---:|---:|
| 1 | 0 | 0 | 0.00000224 |
| 2 | 0 | 0 | 0.00000584 |
| 3 | 0 | 0 | 0.00003806 |
| 4 | 0 | 0 | 0.00001121 |
| 5 | 0 | 0 | 0.00000983 |
| 6 | 0 | 0 | 0.00001186 |

All six post-reinitialization fields now have identical sign masks. The final
smoothed surface passes the production JIST gate: Dice `1.0`, symmetric mean
distance `1.05e-7 mm`, Hausdorff `1.37e-5 mm`. The full held-fixed TGDM
comparison also passes for initial-inner, central, recomputed inner, and outer
surfaces; see the [Java comparison report](plus-aligned-fine-grid-tgdm-java-2026-10-10.md).

## Reproduction

Generate the Java snapshots with the reference harness command in the [TGDM
comparison report](plus-aligned-fine-grid-tgdm-java-2026-10-10.md), then run:

```sh
./.venv/bin/python scripts/compare_plus_tgdm_smoothing_passes.py \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --report docs/reports/plus-aligned-fine-grid-tgdm-smoothing-passes-2026-10-10.json
```

The [machine-readable result](plus-aligned-fine-grid-tgdm-smoothing-passes-2026-10-10.json)
contains all per-pass voxel metrics and the final surface measurements.
