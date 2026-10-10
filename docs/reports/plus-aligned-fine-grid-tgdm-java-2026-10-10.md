# Fine-grid synthetic plus TGDM Java comparison — 2026-10-10

## Scope

This is a held-fixed comparison of the selected `MaGenericTGDM` wrapper with
the saved Python TGDM outputs from the iteration-6 checkpoint of the labeled
81×101×77 (~2 mm) synthetic plus case. The comparison is on the workflow's
73×91×63 RAS crop. Hash-checked crop inputs were used to reconstruct the
tissue, ACE, and GVF fields. Java and Python start from the same held-fixed
initial-membership array and independently construct the signed-distance
seed. The seed fields have identical signs and differ by at most `4.77e-7`.
TOADS/LesionTOADS and all other upstream stages remain held fixed. The Java
process peaked at **0.55 GB sampled RSS** under an 8 GB
watchdog. This is synthetic component source-comparison evidence, not
full-workflow equivalence or real-MRI validation.

## Results

The production TGDM replay now uses the source-ordered indexed-heap fast-marching
routine for signed-distance reinitialization. The Java and Python outputs have
zero sign differences at every compared stage. Production distances below use
the source-matched JIST connectivity-18 extractor; marching-cubes values are
retained in the JSON as diagnostics only. All five surfaces pass Dice ≥0.97,
mean distance ≤0.01 mm, and Hausdorff ≤1 mm:

| Surface | Sign-different voxels | Dice at zero | Mean surface distance | Hausdorff | Gate |
|---|---:|---:|---:|---:|---|
| Smoothed initial | 0 | 1.000000 | 0.000000105 mm | 0.0000137 mm | Pass |
| Initial inner | 0 | 1.000000 | 0.000000010 mm | 0.0000119 mm | Pass |
| Central | 0 | 1.000000 | 0.000000008 mm | 0.0000076 mm | Pass |
| Inner | 0 | 1.000000 | 0.000000009 mm | 0.0000133 mm | Pass |
| Outer | 0 | 1.000000 | 0.000000010 mm | 0.0000221 mm | Pass |

The level-set arrays are not bit-identical, but the largest voxel error among
these final surfaces is `1.97e-5` level-set units, and every zero-level mask
matches Java exactly. This clears the selected TGDM component gates on the
held-fixed crop. It does not establish full MaCRUISEplus parity, final ribbon
quality, or real-MRI performance; the labeled full-plus checkpoint still needs
to be regenerated with this corrected fast-marching backend and re-audited.

As a diagnostic, Python was rerun from the exact Java smoothed level set through
the initial inner, central, recomputed inner, and outer stages. All four surfaces pass the Java thresholds with zero sign differences and Dice 1.0 using production JIST surfaces:

| Surface | Mean surface distance | Hausdorff |
|---|---:|---:|
| Initial inner | 0.00000062 mm | 0.001588 mm |
| Central | 0.00000479 mm | 0.003462 mm |
| Inner | 0.00000058 mm | 0.000453 mm |
| Outer | 0.00000437 mm | 0.002890 mm |

The diagnostic uses source-matched JIST mesh extraction for its primary gate;
marching-cubes values remain secondary diagnostics. It deliberately substitutes
a Java intermediate and is not an end-to-end Python-vs-Java result. The normal
held-fixed replay already passes the same four downstream stages, so this
control adds consistency evidence only. The focused TGDM and Dura reference
regressions passed **27 tests** after the source-ordered fast-marching change.
The downstream diagnostic is reproducible with the committed
[comparison script](../../scripts/compare_plus_surfaces_from_java_smoothing.py)
and its [machine-readable report](plus-aligned-fine-grid-java-smoothing-downstream-diagnostic-2026-10-10.json).

Per-pass pre/post-reinitialization snapshots showed that the prior heapq-based
FMM path first diverged in sign during the pass-2 curvature update. Reusing the
existing source-ordered indexed-heap FMM implementation for unit-speed TGDM
reinitialization removes that mismatch: all six passes now have zero sign
differences, and the final smoothed surface has Dice 1.0, mean distance
`1.05e-7 mm`, and Hausdorff `1.37e-5 mm`. See the
[pass-localization report](plus-aligned-fine-grid-tgdm-smoothing-pass-localization-2026-10-10.md)
and current [machine-readable pass report](plus-aligned-fine-grid-tgdm-smoothing-passes-2026-10-10.json).

## Reproduction

The held-fixed inputs and Java products are generated under ignored `build/`
directories. Recreate the inputs and compile/run the selected Java harness:

```sh
mkdir -p build/plus-vandy-atlas-aligned-fine/java-reference/classes

./.venv/bin/python scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/plus-surfaces \
  --context-manifest build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/synthetic-plus-checkpoint-surface-composition.json \
  --output-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --allow-python-stage-differences

javac -d build/plus-vandy-atlas-aligned-fine/java-reference/classes \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java

./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 8 \
  --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx6g \
  -cp build/plus-vandy-atlas-aligned-fine/java-reference/classes \
  PlusTGDMReferenceHarness \
  build/plus-vandy-atlas-aligned-fine/java-inputs/arrays.bin \
  build/plus-vandy-atlas-aligned-fine/java-inputs/java \
  > build/plus-vandy-atlas-aligned-fine/java-inputs/java-watchdog.log 2>&1

./.venv/bin/python scripts/compare_plus_tgdm_java.py \
  --workflow-root build/plus-vandy-atlas-aligned-fine/surface-checkpoint-0006/plus-surfaces \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --report docs/reports/plus-aligned-fine-grid-tgdm-java-2026-10-10.json

./.venv/bin/python scripts/compare_plus_surfaces_from_java_smoothing.py \
  --reference-dir build/plus-vandy-atlas-aligned-fine/java-inputs \
  --report docs/reports/plus-aligned-fine-grid-java-smoothing-downstream-diagnostic-2026-10-10.json
```

The complete measurements and hashes are in the adjacent machine-readable
[JSON report](plus-aligned-fine-grid-tgdm-java-2026-10-10.json).
