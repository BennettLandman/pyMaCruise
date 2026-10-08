# Current synthetic base-run core: Java TGDM comparison

## Scope

This comparison closes a fixture-alignment gap in the earlier TGDM report. That
report's Java inputs came from base-run attempt-02. The current corrected,
restartable six-stage base workflow saved its CRUISE-core inputs in attempt-03.
This report uses fresh selected-Java `MaGenericTGDM` outputs computed from the
exact six attempt-03 arrays. It compares the saved attempt-03 Python outputs on
the full native RAS grid. This is a selected synthetic source comparison, not
whole MATLAB/JIST workflow equivalence.

The comparison uses the agreed gates: per-surface mask Dice >= 0.97,
symmetric mean nearest-surface distance <= 0.01 mm, and symmetric Hausdorff
distance <= 1 mm. A source-faithful Python segmentation-refinement operation
was also applied to the Java surfaces and compared with the saved Python
workflow's `SEG_refine.nii.gz`; this is not a direct MATLAB refinement output.

## Results

| Output | Dice | Symmetric mean surface distance (mm) | Hausdorff (mm) |
|---|---:|---:|---:|
| Inner | 1.000000000 | 0.000003533 | 0.002414 |
| Central | 1.000000000 | 0.000001369 | 0.003286 |
| Outer | 0.999989543 | 0.000147829 | 0.200093 |

All three outputs pass the agreed surface gates. Source-faithful refinement
compared with the saved Python labels changed 7 voxels total: labels 1 and 4
were exact, and label 100 had Dice 0.999975151. This passes the 0.97
segmentation threshold. Both output sets use the same `(172, 216, 164)` RAS
grid with 0.9375 mm spacing.

The NIfTI fixtures and their SHA-256 hashes are in
[`tests/data/synthetic/base_core_source_match`](../../tests/data/synthetic/base_core_source_match).
The focused regression checks fixture hashes, RAS geometry, and the recorded
acceptance metrics.

## Reproduction

The following commands assume the corrected synthetic base workflow has
already produced `attempt-03` under the ignored `build/` directory. The Java
harness itself used the saved core inputs; its peak sampled RSS was 1.611 GB.

```sh
cd pymacruise
./.venv/bin/python scripts/prepare_tgdm_core_java_inputs.py \
  build/ace-indexed-heap-2026-10-07/base-run/stages/cruise-core/attempt-03 \
  --target build/tgdm-current-base-run/java_inputs
mkdir -p build/java-reference/tgdm-current-core
javac -d build/java-reference/tgdm-current-core \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/TGDMCoreReferenceHarness.java
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  java -Xms512m -Xmx8g -cp build/java-reference/tgdm-current-core \
  TGDMCoreReferenceHarness build/tgdm-current-base-run/java_inputs \
  build/tgdm-current-base-run/java_reference
./.venv/bin/python scripts/convert_tgdm_core_java_reference.py \
  --input-root build/tgdm-current-base-run/java_inputs \
  --reference-root build/tgdm-current-base-run/java_reference
./.venv/bin/python scripts/compare_tgdm_core_saved_java_reference.py \
  --core-dir build/ace-indexed-heap-2026-10-07/base-run/stages/cruise-core/attempt-03 \
  --reference-root build/tgdm-current-base-run/java_reference \
  --report build/tgdm-current-base-run/saved-output-comparison.json
./.venv/bin/python -m pytest -q tests/test_base_core_current_java_match.py
```

The source-faithful Java-surface refinement comparison used
`refine_segmentation_by_surfaces(..., source_faithful=True)` on the original
synthetic MaSeg input and the fresh Java surfaces, then compared against the
saved Python refinement output. The committed regression records the resulting
label metrics and hashes.

## Real-data status and limits

The user-designated Singularity output directory
`/Users/landmaba/dev/mipav/sample-macruise/output` currently contains zero
files. The input directory contains `T1.nii.gz` and
`orig_target_seg.nii.gz`, but no brain mask. No human-data run or comparison
was attempted. The selected Python interface requires a supplied brain mask
that matches the T1 grid; the target segmentation is not used as a substitute.

This evidence does not close P37: other CRUISE profiles, the direct historical
MATLAB refinement output, plus-workflow parity, remaining source deviations,
and real-data validation remain open.
