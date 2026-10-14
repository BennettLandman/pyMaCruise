# Human PLUS held-fixed Java comparisons: case 100260

**Evidence category:** real-human held-fixed source comparison after a fresh Python stage replay. The selected Java MaGenericTGDM wrapper sequence matches all agreed surface limits on the same crop-ready Python stage inputs. This is not an independently composed full PLUS Java workflow.

**Correction, 2026-10-14:** The earlier TGDM divergence results in this report were invalid. That comparison reused a saved Python `tissue_sign` field but reconstructed Java's field from current WM/ACE inputs; they differed at **15,707 voxels**. A separate rerun also passed the fresh membership image to Java in precomputed-level-set mode. Both issues are corrected below. Do not use the earlier central/inner/outer mismatch metrics or the old per-voxel trace.

## Workflow case

The comparison uses the axis-only RAS case from the [human PLUS integration report](in-vivo-100260-plus-full-run-2026-10-14.md). The crop is `140×178×134` on a RAS grid. The crop inputs and settings are read from the saved workflow manifest; tissue update, topology, ACE, GVF, and TGDM are freshly recomputed by the current Python implementation. Java receives those same current fields and the source membership image, then follows the selected wrapper sequence. All processes used RSS watchdogs.

The fresh Python replay matches the saved workflow's `initial_level_set` and `initial_inner`. Its `tissue_sign`, delivered `inner`, `central`, and `outer` differ from the saved outputs at 15,707, 1,735,466, 1,717,832, and 2,253,466 voxels respectively. This means the older saved full-workflow output needs a separate provenance/reproducibility investigation; it does not invalidate the matched fresh Python-to-Java comparison below.

## Results

| Selected source stage | Held-fixed result | Agreed criteria |
|---|---|---|
| `CombineMAandTOADS` tissue fusion | GM, WM, and mask each match all **3,339,280 voxels exactly** | Pass |
| MaCRUISE ACE mode 3 | Enhanced GM, skeleton, and thinned skeleton each match all **3,339,280 voxels exactly** | Pass |
| VandyCRUISE GVF (`lambda=0.2`, 200 iterations) | Mean absolute error `2.97e-8`; maximum absolute error `5.22e-7` | Pass current bounds (`1e-5` mean, `5e-5` maximum) |
| TGDM initial signed-distance seed | Mean absolute voxel error `1.40e-10`; max `4.77e-7`; zero sign differences | Pass |
| TGDM initial smoothing | Zero sign changes; Dice `1.0`; mean distance `3.65e-7 mm`; Hausdorff `0.00642 mm` | Pass |
| TGDM first inner surface | Zero sign changes; Dice `1.0`; mean distance `2.18e-7 mm`; Hausdorff `0.00507 mm` | Pass |
| TGDM central surface | Zero sign changes; Dice `1.0`; mean distance `1.51e-7 mm`; Hausdorff `0.00432 mm` | Pass |
| TGDM delivered inner surface | Zero sign changes; Dice `1.0`; mean distance `7.65e-8 mm`; Hausdorff `0.000772 mm` | Pass |
| TGDM outer surface | 14 sign changes; Dice `0.999991`; mean distance `0.000144 mm`; Hausdorff `0.710 mm` | Pass |

All surface results meet the agreed Dice ≥0.97, symmetric mean distance ≤0.01 mm, and Hausdorff ≤1 mm criteria. Tissue fusion, MaCRUISE ACE, and GVF also match the selected Java operations within the stated bounds above. This validates selected stages on a fresh common input chain, not the full PLUS workflow or old saved outputs.

## Smoothing pass localization

The six outer passes of the source `SmoothCortexProfile` wrapper were compared separately, including the curvature update before each reinitialization and the reinitialized field after it. The input fast-marching fields had zero sign differences (mean error `3.67e-11`, maximum `4.77e-7`). None of the six pre- or post-reinitialization pass fields had a sign difference. Post-reinitialization mean absolute voxel error ranged from `4.36e-9` to `1.15e-6`; the final smoothing mesh passed (Dice `1.0`, mean distance `3.65e-7 mm`, Hausdorff `0.00642 mm`). The final pass has a localized maximum level-set value difference of `0.0692 mm`, without changing zero-level-set signs. This is consistent with the fresh matched TGDM result below; it does not prove smoothing values are identical everywhere.

The focused replay exited successfully with peak sampled RSS **1.642 GB** under an 8 GB watchdog. Its machine-readable output is in ignored `build/in-vivo-100260-plus-java-ace-gvf-2026-10-14/smoothing-passes-diagnostic.json`. The report correctly labels this as real-data held-fixed source evidence. Reproduce it with:

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_plus_tgdm_smoothing_passes.py \
  --reference-dir build/in-vivo-100260-plus-java-ace-gvf-2026-10-14 \
  --report build/in-vivo-100260-plus-java-ace-gvf-2026-10-14/smoothing-passes-diagnostic.json
```

The pass summaries below retain the per-outer-iteration error progression (pre-reinitialization mean/max, then post-reinitialization mean/max, all in mm); every field had zero sign differences:

| Pass | Before reinit mean / max | After reinit mean / max |
|---:|---:|---:|
| 1 | `1.42e-8` / `2.19e-4` | `4.36e-9` / `1.85e-6` |
| 2 | `1.99e-8` / `1.07e-4` | `7.88e-9` / `1.75e-5` |
| 3 | `2.43e-8` / `3.33e-4` | `1.23e-8` / `4.51e-4` |
| 4 | `2.81e-8` / `1.62e-4` | `1.53e-8` / `5.73e-5` |
| 5 | `3.18e-8` / `7.36e-4` | `2.76e-8` / `6.39e-4` |
| 6 | `7.75e-8` / `7.42e-2` | `1.15e-6` / `6.92e-2` |

The superseded comparison had two input-chain mismatches: the saved `tissue_sign` was not regenerated from the current held-fixed tissue/ACE arrays, and the Java process was initially invoked as though its input array contained a signed-distance level set. The fresh comparison uses the current Python replay, verifies the sign-field formula exactly on the held-fixed arrays, and invokes Java in membership mode so Java performs the source initialization itself. This removed the apparent central-surface divergence.

## Reproduction

From the `pymacruise` repository root, freshly replay the crop-ready Python chain from the saved case inputs, then run the selected Java TGDM wrapper on those exact current inputs. Run products and metrics stay under ignored `build/`; no human images are committed.

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/in-vivo-100260-plus-no-regrid-lesionatlas-2026-10-14/plus-surface \
  --output-dir build/in-vivo-100260-plus-java-fresh-replay-2026-10-14 \
  --allow-python-stage-differences

REF=build/in-vivo-100260-plus-java-fresh-replay-2026-10-14
mkdir -p "$REF/classes"
javac -d "$REF/classes" \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java
PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx6g -cp "$REF/classes" \
  PlusTGDMReferenceHarness "$REF/arrays.bin" "$REF/java"

PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python \
  scripts/compare_plus_tgdm_java.py \
  --workflow-root build/in-vivo-100260-plus-no-regrid-lesionatlas-2026-10-14/plus-surface \
  --reference-dir "$REF" \
  --report "$REF/tgdm-comparison.json"

scripts/compare_plus_tissue_update_saved_java.sh \
  build/in-vivo-100260-plus-no-regrid-lesionatlas-2026-10-14/plus-surface \
  build/in-vivo-100260-plus-java-tissue-2026-10-14 \
  build/in-vivo-100260-plus-java-tissue-2026-10-14/comparison.json
```

The fresh TGDM Java run peaked at **1.264 GB**, and the Python mesh-distance comparison peaked at **0.599 GB**, each under an 8 GB watchdog. The selected TGDM surface comparison exited successfully with all five surfaces passing. The saved full-workflow stage differences above remain open, so this report does not close P44 or establish whole-workflow equivalence. No full pytest suite was run for this comparison cycle; the source replay and the comparison command are its focused validation.
