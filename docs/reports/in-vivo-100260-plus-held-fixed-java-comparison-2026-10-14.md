# Human PLUS held-fixed Java comparisons: case 100260

**Evidence category:** real-human held-fixed source comparisons. These results cover selected Java stages on the Python workflow's saved crop inputs. They do not establish independent full-workflow Java/MATLAB equivalence; the supplied output tree contains no FLAIR or PLUS products.

## Workflow case

The comparison uses the axis-only RAS case from the [human PLUS integration report](in-vivo-100260-plus-full-run-2026-10-14.md). The crop is `140×178×134` on a RAS grid. LesionTOADS and earlier stages are held fixed. The stage inputs are hash-checked from the saved workflow. The selected Java stages were run with RSS watchdogs; all stage harnesses exited successfully.

## Results

| Selected source stage | Held-fixed result | Agreed criteria |
|---|---|---|
| `CombineMAandTOADS` tissue fusion | GM, WM, and mask each match all **3,339,280 voxels exactly** | Pass |
| MaCRUISE ACE mode 3 | Enhanced GM, skeleton, and thinned skeleton each match all **3,339,280 voxels exactly** | Pass |
| VandyCRUISE GVF (`lambda=0.2`, 200 iterations) | Mean absolute error `2.97e-8`; maximum absolute error `5.22e-7` | Pass current bounds (`1e-5` mean, `5e-5` maximum) |
| TGDM initial signed-distance seed | Exact voxel match | Pass |
| TGDM initial smoothing | Zero sign changes; Dice 1.0; symmetric mean distance `3.65e-7 mm`; Hausdorff `0.00642 mm` | Pass (Dice ≥0.97, mean ≤0.01 mm, Hausdorff ≤1 mm) |
| TGDM first inner surface | Zero sign changes; Dice 1.0; symmetric mean distance `2.18e-7 mm`; Hausdorff `0.00507 mm` | Pass |
| TGDM central surface | Dice `0.998289`; symmetric mean distance `0.01539 mm`; Hausdorff `2.0361 mm`; 2,307 sign changes | **Misses mean and Hausdorff limits** |
| TGDM delivered inner surface | Dice `0.999920`; symmetric mean distance `0.000811 mm`; Hausdorff `2.0053 mm`; 96 sign changes | **Misses Hausdorff limit** |
| TGDM outer surface | Dice `0.997590`; symmetric mean distance `0.02520 mm`; Hausdorff `2.8050 mm`; 3,768 sign changes | **Misses mean and Hausdorff limits** |

To determine how much the small smoothing-value difference contributed, the later Python evolution was replayed with Java's saved `smoothed_initial` substituted. The first inner surface then matched at Dice 1.0 and Hausdorff `0.000012 mm`; delivered inner passed all three limits (Dice `0.999989`, mean `0.000097 mm`, Hausdorff `0.2971 mm`). Central remained just outside the Hausdorff limit (`1.4220 mm`, with mean distance `0.001478 mm`); outer remained just outside it (`1.2546 mm`, with mean distance `0.000642 mm`). This seed-substitution run is diagnostic only; it is not a Python-vs-Java equivalence pass.

The human comparison shows that selected tissue fusion, ACE, and GVF are not the cause of the observed TGDM disagreement. A difference below the surface thresholds after initial smoothing grows in later evolution; the Java-seed substitution reduces that drift substantially, while central and outer maximum-distance outliers remain. This is now the focused source-localization target. No tolerance was changed, and no claim is made that the absolute ribbon-quality flags imply a source mismatch.

## Smoothing pass localization

The six outer passes of the source `SmoothCortexProfile` wrapper were compared separately, including the curvature update before each reinitialization and the reinitialized field after it. The input fast-marching fields had zero sign differences (mean error `3.67e-11`, maximum `4.77e-7`). None of the six pre- or post-reinitialization pass fields had a sign difference. Post-reinitialization mean absolute voxel error ranged from `4.36e-9` to `1.15e-6`; the final smoothing mesh still passed (Dice `1.0`, mean distance `3.65e-7 mm`, Hausdorff `0.00642 mm`). The final pass has a localized maximum level-set value difference of `0.0692 mm`, without changing the zero-level-set signs. This further narrows the human mismatch to the later tissue-guided surfaces; it does not prove smoothing is mathematically identical everywhere.

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

## Tissue-guided evolution localization

The selected Java harness and the Python replay were instrumented at each outer pass, both before and after fast-marching reinitialization. With the Java smoothed seed held fixed, the four-pass initial-inner evolution had no zero-level-set sign differences in any pass. The first sign changes appear in central-surface pass 1, before reinitialization: **2 voxels**. They grow to 129 by the delivered central surface. Later counts include changes inherited from their input surfaces, so they are not independent failure counts.

| Evolution | Outer passes | First pass with sign differences | Final sign differences |
|---|---:|---|---:|
| Initial inner | 4 | None | 0 |
| Central | 7 | Pass 1, before reinitialization: 2 | 129 |
| Central shrink | 1 | Pass 1, before reinitialization: 249 (starts from the already-divergent central field) | 249 |
| Delivered inner | 4 | Pass 1, before reinitialization: 72 (starts from the already-divergent central/shrink path) | 13 |
| Outer | 7 | Pass 1, before reinitialization: 91 (starts from the already-divergent central field) | 72 |

This localizes the first topology-changing discrepancy to central-surface evolution, after the initial-inner result already passes the agreed surface limits. Java tissue fusion and ACE are exact, the sign-force field is assembled from those same held-fixed arrays, and the GVF field is within its numerical bounds. The remaining likely causes are within the central TGDM force/update/topology semantics or the interaction between their floating-point results and the source's topology test. That is the next source-comparison target; these held-fixed results do not establish independent full PLUS equivalence.

The instrumented Java harness completed under an 8 GB watchdog at **1.306 GB** peak sampled RSS. The Python diagnostic completed at **3.153 GB** under an 8 GB watchdog; its nonzero exit is expected because the already-known central and outer surface limits remain unmet. The machine-readable per-pass output is in ignored `build/in-vivo-100260-plus-java-ace-gvf-2026-10-14/java-smoothed-downstream-with-passes.json`. Reproduce the added checkpoints and Python comparison with:

```bash
REF=build/in-vivo-100260-plus-java-ace-gvf-2026-10-14
javac -d "$REF/classes" \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java
PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx6g -cp "$REF/classes" \
  PlusTGDMReferenceHarness "$REF/arrays.bin" "$REF/java" \
  precomputed_initial_level_set
PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/compare_plus_surfaces_from_java_smoothing.py \
  --reference-dir "$REF" \
  --report "$REF/java-smoothed-downstream-with-passes.json"
```

The Python-seed-substitution command returns exit status 1 when the central or outer surfaces fail the agreed distance limits; this is the expected scientific result, while the JSON is still written and the pass comparison completes. The stage-final surface results in this instrumented replay are unchanged from the earlier seed-substitution run: initial inner and delivered inner pass; central and outer Hausdorff distances are `1.4220 mm` and `1.2546 mm`.

## Reproduction

From the `pymacruise` repository root, prepare held-fixed fields from the saved no-regrid workflow. The Java harness comparisons can then be repeated with the scripts below. The run products and JSON metrics are under ignored `build/` and do not add human images to Git.

```bash
PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python \
  scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/in-vivo-100260-plus-no-regrid-lesionatlas-2026-10-14/plus-surface \
  --output-dir build/in-vivo-100260-plus-java-ace-gvf-2026-10-14 \
  --use-saved-tgdm-outputs

PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python \
  scripts/compare_plus_tgdm_java.py \
  --workflow-root build/in-vivo-100260-plus-no-regrid-lesionatlas-2026-10-14/plus-surface \
  --reference-dir build/in-vivo-100260-plus-java-ace-gvf-2026-10-14 \
  --report build/in-vivo-100260-plus-java-ace-gvf-2026-10-14/tgdm-comparison.json

PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python \
  scripts/compare_plus_surfaces_from_java_smoothing.py \
  --reference-dir build/in-vivo-100260-plus-java-ace-gvf-2026-10-14 \
  --report build/in-vivo-100260-plus-java-ace-gvf-2026-10-14/java-smoothed-downstream-diagnostic.json

scripts/compare_plus_tissue_update_saved_java.sh \
  build/in-vivo-100260-plus-no-regrid-lesionatlas-2026-10-14/plus-surface \
  build/in-vivo-100260-plus-java-tissue-2026-10-14 \
  build/in-vivo-100260-plus-java-tissue-2026-10-14/comparison.json
```

The full machine-readable metrics remain in the ignored `build/` directory. The TGDM source run peaked at 1.344 GB Java RSS under an 8 GB watchdog; ACE and GVF peaked at 0.649 GB and 1.331 GB; tissue fusion peaked at 0.254 GB. No full pytest suite was run for these evidence-label and comparator-script changes; the Java/Python comparison commands above are the focused validation for this report.
