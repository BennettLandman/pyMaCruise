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
