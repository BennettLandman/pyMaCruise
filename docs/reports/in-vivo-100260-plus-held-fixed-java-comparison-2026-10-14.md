# Human PLUS held-fixed Java comparisons: case 100260

**Evidence category:** real-human held-fixed source comparison after a fresh Python stage replay. The selected Java MaGenericTGDM wrapper sequence matches all agreed surface limits on the same crop-ready Python stage inputs. This is not an independently composed full PLUS Java workflow.

**Correction, 2026-10-14:** Earlier TGDM divergence results in this report were invalid for three reasons. One comparison reused a saved Python `tissue_sign` field but reconstructed Java's field from current WM/ACE inputs; they differed at 15,707 voxels. A separate rerun passed the membership image to Java in precomputed-level-set mode. Finally, the Python replay helper passed TOADS WM fill as the boundary-cleanup input, while the full PLUS workflow passes TOADS WM mask. With the source-selected input order corrected, a fresh full Python run reproduces all six saved TGDM outputs bit-for-bit, and the selected Java comparison below passes all surface gates. Do not use the superseded mismatch metrics or point trace.

## Workflow case

The comparison uses the axis-only RAS case from the [human PLUS integration report](in-vivo-100260-plus-full-run-2026-10-14.md). The crop is `140×178×134` on a RAS grid. The crop inputs and settings are read from the saved workflow manifest; tissue update, topology, ACE, GVF, and TGDM are freshly recomputed by the current Python implementation. Java receives those same current fields and the source membership image, then follows the selected wrapper sequence. All processes used RSS watchdogs.

The new raw-input no-regrid PLUS run completed all nine stages and naturally converged at LesionTOADS iteration 13. Its run signature is unchanged from the previous no-regrid run (`5061e1a6…`). Seventy manifest-referenced NIfTI/VTK products across the workflow's stage manifests have identical hashes between runs. The six crop-grid TGDM outputs (`tissue_sign`, `initial_level_set`, `initial_inner`, `inner`, `central`, and `outer`) match exactly. The corrected crop-stage replay also reproduced all six Python TGDM products exactly (zero differing voxels per stage), confirming that the earlier mismatch came from the replay helper's boundary input order, not from full-workflow nondeterminism.

## Results

| Selected source stage | Held-fixed result | Agreed criteria |
|---|---|---|
| `CombineMAandTOADS` tissue fusion | GM, WM, and mask each match all **3,339,280 voxels exactly** | Pass |
| MaCRUISE ACE mode 3 | Enhanced GM, skeleton, and thinned skeleton each match all **3,339,280 voxels exactly** | Pass |
| VandyCRUISE GVF (`lambda=0.2`, 200 iterations) | Mean absolute error `2.97e-8`; maximum absolute error `5.22e-7` | Pass current bounds (`1e-5` mean, `5e-5` maximum) |
| TGDM initial signed-distance seed | 983 differing voxels from float rounding; mean absolute voxel error `1.40e-10`; max `4.77e-7`; zero sign differences | Pass |
| TGDM initial smoothing | Zero sign changes; Dice `1.0`; mean distance `3.65e-7 mm`; Hausdorff `0.00642 mm` | Pass |
| TGDM first inner surface | Zero sign changes; Dice `1.0`; mean distance `2.18e-7 mm`; Hausdorff `0.00507 mm` | Pass |
| TGDM central surface | Zero sign changes; Dice `1.0`; mean distance `1.51e-7 mm`; Hausdorff `0.00432 mm` | Pass |
| TGDM delivered inner surface | Zero sign changes; Dice `1.0`; mean distance `7.68e-8 mm`; Hausdorff `0.000766 mm` | Pass |
| TGDM outer surface | 10 sign changes; Dice `0.999994`; mean distance `0.000110 mm`; Hausdorff `0.2255 mm` | Pass |

All surface results meet the agreed Dice ≥0.97, symmetric mean distance ≤0.01 mm, and Hausdorff ≤1 mm criteria. Tissue fusion, MaCRUISE ACE, and GVF also match the selected Java operations within the stated bounds above. The full Python workflow is reproducible against its previous saved run; the selected Java comparisons remain held-fixed stage evidence, not a complete independent PLUS Java/MATLAB run.

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

The superseded comparisons mixed saved and freshly derived fields, once selected the wrong Java input mode, and once passed WM fill instead of WM mask to the boundary stage. The corrected replay now follows the full PLUS workflow's eight-input order and enforces exact agreement with its saved TGDM products before the Java comparison. Java receives the same membership image and performs the source distance initialization itself.

## Reproduction

From the `pymacruise` repository root, freshly run the axis-only RAS full PLUS workflow on the case inputs, then replay its crop-ready Python chain with exact saved-output checks, and run the selected Java TGDM wrapper on those same inputs. Products and metrics stay under ignored `build/`; no human images are committed.

```bash
MPLCONFIGDIR=/private/tmp/mpl-pymacruise-fresh-plus \
PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --report-interval-seconds 30 -- \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-mp2rage_T1w.nii.gz \
  --flair /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-highres_FLAIR.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-mp2rage_T1w_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/100260/sub-010260/ses-01/SLANT-TICVv1.2/post/FinalResult/sub-010260_ses-01_T1w_seg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14 \
  --no-regrid-inputs

PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --report-interval-seconds 30 -- \
  ./.venv/bin/python scripts/prepare_plus_tgdm_java_inputs.py \
  --workflow-root build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/plus-surface \
  --output-dir build/in-vivo-100260-plus-java-corrected-replay-2026-10-14

REF=build/in-vivo-100260-plus-java-corrected-replay-2026-10-14
mkdir -p "$REF/classes"
javac -d "$REF/classes" \
  -sourcepath tests/java_reference/stubs:tests/java_reference/upstream \
  tests/java_reference/PlusTGDMReferenceHarness.java
set -o pipefail
PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  /usr/bin/java -Xms256m -Xmx6g -cp "$REF/classes" \
  PlusTGDMReferenceHarness "$REF/arrays.bin" "$REF/java" 2>&1 | tee "$REF/java-watchdog.log"

PYTHONNOUSERSITE=1 PYTHONPATH=src ./.venv/bin/python \
  scripts/compare_plus_tgdm_java.py \
  --workflow-root build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14/plus-surface \
  --reference-dir "$REF" \
  --report "$REF/tgdm-comparison.json"

scripts/compare_plus_tissue_update_saved_java.sh \
  build/in-vivo-100260-plus-no-regrid-lesionatlas-2026-10-14/plus-surface \
  build/in-vivo-100260-plus-java-tissue-2026-10-14 \
  build/in-vivo-100260-plus-java-tissue-2026-10-14/comparison.json

MPLCONFIGDIR=/private/tmp/mpl-pymacruise-fresh-plus \
PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 PYTHONPATH=src \
  ./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --report-interval-seconds 30 -- \
  ./.venv/bin/python -m pymacruise.cli plus-full-run \
  --t1 /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-mp2rage_T1w.nii.gz \
  --flair /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-highres_FLAIR.nii.gz \
  --brain-mask /Users/landmaba/dev/mipav/sample-macruise/100260/input/ses-01/anat/sub-010260_ses-01_acq-mp2rage_T1w_bet_mask.nii.gz \
  --maseg /Users/landmaba/dev/mipav/sample-macruise/100260/sub-010260/ses-01/SLANT-TICVv1.2/post/FinalResult/sub-010260_ses-01_T1w_seg.nii.gz \
  --atlas-descriptor data/atlases/lesionToads-atlas-2012/cruise-atlas-12obj-lesiontoads2012.txt \
  --output-dir build/in-vivo-100260-plus-fresh-current-no-regrid-2026-10-14 \
  --no-regrid-inputs --resume
```

The corrected crop replay peaked at **2.599 GB** under an 8 GB watchdog; all six saved TGDM outputs matched exactly. The selected Java TGDM run peaked at **1.274 GB**, and the mesh-distance comparison peaked at **0.597 GB**, each under an 8 GB watchdog. All five surface outputs pass the agreed limits. This remains selected-stage source evidence, not an independent full PLUS Java/MATLAB comparison; P44 stays open for that workflow-level source evidence. The focused input-order regression passed **1 test in 0.54 s**; no full pytest suite was run for this comparison cycle.
