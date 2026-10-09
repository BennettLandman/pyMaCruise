# Dura source comparison on the full-registration synthetic base run — 2026-10-08

## Scope

This runs the selected Java Dura wrapper on the exact N4 T1 and first-pass TOADS white-matter fill, gray-matter, and white-matter-mask inputs from the full-registration base run. The inputs are synthetic, RAS-oriented, and share the `(172,216,164)` grid at `0.9375 mm` spacing. The Java wrapper performs one outer Dura pass with 30 internal TGDM updates. Python is recomputed with the same held-fixed inputs and one-pass profile, then checked against both Java and the saved workflow stage.

## Reproduction

```sh
mkdir -p build/base-workflow-natural-registration-dura-java-2026-10-08/inputs
cp build/base-workflow-2026-10-08-natural-registration/stages/n4/attempt-01/t1-n4.nii.gz \
  build/base-workflow-natural-registration-dura-java-2026-10-08/inputs/t1.nii.gz
cp build/base-workflow-2026-10-08-natural-registration/stages/toads-first/attempt-01/dura-inputs/t1-n4_wmfill.nii.gz \
  build/base-workflow-natural-registration-dura-java-2026-10-08/inputs/white_matter.nii.gz
cp build/base-workflow-2026-10-08-natural-registration/stages/toads-first/attempt-01/dura-inputs/t1-n4_gm.nii.gz \
  build/base-workflow-natural-registration-dura-java-2026-10-08/inputs/gray_matter.nii.gz
cp build/base-workflow-2026-10-08-natural-registration/stages/toads-first/attempt-01/dura-inputs/t1-n4_wmmask.nii.gz \
  build/base-workflow-natural-registration-dura-java-2026-10-08/inputs/white_matter_mask.nii.gz
./.venv/bin/python scripts/prepare_dura_base_java_input.py \
  build/base-workflow-natural-registration-dura-java-2026-10-08/inputs \
  build/base-workflow-natural-registration-dura-java-2026-10-08/dura-base-input.bin
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 20 -- \
  env DURA_WRAPPER_JAVA_HEAP=4g ./scripts/run_java_dura_base_reference.sh \
  build/base-workflow-natural-registration-dura-java-2026-10-08/dura-base-input.bin \
  build/base-workflow-natural-registration-dura-java-2026-10-08/java-output
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 20 -- \
  ./.venv/bin/python scripts/compare_dura_saved_java.py \
  --input-dir build/base-workflow-natural-registration-dura-java-2026-10-08/inputs \
  --java-output-dir build/base-workflow-natural-registration-dura-java-2026-10-08/java-output \
  --base-stage build/base-workflow-2026-10-08-natural-registration/stages/dura/attempt-01 \
  --report build/base-workflow-natural-registration-dura-java-2026-10-08/comparison.json
```

The Java harness exited 0 under its 8 GB watchdog. Python's direct recomputation/comparison exited 0 at **2.106 GB** peak sampled RSS and reproduced all eight Dura products in the saved stage byte-for-byte.

## Results

The WM topology, GM topology, and WM+GM union are exact against Java. Initial distance fields have no sign disagreements; the inner maximum absolute error is `9.54e-07`, and the outer mean absolute error is `1.61e-05`. The one-pass evolving level set has one sign disagreement in 6,092,928 voxels; zero-threshold Dice is `0.99999943`.

The final brain mask differs at **5 voxels**. Foreground Dice is **0.99999768**, symmetric mean nearest-surface distance is **0.00004772 mm**, and Hausdorff distance is **0.66291 mm**. These are within the agreed Dice ≥0.97, mean distance ≤0.01 mm, and Hausdorff ≤1 mm gates. Changes in stripped T1 and corrected WM/GM are limited to those same mask-boundary voxels.

Input hashes, all 12 Java/Python output metrics, workflow-reproduction flags, and mesh counts are saved in the linked [full-registration machine-readable evidence](base-workflow-natural-registration-2026-10-08.json). This is a selected Java Dura wrapper comparison on held-fixed synthetic inputs, not whole-workflow equivalence or real-MRI validation.
