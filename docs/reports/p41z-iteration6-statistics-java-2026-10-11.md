# P41z iteration-6 statistics refresh: selected Java comparison

## Result

The post-thinning statistics refresh was compared on the committed 71×89×67
RAS synthetic iteration-6 fixture. The saved Java-matched thinning
segmentation is fixed input, as are the iteration-5 memberships, auxiliary
maps, registered T1/FLAIR images, and iteration-5 centroid state.

The Python refresh and a selected Java excerpt of
`computeCentroids` → `computeVariances` → `computeChannelWeights` match
bit-for-bit for all seven output arrays: class, lesion, and black-hole
centroids; class, lesion, and black-hole standard deviations; and channel
weights. Both report **0 NaNs**. The Java method excerpt is pinned to
`LesionToadSegmentation.java` SHA-256
`545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`; the
selected wrapper is pinned to SHA-256
`556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`.

This comparison exposed a Python behavior gap for `OPTIMIZED_CSF`: selected
Java applies its CSF PowerTable in every modality channel, including FLAIR.
The P41z loop and branch now opt into `optimized_csf_all_modalities`; the
existing default remains unchanged for older source comparisons. With the
production settings—unit tissue powers, lesion powers `(0.4, 1.2)`, no
correction fields, `norm_type=Identity`, and outlier distance `0.1`—the full
volume comparison is exact.

## Reproduction

The default test suite runs two small synthetic comparator cases and skips the
volume-scale comparison unless explicitly requested:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/pytest -q --tb=short \
  tests/test_p41z_next_iteration_statistics_java.py \
  tests/test_lesion_mixed_centroids.py \
  tests/test_lesion_post_thinning_statistics.py
```

Focused result after the correction: **5 passed, 1 skipped**. The opt-in
full-volume check passed as **1 passed, 2 deselected**, under a 16 GB RSS
watchdog with **0.237 GB** peak sampled RSS. It compares raw float32 bits and
counts NaNs explicitly; it does not silently accept NaN equality.

The fixture and runnable command are documented in
[`tests/data/synthetic/p41z_iteration6_statistics_java`](../../tests/data/synthetic/p41z_iteration6_statistics_java/README.md).

## Scope

This is a held-fixed comparison of selected Java methods on one synthetic
iteration checkpoint. It follows the exact-matched thinning output but stops
before the next membership update, energy-ratio decision, registration, or
growth. It is not a complete LesionTOADS iteration, whole-workflow
equivalence, or real-MRI validation.
