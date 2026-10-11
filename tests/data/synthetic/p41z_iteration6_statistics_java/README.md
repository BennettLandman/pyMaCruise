# P41z synthetic iteration-6 statistics comparison

This is a self-contained RAS `.nii.gz` fixture for the held-fixed statistics
refresh immediately after iteration-6 thinning. No ignored `build/` files or
external source checkout are needed to run the comparison.

The fixture combines the Java-matched thinning result in
`plus_iteration_6_thinning_result` with iteration-5 memberships and auxiliary
maps from `plus_iteration_6_thinning`, plus the registered T1/FLAIR recrop.
The iteration-5 `refreshed_centroids` initialize the next pass. `case.json`
records hashes for the source checkpoint, thinning state, recropped images,
each committed image, and the selected Java source.

## Files

| File | Contents |
| --- | --- |
| `images_t1_flair.nii.gz` | Registered float32 T1_MPRAGE and FLAIR channels |
| `class_memberships_13ch.nii.gz` | Channels 0–12 from the 17-channel iteration-5 state |
| `lesion_membership.nii.gz` | Channel 13 |
| `black_hole_membership.nii.gz` | Channel 14 |
| `outlier_membership.nii.gz` | Channel 15 |
| `thinned_segmentation.nii.gz` | Java-matched iteration-6 thinning result segmentation |
| `thinned_skeleton.nii.gz` | Java-matched iteration-6 thinning result skeleton |

The registered grid is `71 × 89 × 67`. Production settings are centroid mode
`none`, no correction fields, unit tissue powers, lesion powers `(0.4, 1.2)`,
`norm_type=Identity`, and `outlier_distance=0.1`. The selected VandyCRUISE source applies the optimized-CSF PowerTable independently for every channel (including FLAIR; see `LesionToadSegmentation.java:2477-2486`), so the P41z path explicitly enables `optimized_csf_all_modalities`. The default helper setting remains available for older-source comparisons. The Java harness is a selected
excerpt of `computeCentroids`, `computeVariances`, and
`computeChannelWeights`, pinned to the source SHA in `case.json`. It compares
raw float32 output bits and reports NaN counts; it does not treat NaNs as equal
implicitly.

## Run

The full-volume check is opt-in so the default suite stays short:

```sh
PYMACRUISE_RUN_P41Z_ITER6_STATS=1 PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/pytest -q --tb=short tests/test_p41z_next_iteration_statistics_java.py \
  -k held_fixed
```

The standalone reproducer prints the comparison report:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/compare_p41z_next_iteration_statistics_java.py \
  tests/data/synthetic/p41z_iteration6_statistics_java
```

Evidence from this fixture is a held-fixed source comparison on one synthetic
checkpoint. It is not full-loop, whole-workflow, or human-MRI validation.
