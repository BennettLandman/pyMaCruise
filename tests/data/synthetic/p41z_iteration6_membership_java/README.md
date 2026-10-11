# P41z iteration-6 held-fixed membership comparison

This synthetic case compares one ordered 13-class membership update after the
iteration-6 thinning result. The fixed iteration-5 retained relation state is
recovered from the hash-verified iteration-5 manifest and summary; relation
flags are looked up from the retained relation-membership table. They are not
reconstructed from the thinned segmentation.

Every image needed by the comparison is materialized here as RAS `.nii.gz` and
is listed with its hash, shape, and datatype in `case.json`. `case.json` also
records the iteration-5 and iteration-6 checkpoint hashes, selected Java
source and wrapper hashes, and the exact profile. The image channels are T1
MPRAGE and FLAIR. The comparison uses Identity channel weights, outlier energy
0.1, lesion weighting on FLAIR with exponent 1.2, unit shape and smoothing
priors, GMflag enabled, and no correction fields. Statistics are refreshed
first using the previously validated P41z iteration-6 statistics profile.

The source method's class loop is guarded by `if (shape[k] > 0)`; its later
`shape[k] == 0` branch is unreachable in this source version. The outer
`else` clears the skipped class membership, while the class-energy slot and
lesion auxiliary maps remain untouched. The harness follows that behavior.
The wrapper ratio step is now included after the exact membership objective.
Its previous objective, `3727.62158203125`, comes from the hash-verified
iteration-5 integration summary recorded in `case.json`; it is the preceding
synthetic run's saved objective, not an independently Java-validated energy.
Registration now runs against the 13 bundled atlas priors and matches the
selected Java single-rigid optimizer trace and positive-support bounds.
The selected wrapper branch and energy ratio also match the Java reference.
Its isolated branch check uses an empty skeleton, so it is not growth evidence.

The saved thinned segmentation's `computeRelations` map and relation table
match the selected Java method for a 14-pass diffusion. This relation check is
pre-growth. The fixture does not include the matching pre-growth ordering and
class-volume arrays needed to run source-faithful `initGrowingLabels` and
`propagateGrowing`; the saved registered ordering/volumes are post-growth and
cannot substitute for those inputs without inventing stage history.

## Run

The full volume is opt-in. Run it under the 16 GB RSS watchdog:

```sh
PYMACRUISE_RUN_P41Z_ITER6_MEMBERSHIP=1 PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/pytest -q --tb=short tests/test_p41z_iteration6_membership_java.py
```

To print the complete comparison report directly (also source-compares the
statistics stage first):

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/compare_p41z_iteration6_membership_java.py
```

The comparator checks raw float32 bits and separately counts IEEE NaNs on
each side. This is a held-fixed synthetic stage comparison, not evidence of a
complete workflow run or real-MRI validation.
