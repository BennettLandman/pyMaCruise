# P41z synthetic iteration-6 thinning: selected Java comparison

## Result

The completed iteration-5 checkpoint was used as the input state for the
following `initThinningLabels` + `propagateThinning` pass. The registered
13-class shape priors were reconstructed from the iteration-5 transform and
support bounds. Their float32 voxel values are bit-exact to the workflow's
saved registered shape volume. The committed input fixture includes the
segmentation, skeleton, 17-channel membership state, 13-channel availability,
ordering field, and shape maps as RAS `.nii.gz` files; its provenance records
source checkpoint and data hashes.

The focused comparison passed against the committed
`Lesion13ClassNextThinningReferenceHarness` with exact serialized equality.
The result comparison covered every available selected-method output:

- 98,325 initialization boundary entries and float32 class/minimum volumes.
- All 267,638 popped events, including processed flags, candidate labels,
  topology outcomes, and inserted neighbors.
- Final segmentation, skeleton, availability, ordering, class volumes, and
  the empty remaining queue.

The validated Python post-thinning arrays are preserved separately in
[`tests/data/synthetic/plus_iteration_6_thinning_result`](../../tests/data/synthetic/plus_iteration_6_thinning_result/README.md)
for the next statistics comparison. The selected source file hashes are
`LesionToadSegmentation.java` `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`,
`AlgorithmLesionToads.java` `556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`,
and the Java harness `101aa12bfe874d937b4a0f8f208471dfe38139aff739f2f1d61cda16a159b83f`.

## Reproduction

```sh
JAVA_TOOL_OPTIONS=-Xmx4g \
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 \
  --interval-seconds 1 --report-interval-seconds 60 -- \
  .venv/bin/pytest tests/test_plus_iteration6_thinning_java_reference.py -q
```

Result: **1 passed in 10.81 seconds**, peak sampled RSS **0.655 GB** under a
16 GB watchdog; Java heap was capped at 4 GB. A prior attempt with the NaN
source-faithful guard enabled was stopped after it repeatedly scanned the full
volume on each queue event. The fixture inputs have no nonfinite membership,
ordering, or shape values, so the passing comparison uses the ordinary finite
path after explicit one-time finiteness assertions. NaN-specific behavior is
outside this case.

## Scope

This is synthetic, fixed-input agreement for the selected Java methods and a
hand-coded reference harness with copied heap/topology helpers. It does not
validate the Java wrapper's complete iteration, registration, membership
updates, or the full MaCRUISEplus workflow; it is not real-MRI validation.
