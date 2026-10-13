# Corrected synthetic PLUS Dura output comparison — 2026-10-12

## Result

On the saved corrected atlas-aligned synthetic PLUS checkpoint, the selected
Java `LesionToadSegmentation.generateDuraRemovalOutputs` loops followed by the
source wrapper's `scaleFloatValue(255)` produce bit-identical float32 CSF, GM,
WM-fill, and WM-mask outputs to Python. Each channel contains 5,471,360
voxels; every channel has zero changed voxels and zero maximum absolute error.
Java and Python used the identical hash-checked registered memberships and
segmentation on the RAS grid `(166, 206, 160)`. Peak sampled RSS was 1.240 GB
under an 8 GB watchdog.

This localizes no error to the Dura-output adapter on this checkpoint. It does
not validate the registered lesion memberships, TOADS/LesionTOADS evolution,
or transfer back to native space; those remain upstream and downstream
comparison points.

## Reproduction

```sh
PYTHONNOUSERSITE=1 PYTHONPATH=src .venv/bin/python \
  scripts/run_with_rss_watchdog.py --limit-gb 8 --interval-seconds 2 \
  --report-interval-seconds 30 -- .venv/bin/python \
  scripts/compare_plus_dura_checkpoint_java.py \
  --workflow-root build/plus-atlas-consistent-corrected-full-2026-10-12 \
  --output-dir build/plus-atlas-consistent-corrected-dura-java-2026-10-12 \
  --report docs/reports/plus-atlas-consistent-corrected-dura-java-2026-10-12.json
```

The report JSON records the checkpoint run signature, input hashes, source
method hash, per-channel voxel errors, and synthetic-only scope. The selected
Java method is provided by
[`LesionToadsDuraCheckpointReferenceHarness.java`](../../tests/java_reference/LesionToadsDuraCheckpointReferenceHarness.java).
