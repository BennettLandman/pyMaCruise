# In-vivo TOADS first-thinning border state — 2026-10-11

## Result

On the first human-data thinning pass, Java and Python produce identical
segmentation and ordering arrays. The skeleton arrays also match throughout
the interior. Their only difference is all **170,872 border voxels**: the
selected Java wrapper leaves the border at zero because
`MultipleToadSegmentation` allocates `skeleton` as a Java `byte[][][]` and
`initThinningLabels` assigns only the interior (`MultipleToadSegmentation.java`,
lines 274–275 and 505–541). The Python initializer had filled the entire array
with its removed-voxel sentinel `-1` before assigning the interior.

The Python default now matches Java's zero-initialized constructor state.
Callers that replay a separately pre-seeded Java state can still pass that
initial skeleton explicitly. The existing small Java harness seeds its border
with `-1`, so its test now states that setup explicitly instead of treating it
as the default constructor state.

The one-cycle energy values also match at float32 precision: Java reports
`Eprev=22871.549` and `energy=16610.645`; Python reports `22871.548828125` and
`16610.64453125`. The previously observed 3,234 classification differences
after one iteration have **not yet been rerun** with the corrected border.
This identifies and fixes the earliest observed state mismatch, but does not
yet show whether it caused the later classification differences.

## Validation

Focused tests ran on implementation commit `c3c0478`:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/pytest -q --tb=short \
  tests/test_toads_initialized_thinning_pass.py \
  tests/test_toads_thinning_queue.py \
  tests/test_toads_registered_field_refit.py \
  tests/test_toads_registered_first_iteration_cycle.py
```

Result: **13 passed in 9.57 s**. The full suite was not run.

The human checkpoint reproduction is kept under ignored `build/`:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/checkpoint/python_checkpoint.py
```

The corresponding temporary Java harness invocation used a one-iteration
ceiling and exported its post-thinning segmentation, skeleton, and ordering
to
`build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/checkpoint/out/`.
Its runnable command was:

```sh
CP="build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/checkpoint/classes:build/java-reference/toads-selected-wrapper-harness:$(cat build/java-reference/toads-selected-wrapper-harness/classpath.txt)"
java -Xmx12g \
  -Dtoads.checkpoint.dir=build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/checkpoint/out \
  -cp "$CP" ToadsWorkflowSourceHarness \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/first-input-masked.raw \
  /Users/landmaba/dev/mipav/toads-cruise-scm-2026-09-26/VandyCRUISE/src/Atlas/cruise-atlas-10obj-toads2008.txt \
  208 309 305 1 1 1 1 \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/checkpoint/out/classification.raw
```

Human data, logs, and intermediate arrays remain ignored and are not committed.

## Scope

This is a held-fixed first-thinning state comparison on one in-vivo input. It
does not establish whole-workflow equivalence, final natural-convergence
parity, or clinical validity. The next check reruns the human Python
first-cycle comparison with the corrected border and measures whether the
later segmentation difference decreases.
