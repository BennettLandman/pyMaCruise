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

The corrected one-cycle human replay now matches Java exactly. Segmentation,
skeleton (including the border), and ordering have zero differing voxels over
4,847,040 values each. The final first-iteration classification matches in
all 19,602,960 voxels, with identical SHA-256
`713ed9427c1d3572ed4101af07e04e8fad18ad704721356b4b7e89999fae4651`; every
one of the 11 labels has Dice 1.0 and identical counts. This confirms the
border initialization difference caused the previously observed first-cycle
classification mismatch for this held-fixed in-vivo replay.

The one-cycle energy values also match at float32 precision: Java reports
`Eprev=22871.549` and `energy=16610.645`; Python reports `22871.548828125` and
`16610.64453125`. Python's summary records an energy ratio of `0.31715077`,
1,428,157 thinning pops, and zero growth events.

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

The corrected one-iteration human checkpoint rerun completed in approximately
**6.5 minutes** with **4.789 GiB peak sampled RSS** under a 12 GiB watchdog.
The corrected two-iteration Python replay completed in approximately **6m40s**
at **5.771 GiB peak sampled RSS**; it differs from the existing Java
two-iteration classification at **2,985 of 19,602,960 voxels**. Every label
has Dice above 0.993; the minimum is label 6 at **0.993825749**. Full label
metrics are preserved under the ignored two-iteration checkpoint directory.

An instrumented two-iteration Java replay completed in about **35s** at
**3.519 GiB peak sampled RSS**. Its output is byte-identical to the preserved
Java baseline, so the baseline remains valid. The first-cycle energy values
agree at float32 precision. At the second cycle Java reports `Eprev=16610.645`,
`energy=13595.354`, and ratio **0.1996485**; Python reports ratio
**0.0643060878**. Thus high label Dice does not mean source-equivalent
convergence: a second-cycle energy divergence remains. The current evidence
does not localize that divergence to energy accumulation versus a changed
membership, growth, or registration state; the next comparison should inspect
the saved second-cycle components in call order. The checkpoint reproduction
is kept under ignored `build/`:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/checkpoint/python_checkpoint.py
```

The post-fix run used the repository RSS watchdog:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 1 --report-interval-seconds 60 -- \
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

This is a held-fixed human comparison through two iterations. The first
iteration matches Java exactly, but iteration two differs by 2,985 voxels and
its energy ratio differs substantially. Natural-convergence parity still has
to be rerun from the corrected border state. This does not establish
whole-workflow equivalence or clinical validity.
