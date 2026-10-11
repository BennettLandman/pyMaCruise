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
The previously reported two-iteration Python replay used `relation_limit=1`;
its 2,985 differing classifications, minimum Dice 0.993825749, and energy
ratio 0.0643060878 are superseded and are not source-parity evidence. The
selected Vandy `AlgorithmMultipleToads` wrapper passes `marchingLimit=3.0` to
`MultipleToadSegmentation`, where `limit` controls relation and closest-label
diffusion (Vandy `AlgorithmMultipleToads.java:603–606`;
`MultipleToadSegmentation.java:200, 1409, 1500`). At limit 1, the saved
post-growth closest-label slots 1–3 differ from Java in 772,852, 1,348,596,
and 1,051,653 voxels. A held-fixed limit-3 replay matches all four slots
exactly over all **4,847,040** voxels, with identical per-slot hashes.

The corrected two-cycle Python replay used limit 3 and matches Java exactly
on pre-cycle-2 segmentation, ordering, class volumes, relation map, all 935
relation-membership entries, translated skeleton labels, all 11 availability
channels, all 11 boundary-distance channels, all four closest-label channels,
and terminal cropped classification. The terminal classification SHA-256 is
`dcadcad78f5a4388cbcd5a4cdee7d210470608a66936c9348c84882646ba24ad` for both
implementations. Cycle-2 energy ratios are **0.1996484995** in Python and
**0.1996485** in Java. The corrected replay peaked at **6.280 GiB RSS** under
the 8 GiB watchdog. Java's instrumented two-cycle replay completed in about
35s at **3.519 GiB peak sampled RSS** and is byte-identical to the preserved
Java baseline.

The Python algorithm and all comparison arrays completed; an optional
post-run statistics exporter then stopped on an invalid attribute name
(`class_standard_deviations` instead of `standard_deviations`). The report
uses the arrays saved before that exporter. The config default was corrected
to relation limit 3, consistent with the selected wrapper and CLI. Its focused
default test passed **1 test in 0.69 s**. This held-fixed two-cycle result is
source-equivalence evidence for this human input state, not natural-convergence
or whole-workflow equivalence. The corrected replay and exact command are
under ignored `build/`:

```sh
PYTHONPATH=src .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 8 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python \
  build/invivo-base-2026-10-10-ras-grid/source-selected-dura03/toads-second/java-reference/checkpoint-localize-iter2/out-limit3/python_localize_iter2.py
```

The Java two-cycle checkpoint harness and raw comparison arrays are also
preserved in ignored `build/` under the same `checkpoint-localize-iter2/`
directory. No human data or checkpoint arrays are committed.

Human data, logs, and intermediate arrays remain ignored and are not committed.

## Scope

This is a held-fixed human comparison through two iterations. The earlier
cycle-2 discrepancy came from relation limit 1 in the diagnostic replay; the
selected limit-3 replay matches the two-cycle state and classification exactly.
Natural-convergence parity still needs a corrected-border run. This does not
establish whole-workflow equivalence or clinical validity.
