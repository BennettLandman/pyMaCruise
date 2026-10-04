# LesionTOADS registered initial thinning — 2026-10-04

## Scope

This synthetic slice composes the bundled 13-class descriptor alignment with
the selected topology/shape resampling, constructor initialization, and first
thinning queue transition. The synthetic subject patch is generated from the
bundled atlas shape/intensity maps; no human MRI is used.

The composition returns the transform/bounds, RAS topology, transformed class
maps, constructor state, and first thinning state in one result. It does not
run normalization, preprocessing-transform estimation, or Java's
post-registration `uncropImages` / `findCroppingBoundaries` / `cropImages`
sequence, and it does not run later LesionTOADS iterations. Those crop stages
are explicitly listed as excluded in the result because they can change the
grid and are tracked next in P39zzzzzzbk.

This result is a synthetic integration of source-checked component adapters,
not a composite Java output comparison or full LesionTOADS equivalence claim.

## Focused reproduction

```bash
PYTHONPATH=src ./.venv/bin/pytest -q \
  tests/test_lesion_registered_initial_thinning.py \
  tests/test_lesion_shapes.py \
  tests/test_lesion_loaded_atlas_first_thinning.py \
  tests/test_lesion_initial_segmentation_state.py
```

Result after adding explicit dtype checks: **21 passed, 0 failed, 0 skipped in
6.72s**. The committed implementation under evaluation is `abe5bc8`
(`abe5bc8d3afdd8e6104cdca4b947b5ba2c71d11b`). The extra dtype assertions were
rerun in this focused group after the full-suite run; they did not change
implementation code. The constructor label map is int8; the thinning queue's
mutable segmentation work array is int16.

The new integration test checks RAS topology orientation, atlas and subject
grid sizes, class order, transformed-template/shape dtypes, constructor and
queue state dtypes, stage order, exact equality to direct calls of the same
component APIs, zero constructor lesion membership, and input immutability.
Java comparisons for the constituent selected methods remain in their own
component fixtures/reports; the integration test does not relabel composition
as a Java match.

## Full-suite evidence

Command:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  env PYTHONPATH=src ./.venv/bin/pytest -q
```

Tested implementation commit: **`abe5bc8d3afdd8e6104cdca4b947b5ba2c71d11b`**.
Result: **1,160 passed, 0 failed, 0 skipped, 11 expected RAS warnings in
887.61s**. Peak sampled RSS was **13.657 GB**; the watchdog did not terminate
the run. The later plan and report commits do not change the implementation
tested here.

This verifies a bounded synthetic integration. It does not establish
post-registration recrop equivalence, whole LesionTOADS equivalence, plus
integration, or real-data validation.
