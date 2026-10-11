# Full pytest: 2026-10-11 TOADS skeleton state

## Full-suite result

- **Tested implementation commit:** `9cade54` (`Preserve carried TOADS skeleton state`).
- **Command:** `PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- .venv/bin/pytest -q --tb=short`
- **Result:** 1,531 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 801.22 seconds (13:21).
- **Peak sampled RSS:** 13.450 GB; the 16 GB watchdog did not terminate the run.
- **Selected TOADS source-method estimate:** unchanged at 3,511–4,911 of 13,651 lines (26–36%); this cycle added 0 mapped Java method lines. This estimate is not Python runtime coverage or end-to-end completion.

## What changed and what the evidence means

The source constructor allocates a zero-filled Java skeleton grid. `initThinningLabels()` writes the interior and preserves the border. Continuing wrapper calls carry the skeleton returned by `propagateGrowing()` into the following thinning pass. Python low-level calls retain the fresh-constructor zero default; selected continuation adapters now pass the actual carried array. The post-registration Java harness was corrected to consume and preserve the same supplied state rather than resetting its full skeleton, and first-pass reference tests now seed their stated harness state explicitly.

The full suite clears the preceding 19 TOADS thinning/iteration failures on the corrected implementation. It is broad Python regression and selected synthetic Java-method evidence, not whole-workflow equivalence or human-MRI acceptance.

A subsequent read-only source audit identified two thin test-oracle paths: the multichannel next-loop test had supplied a fresh zero skeleton, and the next-default-membership fixture had compared only skeleton interiors. Both were tightened to provide the Java growth skeleton and assert full-grid Java/Python skeleton equality. Their focused tests passed after those follow-up changes. The full suite result above is for `9cade54`, before those test-only oracle updates; the production implementation is unchanged by them.

## Focused validation

- The affected TOADS state/continuation test group passed **18 tests in 35.40 seconds** on the working tree containing implementation commit `9cade54`.
- The strengthened multichannel carried-skeleton Java comparison passed **1 test in 6.71 seconds**.
- The regenerated next-default-membership Java reference and full-grid skeleton comparison passed **1 test in 5.43 seconds**.

Both stricter tests use deterministic synthetic RAS `.nii.gz` fixtures and source-derived Java harnesses. They do not establish full registered-loop convergence.
