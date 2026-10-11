# Full pytest: 2026-10-11 V3 boundary integration

## Result

- **Tested implementation commit:** `f3b11c6` (`Update human base status with boundary results`). Subsequent commits before and during the run changed only `plan.md`; no Python implementation changed.
- **Command:** `PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl .venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 1 --report-interval-seconds 60 -- .venv/bin/pytest -q --tb=short`
- **Outcome:** 1,512 passed, 19 failed, 0 skipped, 11 warnings; pytest elapsed time 781.81 seconds (13:01).
- **Peak sampled RSS:** 13.535 GB under the 16 GB watchdog; the run was not terminated for memory.
- **Selected TOADS source-method estimate:** unchanged at 3,511–4,911 of 13,651 lines (26–36%), with 0 mapped selected Java method lines added by the boundary wiring. This is a coarse source-method estimate, not Python runtime line coverage or end-to-end completion.

## Failure scope

All 19 failures are in selected TOADS thinning/iteration comparisons against committed Java or Python snapshots. They cover constructor border/skeleton state, downstream growing/thinning behavior, and carried-degree energy ratios. No ACE-boundary, base-workflow restart, human-tail, or non-TOADS test failed. Two representative cases were rerun individually and reproduced:

- `tests/test_toads_default_first_thinning.py::test_selected_exact_initial_memberships_feed_first_default_thinning_pass`
- `tests/test_toads_continuing_iteration_slice.py::test_folded_asymmetric_continuing_slice_matches_fieldwise_snapshot`

The failures were introduced after the previous full-suite run by `c3c0478` (`Match Java TOADS constructor border state`), which changed the omitted-initial-skeleton default from the empty sentinel to Java's zero-initialized byte-array border state. The full suite establishes the regression in the current test/reference contract; a source/harness audit is in progress to decide which comparison fixtures must supply their original explicit initial skeleton and whether any production behavior needs correction. Do not treat this run as a passing full-suite gate or as evidence that the boundary integration caused the TOADS failures.

## Boundary integration evidence remains separate

The boundary change itself has focused synthetic coverage (19 passed in 26.36 seconds on `ba68467`) and a corrected human tail run from a verified TOADS checkpoint. See [ACE boundary integration](base-ace-boundary-integration-2026-10-11.md) and [corrected human tail](in-vivo-base-tail-with-v3-boundary-2026-10-11.md). The human ribbon support remains below its acceptance target, and there is no matched historical full-pipeline reference.
