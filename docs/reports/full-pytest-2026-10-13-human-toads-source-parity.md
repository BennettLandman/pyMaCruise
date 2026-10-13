# Full pytest — 2026-10-13 human TOADS source-parity checkpoint

## Run evidence

- Tested repository commit: `703b8396` (`Document human TOADS orientation parity findings`).
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 5 --report-interval-seconds 60 ./.venv/bin/pytest -q`
- Result: **1,569 passed, 4 skipped, 0 failed**.
- Pytest elapsed time: **705.55 seconds**.
- Peak sampled RSS: **15.764 GB**; the 16 GB watchdog did not terminate the run.
- Warnings: **18 expected orientation warnings**.

The report's own source-orientation explanatory paragraph was added after this
run; no Python implementation or test code changed during the run. The tested
commit includes the Java harness orientation option and the matched human
comparison artifacts. The full suite confirms repository-wide regression
status at that checkpoint; it does not test the human Java comparison itself,
which is reported separately in the paired JSON evidence.

## Interpretation

This suite result is synthetic/regression evidence, not whole-workflow source
equivalence or human anatomical acceptance. The earlier matched human TOADS
case has separate exact Java comparisons. Case 100249 remains an open
orientation-sensitive source mismatch: Java reproduces its saved first-pass
labels in LPS source order, while the current Python RAS route differs. See
[the matched comparison](human-toads-matched-java-comparison-2026-10-13.md)
and [the 100249 replay analysis](human-toads-100249-lps-source-replay-2026-10-13.md).
