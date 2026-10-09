# Full pytest suite after the labeled folded-cortex milestone

- Tested implementation commit: `eb045989fc9be780186a3ee2d019df225821835e`
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 16 --interval-seconds 0.2 --report-interval-seconds 60 -- ./.venv/bin/pytest -q --tb=short`
- Result: **1,505 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings**
- Pytest duration: **1,831.87 seconds** (30:31)
- Peak sampled RSS: **12.464 GB**, below the 16 GB watchdog limit
- Watchdog exit: 0; process was not terminated for memory

The tested changes add a binary brain-mask-to-signed-bound adapter contract and
its regressions, the labeled folded-cortex fixture and full-profile TGDM Java
reference harness, plus the TOADS Dura-output comparator harness. The follow-up
commit `e0d55ed` fixes float32 serialization in the synthetic comparison
runner and records the generated reports; it does not change the tested TGDM
implementation. Its corrected runner was separately executed under an 8 GB
watchdog, and all saved float outputs were verified as float32.

This full-suite result is regression evidence. It does not supersede the
row-specific source-comparison reports or imply whole-workflow equivalence.
