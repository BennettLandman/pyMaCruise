# Full pytest: restartable plus file workflow

- Tested implementation commit: `4ae3b55` (`Add restartable file-based plus continuation`)
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,158 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 844.58s (14m 04s)**
- Peak sampled RSS: **14.536 GB**; the 24 GB watchdog did not terminate the run.

The four-test increase over the previous repository suite includes the plus FLAIR checkpoint API/CLI tests and the restartable plus continuation API/CLI tests. This validates repository regression behavior for the synthetic code on this commit. It does not establish whole MaCRUISEplus or Java/JIST equivalence, NiftyReg/Aladin equivalence, or real-MRI behavior.
