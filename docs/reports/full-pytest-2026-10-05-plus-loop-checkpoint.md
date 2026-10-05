# Full-suite evidence: restartable plus registered-loop segments

- **Tested commit:** `b196efa08e23e092b9c04235cb77915d5da89331`
- **Full command:** `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- **Result:** 1,186 passed, 0 failed, 0 skipped, 11 expected RAS-orientation warnings in 937.68 seconds.
- **Memory:** Peak sampled RSS 15.191 GB; not terminated by the 24 GB watchdog.
- **Focused chain:** `PYTHONPATH=src ./.venv/bin/pytest -q tests/test_plus_lesion_registered_branch.py tests/test_lesion_registered_repeated_iteration.py` — 7 passed in 64.98 seconds.
- **CLI:** `PYTHONPATH=src ./.venv/bin/python -m pymacruise.cli plus-lesion-registered-loop --help` lists the persisted checkpoint interface and its settings.

The compact committed synthetic case writes the continued registered state for iteration 2, validates RAS outputs and the manifest, resumes without recomputation, repairs a corrupted segmentation file, and continues from iteration 2 to iteration 3. This validates bounded selected-loop checkpointing, not natural convergence, Java composite equivalence, or real-MRI behavior. The initial registered-branch checkpoint still computes its own first membership/thinning pass from the upstream registered recrop.
