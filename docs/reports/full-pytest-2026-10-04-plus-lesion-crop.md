# Full repository pytest — plus LesionTOADS normalization and crop

- Tested implementation commit: `a00a23d259387210f10b854f33d0ef985a0b2ba9`.
- Command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Result: **1,173 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 872.95s (14:32)**.
- Peak sampled RSS: **13.820 GB**; the 24 GB watchdog did not terminate pytest.

P41f adds a restartable checkpoint from the P41e registered, N4-corrected,
user-mask-prepared modality bundle through selected LesionTOADS intensity
normalization, topology-inclusive crop-bound calculation, zero-padding, and
crop-transform translation. Its synthetic test compares the saved crop exactly
with direct composition of the existing source-checked normalization and crop
adapters, including its RAS voxel geometry. This does not validate initial
atlas alignment, the iterative LesionTOADS workflow, Java whole-workflow
equivalence, or real MRI. Selected TOADS source-line progress is unchanged.
