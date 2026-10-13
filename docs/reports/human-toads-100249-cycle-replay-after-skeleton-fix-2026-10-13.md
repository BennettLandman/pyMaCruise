# 100249 human TOADS replay after constructor-skeleton fix (2026-10-13)

## Finding

The earlier Python/Java disagreement came from missing constructor state at the first registered thinning pass. The selected Java `MultipleToadSegmentation` constructor seeds its skeleton from the initial classification over the complete grid. Its first `initThinningLabels()` call rewrites only interior voxels, preserving those seeded border values. Python had called the low-level thinning pass without that initial skeleton, so its border remained zero. The registered startup wrapper now forwards the constructor skeleton into the thinning operation.

## Evidence

On case 100249's held-fixed RAS input, selected VandyCRUISE atlas descriptor, mask, and wrapper settings:

- The first-cycle output matches Java exactly: **0 of 10,813,440 changed voxels**, Dice **1.0** for all 11 labels.
- A two-iteration run also matches Java exactly: **0 of 10,813,440 changed voxels**, Dice **1.0** for all 11 labels.
- Captured Java/Python startup state matches exactly before and after the first growth: segmentation, skeleton, available voxels, ordering, memberships, centroids, boundary distances, and closest labels all have zero mismatched values.
- The corrected two-iteration run completed in 577.7 seconds with peak sampled RSS **5.603 GB**, below its 8 GB watchdog.
- Focused regression suite: `./.venv/bin/pytest -q --tb=short tests/test_toads_registered_field_refit.py` — **6 passed in 7.72 seconds**.

Per-label counts for the corrected two-iteration comparison and complete array-state mismatch counts are in [machine-readable output metrics](human-toads-100249-cycle-replay-after-skeleton-fix-2026-10-13.json) and [state metrics](human-toads-100249-cycle-state-after-skeleton-fix-2026-10-13.json).

## Limits and next check

The old 5,445 first-cycle and 18,557 two-cycle mismatch figures were generated before this fix and are superseded for the RAS matched-input comparison. Existing natural-run LPS/RAS differences are a separate orientation/cropping question and remain open. This evidence does not validate natural convergence, the historical LPS run, the whole MaCRUISE workflow, or anatomy quality. The P18cycaq row therefore remains Started. Next, rerun the natural RAS replay against Java and compare with the saved source output; then isolate the historical LPS/RAS preprocessing difference.
