# LesionTOADS registered initial membership composition

- Implementation: `run_lesion_toads_registered_initial_membership`
- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_initial_membership.py tests/test_lesion_relations.py tests/test_lesion_statistics_membership_bridge.py tests/test_lesion_channel_weights.py tests/test_lesion_standard13_voxel.py tests/test_lesion_standard13_gm_flag.py tests/test_lesion_standard13_volume.py tests/test_lesion_ventricle_distance.py tests/test_lesion_gm_distance.py tests/test_lesion_brainstem_distance.py tests/test_lesion_initial_loop_entry.py tests/test_lesion_initial_thinning_composition.py`
- Focused result: **30 passed, 0 failed, 0 skipped in 11.67s**.
- Full-suite command: `./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`
- Full-suite result: **1,163 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 872.21s** on tested commit `107cbcf2dafc69a24e92a84ec361a1af41de82e0`.
- Peak sampled RSS: **13.805 GB**; the 24 GB watchdog did not terminate pytest.

The slice starts from Java constructor arrays and centroid/atlas state already
on the post-registration RAS recrop grid. It runs relations, variances,
channel weights, the three selected distance maps, initial standard-13
memberships with `GMflag=true`, loop-entry initialization, and the first
thinning handoff. It records `atlas.updateShapeRegistration` and optional
inhomogeneity correction as supplied/disabled upstream stages. Registration
memberships are never reused as segmentation memberships.

The committed synthetic 9³ RAS fixture checks geometry and byte digests,
default constructor membership values, finite outputs, GM-enabled stage order,
and the zero-lesion-deviation path with the source IEEE mode enabled. It then
independently recomposes the existing relation, statistics, channel-weight,
distance, volume-membership, and thinning components and compares those
results exactly. The component suites include their own pinned Java
source-excerpt/harness evidence. This is **selected-component source evidence
plus synthetic Python composition**, not a composite Java output comparison,
complete LesionTOADS equivalence, or real-image validation. The synthetic
first-thinning queue limit is deliberately high to exercise queued events; it
does not claim to reproduce the selected clinical run's convergence.

Implementation changes only compose already represented selected methods;
the selected TOADS source-line accounting delta is **0**. Full repository
suite evidence is regression evidence for this commit, not whole-workflow
equivalence or real-image validation.
