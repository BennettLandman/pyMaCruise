# Full pytest evidence: P39zzn topology-changing LesionTOADS thinning

- Run: `.venv/bin/pytest -q`
- Result: **607 passed, 2 warnings, in 89.92s**.
- Warnings: the two expected RAS-normalization warnings from tests that deliberately supply LPS-oriented NIfTI images.
- Tested implementation commit: `5c9b8d38833fe1834c96c5e8cc4fe4003076e570`.
- Focused selected-Java comparison: `tests/test_lesion_toads_topology_thinning.py` (**1 passed**).
- Synthetic input and harness details: [topology-changing thinning phantom](../../tests/data/synthetic/lesion_toads_topology_thinning/README.md).

On the committed 11³ RAS phantom, the selected-method Java harness and Python
pass agree exactly on the 706 popped heap items, candidate labels, 9 accepted
and 69 rejected topology checks, 5 follow-on neighbor insertions, final
segmentation/skeleton/availability/ordering arrays, float32 class-volume bits,
and empty remaining heap. The output differs from the initial labels at 9
voxels. The comparison uses the selected source critical-point LUT and copied
Java heap, with fixed synthetic memberships and shape maps.

This validates the selected thinning stage with topology-changing labels. It
does not validate lesion-membership estimation, atlas sampling, registration,
or the complete LesionTOADS segmentation workflow.
