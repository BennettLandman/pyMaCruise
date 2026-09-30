# Full pytest evidence: P39zzm selected LesionTOADS thinning pass

- Run: `.venv/bin/python -m pytest`
- Result: **606 passed, 2 warnings, in 91.07s**.
- Warnings: the two existing expected warnings from the RAS normalization test, where test NIfTI inputs are deliberately stored as LPS and reoriented.
- Tested implementation commit: `dcff6ce0f975f9f03ae650fc957d4d11ce62dee4`.
- Focused selected-Java pass comparison: `tests/test_lesion_toads_thinning_pass.py` (**1 passed**).
- Synthetic input: [selected LesionTOADS thinning-pass phantom](../../tests/data/synthetic/lesion_toads_thinning_pass/README.md).

The selected-method Java harness and Python implementation match the
324-entry initialized frontier, class-volume sums, all three minimum-volume
values, all 710 popped events (371 processed removals), final segmentation,
skeleton, availability and ordering maps, final class volumes, and exhausted
heap. `volumeRatio=0.4` exercises both the ratio-derived minimum and the
100-voxel floor. The fixture keeps each voxel's best membership at its
current label; it does not exercise label changes/topology rejection in the
composed volume pass, nor the upstream membership and transformed-shape
estimation stages.
