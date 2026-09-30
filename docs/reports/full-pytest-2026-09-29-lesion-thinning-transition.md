# Full pytest evidence: P39zzi LesionTOADS thinning event

- Run: `.venv/bin/python -m pytest`
- Result: **601 passed, 2 warnings, in 80.36s**.
- Warnings: the two existing expected warnings from the RAS normalization test, where test NIfTI inputs are deliberately stored as LPS and reoriented.
- Tested implementation commit: `f2c087d331d8199bde9153b18f0d30d19f7c8376`.
- Focused Java comparison: `tests/test_lesion_toads_thinning_transition.py` (**1 passed**).
- Synthetic input: [LesionTOADS thinning-event phantom](../../tests/data/synthetic/lesion_toads_thinning_transition/README.md).

The Java excerpt compares lesion-aware best-membership selection, direct and
shape-fallback candidate choice, LUT topology acceptance, critical-volume
eligibility, segmentation and skeleton updates, availability and ordering,
and float32 class-volume decrement. This is one local removal event, not a
complete LesionTOADS thinning pass or whole segmentation validation.
