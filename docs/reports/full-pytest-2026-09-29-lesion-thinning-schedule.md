# Full pytest evidence: P39zzj LesionTOADS neighbor scheduling

- Run: `.venv/bin/python -m pytest`
- Result: **602 passed, 2 warnings, in 95.25s**.
- Warnings: the two existing expected warnings from the RAS normalization test, where test NIfTI inputs are deliberately stored as LPS and reoriented.
- Tested implementation commit: `f96689b32d6c5589ed0e921fd176b8affdac4eeb`.
- Focused Java/heap comparison: `tests/test_lesion_toads_thinning_schedule.py` (**1 passed**).
- Synthetic input: [LesionTOADS thinning-neighbor phantom](../../tests/data/synthetic/lesion_toads_thinning_schedule/README.md).

The selected LesionTOADS neighbor loop matched Python for every insertion,
availability change, and Java heap pop order. Cases cover ordinary and
lesion processing classes, both lesion scoring modes, six faces, equal-ratio
heap ties, strict `-limit` boundaries, equal/below class-volume cutoffs, and
an unavailable neighbor. This validates neighbor scheduling only, not a
complete LesionTOADS thinning pass.
