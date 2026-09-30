# Full pytest evidence: P39zzl LesionTOADS thinning outer loop

- Run: `.venv/bin/python -m pytest`
- Result: **605 passed, 2 warnings, in 81.69s**.
- Warnings: the two existing expected warnings from the RAS normalization test, where test NIfTI inputs are deliberately stored as LPS and reoriented.
- Tested implementation commit: `ede7763567753332a776ddca616e10c7fa13b8e4`.
- Focused selected-Java queue-run comparison: `tests/test_lesion_toads_thinning_queue_run.py` (**1 passed**).
- Synthetic input: [LesionTOADS thinning-run phantom](../../tests/data/synthetic/lesion_toads_thinning_queue_run/README.md).

The repeated queue run matches the selected Java excerpt for all 22 popped
items (20 processed), final segmentation, skeleton, availability, ordering,
class volumes, and the remaining heap. It covers face-neighbor propagation,
one stale item, one minimum-volume-gated item, and strict cutoff behavior: an
item at `-limit` is processed after a previously-above-limit item, then the
next lower-priority item remains queued. Frontier initialization and the
`volumeRatio` minimum-volume formula are still separate stages.
