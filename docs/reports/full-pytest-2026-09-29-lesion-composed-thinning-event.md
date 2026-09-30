# Full pytest evidence: P39zzk composed LesionTOADS thinning event

- Run: `.venv/bin/python -m pytest`
- Result: **604 passed, 2 warnings, in 125.47s**.
- Warnings: the two existing expected warnings from the RAS normalization test, where test NIfTI inputs are deliberately stored as LPS and reoriented.
- Tested implementation commit: `60ecc975ab4b46cefa009fdd621b56121001c35b`.
- Focused composed-event checks: `tests/test_lesion_toads_composed_thinning_event.py` (**2 passed**).
- Synthetic input: [composed LesionTOADS event phantom](../../tests/data/synthetic/lesion_toads_composed_thinning_event/README.md).

The Java-excerpt comparison covers heap pop, lesion-aware target selection,
topology acceptance, center segmentation/skeleton/availability/ordering and
ordinary-class volume updates, six neighbor insertions, and remaining heap
order. Python controls also confirm that a stale item and a class at its
minimum volume leave image and queue state unchanged after the pop. This is
one composed event, not repeated queue evolution or full segmentation
equivalence.
