# Full pytest evidence: P18bk TOADS relation-map diffusion

- Date: 2026-09-29
- Implementation commit tested: `19ae617f0fca849b11a60211758792acedc25ec1`
- Command: `.venv/bin/python -m pytest -q`
- Result: **554 passed**, 2 expected RAS orientation warnings, in 88.89 seconds.
- Focused Java comparison and validation checks: 5 passed.

Selected Java excerpts for relation-set construction and `computeRelations`
match Python's final map and loop-pass count on synthetic uniform, planar, and
folded cases. This covers exact-membership relation propagation only; it does
not validate approximate membership, topology evolution, or full TOADS.
