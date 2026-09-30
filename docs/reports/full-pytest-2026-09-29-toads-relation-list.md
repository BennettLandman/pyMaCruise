# Full pytest evidence: P18bj TOADS relation-list construction

- Date: 2026-09-29
- Implementation commit tested: `3ebf1248ef891ba5ef772fe4c8232071e2fc412f`
- Command: `.venv/bin/python -m pytest -q`
- Result: **549 passed**, 2 expected RAS orientation warnings, in 77.03 seconds.
- Focused Java comparison and validation checks: 7 passed.

Three selected-method Java synthetic references compare relation IDs and
tables on uniform, planar, and folded multi-label volumes. This is evidence
for relation-set construction only, not the source's subsequent relation-map
diffusion or complete TOADS evolution.
