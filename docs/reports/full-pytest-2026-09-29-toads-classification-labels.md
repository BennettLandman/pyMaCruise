# Full pytest evidence: P18bh TOADS label classification transfer

- Date: 2026-09-29
- Implementation commit tested: `b005798`
- Command: `.venv/bin/python -m pytest -q`
- Result: **535 passed**, 2 expected RAS orientation warnings, in 78.80 seconds.
- Focused classification transfer and validation checks: 6 passed.

The selected-method Java harness and Python helper match the
`MultipleToadSegmentation.updateClassificationFromLabels` label lookup and
changed-voxel count on an asymmetric synthetic 3×2×4 categorical volume. This
is method-level evidence only; it does not validate full TOADS evolution or
complete segmentation output equivalence.
