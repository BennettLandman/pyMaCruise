# Full pytest evidence: P18bg TOADS outer-loop control

- Date: 2026-09-29
- Implementation commit tested: `d1e8721e034524062995039724ed93dfbdc18155`
- Command: `.venv/bin/pytest -q`
- Result: **529 passed**, 2 expected RAS orientation warnings, in 77.99 seconds.
- Focused Java loop-control comparison: 1 passed.

The committed selected-source excerpt covers the float32 energy-ratio formula,
convergence check, iteration ceiling, growing on the stopping iteration when
`doGrowing` is true, final classification update, and relation refresh only
while continuing. This validates the post-membership loop decisions, not the
complete numerical iteration or whole TOADS algorithm.
