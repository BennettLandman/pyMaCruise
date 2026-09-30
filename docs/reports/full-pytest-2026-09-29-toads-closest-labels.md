# Full pytest evidence: P18bi TOADS closest-label propagation

- Date: 2026-09-29
- Implementation commit tested: `5f7a1fad86af0be938d8493a17b179b35d6c2e6b`
- Command: `.venv/bin/python -m pytest -q`
- Result: **542 passed**, 2 expected RAS orientation warnings, in 78.03 seconds.
- Focused Java comparison and validation checks: 7 passed.

The selected `computeClosestLabels` method excerpt and Python implementation
match every four-slot output on the committed asymmetric 9×8×7 phantom for
diffusion limits 0, 1, and 3. This is method-level evidence only; relation-map
construction and full TOADS evolution remain open.
