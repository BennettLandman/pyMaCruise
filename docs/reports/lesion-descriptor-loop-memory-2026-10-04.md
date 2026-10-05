# LesionTOADS loop-driver memory behavior

- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_repeated_iteration.py`
- Focused result: **1 passed, 0 failed, 0 skipped in 14.09s**.
- Tested loop driver: commit `dbb6a64741a5befc67487a62d8ff0c0d46f723e0` plus this retention-mode change.

The loop driver now stores compact iteration summaries by default rather than
retaining every full-volume shape, membership, and transition result. Its
`retain_iteration_details` option is explicit and intended for bounded
diagnostics. The committed synthetic integration compares the retained
diagnostic results against manual stage composition and verifies default
records omit the full arrays while retaining transforms, energy ratios, and
stop outcomes.

The final branch necessarily retains final scientific state. The synthetic
fixture does not establish peak memory for a brain-resolution LesionTOADS
run. Full repository regression evidence is pending.
