# LesionTOADS loop-driver memory behavior

- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_repeated_iteration.py`
- Focused result: **1 passed, 0 failed, 0 skipped in 14.09s**.
- Tested implementation commit: `30d41c8470b80b9f48051d7273e1b758a37fe9fe`.

The loop driver now stores compact iteration summaries by default rather than
retaining every full-volume shape, membership, and transition result. Its
`retain_iteration_details` option is explicit and intended for bounded
diagnostics. The committed synthetic integration compares the retained
diagnostic results against manual stage composition and verifies default
records omit the full arrays while retaining transforms, energy ratios, and
stop outcomes.

The final branch necessarily retains final scientific state. The synthetic
fixture does not establish peak memory for a brain-resolution LesionTOADS
run. Full repository suite: **1,167 passed, 0 failed, 0 skipped, 11 expected
RAS warnings in 907.95s**, peak sampled RSS **13.706 GB** under the 24 GB
watchdog. See the [full-suite report](full-pytest-2026-10-04-lesion-descriptor-loop-memory.md).
