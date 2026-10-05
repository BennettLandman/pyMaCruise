# LesionTOADS descriptor loop convergence exit

- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_repeated_iteration.py`
- Focused result: **1 passed, 0 failed, 0 skipped in 10.48s**.
- Implementation under test: registered loop driver from commit `dbb6a64741a5befc67487a62d8ff0c0d46f723e0`.
- Test and report commit: `7e9677a842f306d6e7dd27b0e56eb25d0b4a1ca7`.

The integration test reuses the asymmetric 9³ RAS lesion phantom and sets a
high synthetic energy threshold. It verifies the loop driver stops after its
first new transition with `stop_reason="convergence"`, reports a convergence
stop, and does not misclassify the stop as the iteration cap. A separate run
in the same test verifies the three-update iteration-limit route.

This is synthetic branch evidence for the Python source-order composition,
not a full Java workflow comparison, natural convergence behavior, or
real-image validation. Full suite on the test commit: **1,167 passed, 0
failed, 0 skipped, 11 expected RAS warnings in 874.45s**, peak sampled RSS
**13.384 GB** under the 24 GB watchdog. See the [full-suite report](full-pytest-2026-10-04-lesion-descriptor-convergence-stop.md).
