# Registered LesionTOADS first iteration branch

- Implementation: `run_lesion_toads_registered_first_iteration_branch`
- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_first_iteration_branch.py tests/test_lesion13_finite_grow_thin.py`
- Focused result: **2 passed, 0 failed, 0 skipped in 92.04s**.
- Tested implementation commit: `6488c76d96d7342530a38b23f39f389b471622f8`.

The integration carries P39bm/bn's synthetic first loop state through the
selected `registerShapes` handoff and its subsequent convergence/iteration
branch. The synthetic test enables registration on the asymmetric 9³ RAS
phantom, then exercises the source-configured terminal branch, optional
growth, and classification output. It compares transform and terminal
classification exactly to direct calls of the existing registration and
branch components and verifies registration precedes the convergence check.
The companion 13-class test provides selected Java method-harness comparisons
for component behavior; there is no composite Java run on this fixture.

This is a synthetic integration of source-checked components, not a full
LesionTOADS convergence run, real-image validation, or composite Java output
match. The source-method line-count delta for selected base TOADS is **0**.

Full repository regression evidence is in the [P39bo report](full-pytest-2026-10-04-lesion-registered-first-iteration-branch.md).
