# Registered LesionTOADS continuing first iteration

- Implementation: `run_lesion_toads_registered_first_iteration_branch`
- Focused command: `./.venv/bin/pytest -q tests/test_lesion_registered_first_iteration_branch.py tests/test_lesion13_finite_grow_thin.py`
- Focused result: **3 passed, 0 failed, 0 skipped in 94.22s**.
- Tested implementation commit: `629ec61127b541153c49bb766b16a2b40f9d0466`.

The new synthetic case forces the non-converged branch with a zero energy
threshold. It verifies the selected branch continues, executes growth, and
refreshes relation-map and relation-membership arrays exactly like direct
calls to the existing components. The terminal path remains covered by the
companion case. The existing 13-class continuation test supplies selected
Java source-excerpt/harness comparisons for the component methods; neither
synthetic composition is a composite Java run.

This is branch-control and array-wiring validation on synthetic RAS data, not
a repeated-iteration convergence result, complete LesionTOADS equivalence,
or real-image validation. It adds **0 selected base-TOADS source lines**.

Full repository regression evidence is in the [P39bp report](full-pytest-2026-10-04-lesion-registered-first-iteration-continuation.md).
