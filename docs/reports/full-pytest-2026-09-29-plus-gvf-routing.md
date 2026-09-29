# Full test suite: MaCRUISEplus GVF parameter routing

- Date: 2026-09-29
- Tested implementation commit: `595db90207125764b4b427bac6702aa1cc928ab3`
- Command: `.venv/bin/pytest -q`
- Result: **467 passed in 66.01 seconds**
- Warnings: two expected warnings for synthetic Dura inputs reoriented from LPS to RAS.
- Python: 3.13.2
- pytest: 9.1.1

`run_plus_gvf_stage` passes `PlusPipelineOptions.gvf_lambda` (selected value
0.2) and the MATLAB wrapper's 200-iteration maximum to the existing
`generalized_gvf_image` implementation. On the committed 13³ anisotropic RAS
synthetic scalar phantom, its component-last vector output is exactly equal
to a direct Python adapter call with those settings and differs from the base
default-lambda output. The existing GVF implementation remains the solver;
this step validates configuration routing, not equivalence of SciPy CG to the
source multigrid solver.

Focused plus workflow tests: 12 passed. Fixture provenance and regeneration
instructions are in `tests/data/synthetic/plus_gvf/README.md`.
