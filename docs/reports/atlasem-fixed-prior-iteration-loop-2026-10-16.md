# AtlasEM fixed-prior iteration loop — 2026-10-16

`run_atlas_em_fixed_prior_iterations` assembles the selected no-field, one-channel Atlas-mode wrapper loop after initialization and transformed-prior generation. Each iteration performs `computeCentriods`, `computeCovariance`, and `computeMemberships` in source order, with fixed atlas priors, then stops when maximum membership change is below the configured threshold or the iteration limit is reached. The selected algorithm does not call `computeAtlas()` again in this loop. Source class SHA-256: `b938bbe5631e41d53bf76c56cce5d7ff745fd000ea2b92e5debf1af7233ba35a`.

A deterministic 5×5×5 three-class phantom with a false one-voxel border and fixed synthetic atlas priors matches the Java excerpt float32 bits for every iteration's energy/change metrics, final centroids, covariances, and posterior memberships. It also matches convergence decision and iteration count. A zero-iteration test verifies that source-selected behavior leaves already initialized state untouched.

Focused command/result: `PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/pytest -q tests/test_atlas_em_iteration.py` — **2 passed in 1.48 s**.

This is a synthetic fixed-prior EM-loop comparison. It does not cover AtlasEM image/atlas registration, the absent `ibsr8obj` files, complete atlas-driven initialization, inhomogeneity correction, or pre-LesionTOADS Dura evolution, and it is not full PLUS validation.
