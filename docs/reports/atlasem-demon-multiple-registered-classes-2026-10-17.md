# AtlasEM multiple registered-class comparison — 2026-10-17

## Finding

The direct Java harness derives class count from the shape-prior payload and
accepts an explicit per-class registration mask. A generated three-class
case compares the selected `DemonToadDeformableAtlas` pyramid optimizer and
Python adapter with masks `[false, true, false]` and `[false, true, true]`.
The transforms are nonzero and match bit for bit across Java/Python in both
cases. Enabling the second registered class therefore exercises how its
shape contributes to the shared rigid transform.

The arrays are synthetic, 24³, and keep tested transformed support within
the source image domain. This comparison validates class selection and
multi-class optimizer accumulation, not selected descriptor assets, AtlasEM
EM/Dura behavior, or whole PLUS workflow equivalence.

## Validation

Focused direct source comparison:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_demon_rigid_java_parity.py --tb=short
```

Result: **7 passed in 7.06 s**.

Full pytest on implementation commit
`202939fd224e87298cd1f6cd931292f3e454f38c`: **1,653 passed, 4 skipped,
0 failed, 21 expected orientation warnings in 769.95 s**. Peak sampled RSS
was **12.158 GB** under a 16 GB watchdog; it did not terminate the run. See
the [full-suite report](full-pytest-2026-10-17-atlasem-multiple-registered.md).
