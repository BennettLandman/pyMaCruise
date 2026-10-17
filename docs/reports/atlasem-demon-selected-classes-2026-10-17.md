# AtlasEM selected-class registration comparison — 2026-10-17

## Finding

The direct Java harness now derives its class count from the supplied arrays
and accepts an optional per-class registration mask. A generated three-class
case runs the selected `DemonToadDeformableAtlas` pyramid optimizer with
registration mask `[false, true, false]`. The Java and Python six-parameter
float32 transforms match bit for bit and contain nonzero updates.

The test exercises the optimizer's registered-class selection behavior on
synthetic memberships and priors. It does not validate the selected
`ibsr8obj` resource set or AtlasEM behavior on human data.

## Validation

Focused command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_demon_rigid_java_parity.py --tb=short
```

Result: **6 passed in 7.95 s**.

Full pytest on implementation commit
`f727c51431e19706433f262b8838f1be74d58fe4`: **1,652 passed, 4 skipped,
0 failed, 21 expected orientation warnings in 743.94 s**. Peak sampled RSS
was **13.114 GB** under the 16 GB watchdog; it did not terminate the run.
See the [full-suite report](full-pytest-2026-10-17-atlasem-demon-selected-classes.md)
for the exact command.
