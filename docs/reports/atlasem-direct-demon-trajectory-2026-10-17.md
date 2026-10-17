# AtlasEM direct Demon optimizer trajectory regression — 2026-10-17

## Finding

The direct selected-Java optimizer parity test now rejects unchanged-transform
false passes. It asserts that both the single-scale and three-level pyramid
optimizers produce nonzero transform updates, then requires exact float32
array equality with the selected `DemonToadDeformableAtlas` class.

The pyramid test executes levels 4 and 2 and compares Java and Python at
iteration caps 1, 2, and 5. The five-update case confirms continued optimizer
movement, rather than only matching an initial or first-step state. The
Python API documentation now states the evidence that exists and keeps
descriptor, selected-atlas, and whole-workflow validation gaps explicit.

## Validation

Focused command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_demon_rigid_java_parity.py \
  tests/test_atlas_em_registration_from_descriptor.py --tb=short
```

Result: **9 passed in 6.54 s**.

Full pytest on implementation commit
`9245b50689ff17e317aa22ff2058ed65edc01d09`: **1,651 passed, 4 skipped,
0 failed, 21 expected orientation warnings in 774.40 s**. Peak sampled RSS
was **10.037 GB** under the 16 GB watchdog; it did not terminate the run.
The exact suite command and result are in the
[full-suite report](full-pytest-2026-10-17-atlasem-direct-demon-trajectory.md).

These are synthetic direct optimizer comparisons. They do not validate the
unavailable selected `ibsr8obj` resource set, human AtlasEM convergence, or
complete PLUS workflow equivalence.
