# LesionTOADS selected 13-class registered iteration — 2026-10-04

## Scope

This slice extends the P39zzzzzzbh transition through the selected `registerShapes`
call and the post-registration stop/classification-or-growth branch. It carries
the updated GM-enabled membership result into rigid registration, then uses the
same membership ratio and the current thinned segmentation/frontier for source
loop control. Atlas shape priors remain in template space for registration;
already-transformed shape values feed the thinning and membership stages.

The callable composes one continuing iteration from a prepared branch state.
It does not build the initial atlas, refresh spatial/distance maps or correction
fields, create the next transformed-shape state, iterate to convergence, or run
the MaCRUISEplus file workflow.

## Validation

Focused command:

```bash
PYTHONPATH=src ./.venv/bin/pytest -q tests/test_lesion13_finite_grow_thin.py
```

Result: **1 passed, 0 failed, 0 skipped in 91.14 seconds** on implementation
commit `73056475da515127c43aa64ade613a15cddd3197`. The regression compares the
composed transform, branch stage order, terminal labels, and lesion subtype
against the manually assembled path. The same test runs pinned Java harnesses
for the component thinning, statistics, membership, rigid-registration, and
stop-control states.

The full suite on that implementation commit passed **1,159 tests, 0 failed,
0 skipped, and 11 expected RAS warnings in 862.69 seconds**; peak sampled RSS
was **14.584 GB** under the 24 GB watchdog. See the
[full-suite report](full-pytest-2026-10-04-lesion-registered-iteration.md).

## Interpretation

This is validated synthetic composition across a registered LesionTOADS
iteration boundary. It demonstrates that existing source-checked operations
can run in the selected wrapper order as one callable. It is not complete
LesionTOADS or MaCRUISEplus workflow equivalence and does not use real MRI.
