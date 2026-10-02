# Registered TOADS loop runner (P18cym)

## Ported behavior

Add a callable Python orchestrator that starts from the selected registered
atlas run, initializes unity correction fields and exact memberships, performs
initial thinning, then carries each exact statistics/membership update
through registration and the continuing or terminal branch. Continuing cycles
reuse their transformed shape, relation, availability, and thinning state.
Terminal classification and changed-voxel count are returned with all
intermediate exact iterations and registered cycle records.

The runner composes existing source-checked operations; it adds no generic
segmentation substitute or new numerical kernel. It requires at least two
loop iterations because the existing first-cycle adapter exposes only the
continuing first-iteration path. PSF bias and the full two-pass MaCRUISE/dura
workflow are outside this slice.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 55.86 seconds**.

On the compact synthetic RAS atlas/T1/mask run with a four-iteration limit, the
runner produces four exact membership passes and four registered branches.
Its final membership maps, terminal classification, and changed-voxel count
match the independently assembled stage-by-stage path in the same test. This
validates the selected synthetic configuration only; it is not execution of
the Java wrapper, a full-suite result, or full TOADS equivalence. Full-suite
evidence is pending; P18cym remains Started.
