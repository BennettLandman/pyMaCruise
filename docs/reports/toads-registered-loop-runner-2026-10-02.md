# Registered TOADS loop runner (P18cym)

## Ported behavior

Add a callable Python orchestrator that starts from the selected registered
atlas run, initializes unity correction fields and exact memberships, performs
initial thinning, then carries each exact statistics/membership update
through registration and the continuing or terminal branch. Continuing cycles
reuse their transformed shape, relation, availability, and thinning state.
Terminal crop-space classification and changed-voxel count are returned with all
intermediate exact iterations and registered cycle records. The result also
restores classification to the original native RAS voxel grid using the
selected x-fast uncrop buffer and symmetric crop padding.

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

Result: **1 passed, 0 failed, 0 skipped in 55.94 seconds**.

On the compact synthetic RAS atlas/T1/mask run, the manually assembled
degree-zero path continues to match direct selected-kernel calls. The runner
also executes a four-iteration forced-ramp case with degrees `[0, 1, 2, 3]`,
refitting at each positive degree and reaching terminal classification. The
selected field-fit degrees have separate Java formula-harness evidence in
P18cybb/P18cyba; this runner test checks the carried schedule and terminal
return contract, not each full wrapper state against Java. This is not
execution of the Java wrapper, a full-suite result, or full TOADS equivalence.
Full-suite evidence is pending; P18cym remains Started.
