# TOADS iteration-two membership through registration and wrapper branch (P18cybh)

## Ported behavior

Continue P18cybg's iteration-two exact membership result through the source
float32 energy ratio, actual rigid atlas registration, convergence/iteration
control, and selected stop/continue branch. Derive the preceding loop's state
from the registered cycle: post-thinning segmentation and skeleton, prior
grow-order map, class volumes, availability, refreshed boundary distances,
constructor output labels, and the current rigid transform.

When the selected branch continues, the adapter grows, refreshes relation and
closest-label state, transforms shape priors at the newly registered pose,
and enters the next iteration with thinning. When it stops, it returns the
selected classification result and no next thinning. Registration, control,
growth, relations, shape transformation, and thinning reuse existing kernels.

## Validation

Command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_first_iteration_cycle.py
```

Result: **1 passed, 0 failed, 0 skipped in 20.44 seconds**.

The committed compact RAS case feeds the second exact membership map into the
registration/control adapter by identity and checks the energy ratio against
the selected float32 helper. The test verifies a valid selected branch result;
for a continuing branch it also compares iteration-three thinning against a
direct stage call. Selected Java registration and wrapper-control references
are covered by the related P18cyx and P18cyt tests. This composition is not a
combined Java wrapper run or full TOADS equivalence.

Full-suite evidence is pending. The row adds workflow composition and state
handoff, not newly countable Java algorithm lines.
