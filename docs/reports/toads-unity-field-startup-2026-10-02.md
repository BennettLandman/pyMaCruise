# Source-faithful TOADS unity-field startup (P18cybd)

## Source correction

Inspection of the selected `VandyCRUISE` wrapper changes the interpretation
of the earlier registered degree-three fit adapters. The wrapper constructs
`InhomogeneityCorrection2` before the first exact-membership pass; the
constructor initializes its multiplicative field to `1.0`. The wrapper then
sets `Eratio=1.0` and `deg=0` before entering the main loop. Its first loop
stage is thinning, and `computeCorrectionField(0)` returns without changing
the field. The first fit is scheduled only in a later iteration if the
previous energy ratio passes the degree gate.

Pinned source files:

- `AlgorithmMultipleToads.java`, SHA-256
  `b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d`
- `InhomogeneityCorrection2.java`, SHA-256
  `6f09f2647a2dd4d03bea6bbdfb2979fc09dcc30834bb3b637bd7a615ed6a13e8`

The earlier P18cyav/cyaw/cyax functions remain useful isolated experiments
for a requested field degree and field-injected membership/thinning. They do
not represent the selected wrapper's initial field schedule.

## Ported behavior

`compute_registered_toads_unity_field_startup` now constructs one unity field
per registered image channel, carries degree zero and initial energy ratio
one, computes the registered initial statistics/membership pass with that
field, then runs the first thinning pass. It does not fit a polynomial during
startup.

## Validation

Focused command:

```text
./.venv/bin/pytest -q tests/test_toads_registered_field_refit.py::test_registered_wrapper_startup_keeps_unity_field_and_degree_zero
```

Result: **1 passed, 0 failed, 0 skipped in 21.86 seconds**.

The committed synthetic RAS atlas/T1/mask inputs run through atlas alignment
and startup. The test checks unity fields, degree and ratio state, stage
order, and membership/thinning equality with direct calls to the selected
registered kernels. This validates synthetic composition; it is not a
combined Java wrapper execution or complete TOADS run.
