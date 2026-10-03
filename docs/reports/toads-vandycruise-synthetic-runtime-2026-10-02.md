# VandyCRUISE source-runtime synthetic probe

## Result

The selected VandyCRUISE `AlgorithmMultipleToads` core was compiled and run
headlessly against the committed 13³ synthetic T1, its committed brain mask,
and the synthetic 11-class atlas. Java input T1 was zeroed outside that same
brain mask. With selected padding 10, Java and Python both produce a 35×33×31
working grid, crop bounds `(0,14,1,13,1,11)`, and crop transform
`(0,0,0,-1,-1,0)`.

The first comparison was invalid: the temporary raw atlas descriptor had its
blank separator removed. The legacy Java descriptor parser then skipped the
`Intensity Atlas` section without reporting an error. Restoring the separator
and a following section sentinel made the Java profile visible and produced
source initial registration-membership SHA-256
`03fafc9db300fbf3f287f4a1e154296cbba7bb3e3db9bdad4093360a269c4d53`, matching
Python's pre-registration membership array exactly. Earlier Java energies,
iteration hashes, output hashes, and convergence claims from the malformed
atlas are superseded.

The corrected compact-atlas harness does **not** converge: Java membership
energy becomes NaN by iteration 48 of a 50-iteration request. This is evidence
about this synthetic harness and reduced atlas, not evidence that the source
workflow itself fails on a production atlas. Python's matching-grid 50-cycle
diagnostic took 428.890 seconds, reached its iteration cap, and ended with
ratio `0.375456274` (threshold `1e-4` not reached). Its output was 13³ uint8
RAS labels 1–11, SHA-256
`d73de44b67de5beb4c40d0299365fdba3fbe40d8e229eb7ae4ed0a18c7a32596`.
Python and Java therefore match the registration-membership initialization but
not the later membership/statistics trajectory. No threshold was relaxed.

## Runnable Python diagnostics

The Python state trace is:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_state.py \
  --padding 10 --maximum-iterations 6 \
  --output build/python-vandy-padding-10-state-6.json
```

The slower full loop diagnostic is:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --padding 10 --output build/python-vandy-padding-10-50.nii.gz \
  --maximum-iterations 50 --convergence-threshold 1e-4
```

Both are build artifacts. The JSON labels its evidence as Python-only
synthetic diagnostics, with source equivalence and real-data validation false.
The Java harness and runtime build artifacts are local ignored files under
`build/java-reference/`; they are not a reproducible project command yet.

## Source behavior newly established

- The legacy descriptor parser needs a blank delimiter after the final shape
  block for the next section header to be processed. A malformed descriptor can
  silently omit the intensity profile.
- With that parsing corrected, Java's initialized atlas-registration
  memberships match Python's on the shared padded synthetic grid.
- The first subsequent Java/Python membership/statistics states still differ;
  locating that first difference is the active convergence blocker.

These are direct source-runtime and synthetic observations. They are not full
source equivalence or real-MRI validation.

## Selected wrapper setting and validation record

The registered loop config default `prior_scale` is `0.2`: the selected
`MedicAlgorithmToads08` profile exposes `atlasRange=0.2`, forwarded through
`atlasScale` to `MultipleToadSegmentation.priorScale`.

- Focused trace and registered-first-cycle checks: **4 passed in 19.95s** on
  the current working tree before this report/README edit. The code commit
  `9285095` contains the atlas-range default; the trace change remains
  uncommitted at the time of this report.
- Focused registered-loop and restartable-base checks: **4 passed in 80.77s**
  on `9285095`.
- Latest full suite before the atlas-range default correction:
  **1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,213.79s**
  on `f3b874b`. No full suite has been run after that correction.
- No real MRI validation was performed.

## Next integration milestone

Compare the first exact membership/statistics pass after the matching
registration state, then fix or explicitly bound the cause of the Python
trajectory divergence and 7-minute synthetic runtime. Do not infer source
convergence from the reduced atlas harness. Once the base run converges with
verified required outputs and stage checks, advance the plus workflow using
synthetic FLAIR-to-T1 registration and the LesionTOADS handoff.
