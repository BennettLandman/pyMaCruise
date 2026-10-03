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
workflow itself fails on a production atlas.

The Python CLI's prior factor had also been set to `0.2`; the selected Java
wrapper hard-codes `atlasCoupling=0.5`, forwarded as `priorFactor`. The CLI and
diagnostics now use `0.5`. On the matching 35×33×31 grid, the corrected Python
40-cycle diagnostic took **340.489 seconds**, reached its iteration cap, and
ended with ratio **0.177568302** (threshold `1e-4` not reached). Its output was
13³ uint8 RAS labels 1–11, SHA-256
`2b80d97b8f05953c39675473814ac075a0f820343efdce05d17a7a90c4410acd`. The
25-cycle corrected-factor run was faster (**211.256 seconds**) but ended at a
lower ratio `0.056615934`; thus the ratio is nonmonotonic and 40 cycles do not
show improvement. The former Python 50-cycle measurements used the wrong factor
and are superseded. Python and Java match the registration-membership
initialization but not the later membership/statistics trajectory. No
threshold was relaxed.

## Runnable Python diagnostics

The source-sized Python loop probe is:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --padding 10 --output build/python-vandy-padding-10-coupling-05-40.nii.gz \
  --maximum-iterations 40 --convergence-threshold 1e-4
```

It writes per-cycle energies and timings plus a RAS NIfTI output. The Python
state trace is:

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

The selected wrapper settings are `prior_scale=0.2` (`atlasRange` →
`atlasScale` → `MultipleToadSegmentation.priorScale`) and `prior_factor=0.5`
(`AlgorithmMultipleToads.atlasCoupling` → `priorFactor`).

- Focused synthetic state trace, registered-loop, and restartable-base checks:
  **5 passed in 96.18s** after applying the source wrapper's `prior_factor=0.5`
  on the implementation working tree now committed as `e4e0ce2`.
- The 25- and 40-cycle source-sized diagnostics both completed with capped,
  nonconverged results as above.
- Latest full suite before the atlas-range default correction:
  **1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,213.79s**
  on `f3b874b`. No full suite has been run after that correction.
- No real MRI validation was performed.

## Next integration milestone

Compare the first exact membership/statistics pass after the matching
registration state and identify why the objective/ratio oscillates. The current
source-sized Python loop costs about 8.5 seconds per cycle after startup. Do not
infer source convergence from the reduced atlas harness. Once the base run
converges with verified required outputs and stage checks, advance the plus
workflow using synthetic FLAIR-to-T1 registration and the LesionTOADS handoff.
