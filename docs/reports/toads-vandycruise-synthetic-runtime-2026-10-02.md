# VandyCRUISE source-runtime synthetic probe

## Result

The selected VandyCRUISE `AlgorithmMultipleToads` core was compiled and run
headlessly against the committed 13³ synthetic T1, its committed brain mask,
and the synthetic 11-class atlas. Java input T1 was zeroed outside that same
brain mask. With padding 10, Java and Python use a 35×33×31 working grid, crop
bounds `(0,14,1,13,1,11)`, and crop transform `(0,0,0,-1,-1,0)`.

The first comparison harness had a malformed descriptor: removing the blank
separator after the final shape section made the legacy Java parser silently
skip `Intensity Atlas`. That invalidated earlier Java energies, hashes, and
convergence claims; they are superseded. A second comparison accidentally used
`relation_limit=1` and optimized exponents 2, while the selected wrapper uses
`relation_limit=3` and descriptor `Optimized factors` (the fixture omits that
section, so Java defaults each factor to 1.0). Those earlier 25/40-cycle
measurements also do not describe the selected source settings.

## Source-aligned inputs and first-pass comparison

With the descriptor delimiter and selected wrapper parameters corrected, Java
and Python match the following full-grid input arrays exactly:

- registration memberships: SHA-256
  `03fafc9db300fbf3f287f4a1e154296cbba7bb3e3db9bdad4093360a269c4d53`;
- initial segmentation: `1691bd780ef73c88c800bc1045c225aef7f657522d80e1017e643a993303becd`;
- transformed shape priors: `a2cdad5a49081d4640b385a86485dfa054f99493b43ad19ad49bf390cc267531`;
- closest-label map: `dedbaca1b44da63eb18e4cb49ba4ffdd47603391e321e66e60d672d1b5352630`.

Startup centroids and class deviations also match after accounting for Java's
`Iscale` display multiplier. The exact membership objective matched at startup
(`88.75139`) and the first loop update (`105.55209`). Membership arrays differed
only slightly: startup max absolute difference `1.1921e-7` (mean
`7.37e-11`); after the first exact update max difference `3.3975e-6` (mean
`6.14e-10`). This is direct source-runtime comparison for the first passes,
not full-loop equivalence.

## Convergence and runtime

The source-aligned Python command was:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --padding 10 --relation-limit 3 \
  --output build/python-vandy-padding-10-source-settings-50.nii.gz \
  --maximum-iterations 50 --convergence-threshold 1e-4
```

It completed 50 iterations in **449.390 seconds**. The objective increased
from about 100 at the first update to over `10^12`; the terminal ratio was
`0.873269856`, so it did not reach `1e-4`. The output is 13³ uint8, RAS,
labels 1–11, SHA-256
`2d13318506469adace848ad41fc7848274ef5ead981155a39a18f74860a66196`.

The direct Java source harness was run with the same masked input, atlas
profile, padding, relation limit, and selected constructor values. It also
fails on this compact synthetic atlas: the membership energy becomes NaN by
iteration 48 of a 50-iteration request. This is a stability failure signal for
this reduced fixture. Because the Java harness also fails, it is not evidence
that the production source workflow converges or fails on its complete atlas.
The synthetic base CLI test separately confirms restartable orchestration,
manifest checks, required output names, and stage checks; its deliberately
short capped TOADS settings do not establish convergence.

## Source behavior newly covered

- The Java descriptor parser needs a blank delimiter after the final shape
  block before the next section is reliably processed.
- `relation_limit=3` controls the selected closest-label diffusion; using 1
  changes the map on more than half of this padded grid.
- The wrapper's optimized CSF/GM/WM factors come from the atlas descriptor and
  default to 1.0 when omitted; `atlasRange=0.2` maps to `prior_scale`, while
  wrapper `atlasCoupling=0.5` maps to `prior_factor`.
- The full initial registration/shape/relation state matches across Java and
  Python; later exact membership values have small float differences that
  amplify into unstable trajectories on this reduced fixture.

These are synthetic and source-runtime findings. They are not a whole-loop
source-equivalence claim or real-MRI validation.

## Validation record

- Focused optimized-factor parser, synthetic state trace, and restartable base
  workflow checks: **4 passed in 81.55s** on implementation commit `d90a735`.
- Full-suite baseline before the recent TOADS wrapper and trace corrections:
  **1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,213.79s**
  on `f3b874b`. No full suite has been run since; tests remain focused because
  the recent changes are limited to TOADS settings, parser behavior, and
  diagnostics.
- No real MRI validation was performed.

## Next integration milestone

Build a larger synthetic TOADS atlas with cortex-like folded anatomy and
source-shaped priors, run both Java and Python against it, and locate the first
post-startup membership difference that grows into instability. Profile and
optimize the 6–7 second exact membership pass per cycle without changing the
source threshold. Once the synthetic base workflow is restartable and
convergent with required outputs and stage checks, advance the plus workflow.
