# VandyCRUISE source-runtime synthetic probe

## Result

The selected VandyCRUISE `AlgorithmMultipleToads` core compiled and ran
headlessly against the committed 13³ synthetic T1 and 11-class synthetic
atlas. A request for 50 iterations completed at source iteration 45. Instrumented
only in a build-directory copy to print its existing energy ratio, the source
loop reported **1.4491775e-5** at iteration 45, below the selected
`maxDistance=1e-4` threshold. Its terminal segmentation was 13³ with labels
1–11 and SHA-256
`0aae5a5a7b8187674277f3812b32725ce1c03aae0db8097bf5e794355437c9cc`.

This is direct execution of the VandyCRUISE source core on synthetic data. It
does **not** establish Python/source equivalence: the Java core expands its
working image to **35×34×33**, while the Python registered loop uses
**19×17×15** for the same fixture. The constructor also initializes atlas and
transform state through MIPAV-specific preprocessing. The next comparison
must align these source working grids and compare initial memberships and
first-cycle state before interpreting the differing convergence trajectories.

## Runnable Python diagnostic

The Python diagnostic remains reproducible from the repository root:

```bash
./.venv/bin/python scripts/trace_synthetic_toads_loop.py \
  --output build/toads-prior-range-02.nii.gz \
  --maximum-iterations 50 --convergence-threshold 1e-4
```

With the selected `MedicAlgorithmToads08` atlas range, the Python loop took
61.1 seconds, reached the 50-iteration cap, and ended at ratio **0.13936983**;
the threshold was not met. Output SHA-256 was
`668c5ce6987ecd16008551a872650a2cff96cf8c8442932fd577d8f5d3cf3c8a`, shape
13³, uint8, RAS, labels 1–11. This is Python-only synthetic integration
evidence, not Java equivalence or real-data validation.

The run also corrected the default `prior_scale` in
`ToadsRegisteredLoopConfig` from 0.05 to 0.2. The selected Java wrapper
exposes `atlasRange=0.2`; `AlgorithmMultipleToads` forwards that value as
`atlasScale`, and `MultipleToadSegmentation` stores it as `priorScale`. The
change is source-parameter alignment, but it did not resolve Python loop
convergence.

## Java runtime details and limits

The compiled class was the supplied VandyCRUISE `AlgorithmMultipleToads`; the
referenced MIPAV classes were compiled into the ignored `build/` tree against
the 272 JARs in `mipav-src/src`. The harness supplied only a no-op
`MedicUtilPublic` logger shim. Atlas NIfTI arrays were converted into the raw
little-endian files required by this legacy Java descriptor parser. The
generated run log is `build/java-reference/vandy-core-run-50-traced.log`.
These are local investigation artifacts, not committed runtime dependencies or
a reproducible user-facing command. The source run emitted a MIPAV native
library temporary-directory warning but completed without requiring that
native library.

## Validation record

- Focused registered-loop and restartable-base tests: **4 passed in 80.77s**;
  tested implementation commit is recorded in the next project commit.
- Latest full suite before the `prior_scale` config-default correction:
  **1,130 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,213.79s**
  on `f3b874b`. A full suite has not been rerun for this small source-parameter
  default change.
- No real MRI validation was performed.

## Next milestone

Make the source comparison repeatable and align MIPAV preprocessing, working
grid dimensions, atlas coordinates, and wrapper parameters with the Python
startup. Then compare the startup membership/objective and first two complete
outer-loop states. Keep any code change to convergence behavior behind that
state comparison; resume plus-workflow integration after the TOADS source
trajectory is understood.
