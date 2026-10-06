# Brain-resolution synthetic base workflow after TOADS ordering carry — 2026-10-05

## Result

The restartable six-stage synthetic base workflow completed on the generated
brain-resolution input set `(172,216,164)` using the bundled VandyCRUISE atlas.
Both source-configured TOADS passes reached their energy threshold before the
50-iteration cap: the first pass stopped at iteration 39 with ratio
`6.00817439e-05`; the second stopped at iteration 25 with ratio
`8.62507077e-05`. Dura, CRUISE core, and MaSeg refinement completed their stage
checks. The manifest reports `status=complete`, `converged=true`, and
`source_equivalence_claim=false`.

This is a restartable, end-to-end synthetic base run with converged TOADS
passes. It is not a whole-workflow Java/MATLAB comparison or human-MRI
validation. The iteration trajectory is recorded as Python output; only the
selected TOADS transitions in the ordering-carry report have direct Java
reference comparisons.

## Reproduction

From the repository root:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 120 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/vandy-brain-resolution-ordering-carry/base-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Append `--resume` to verify the manifest and skip cached stages. The completed
resume skipped all six stages, exited 0 in 6.6 seconds, and peaked at 0.127 GB
sampled RSS. The full run exited 0 with **6.887 GB** peak sampled RSS under the
16 GB watchdog; it was not terminated for memory use.

## Outputs and checks

Every stage output passed its runner check. The manifest records hashes and
RAS geometry for outputs. The four final products are finite, RAS, and on the
native `(172,216,164)` grid:

| Output | SHA-256 |
|---|---|
| `inner.nii.gz` | `2a704d8739afd9425f4602ac8ac93064c584063263322a12a7f7a2067ad0c381` |
| `central.nii.gz` | `a9caf2a388c2254ad927d4a17a9e4622e3818ec1836ad7fd8145081af2e050d9` |
| `outer.nii.gz` | `21066e399f92bf6674b6c9e5ac0894d444854d49afb5a7be8b2759216f125540` |
| `SEG_refine.nii.gz` | `ac4492009662430d835dcb7c22d7e014e6f55cbe818c8e164c7448bef80baf7f` |

The final refinement stage recorded `source_faithful=true`. The full pytest
suite for the implementation was run separately on commit
`c8a33033e823d9a2fbce5170ee7ebd984008e801`: **1,202 passed, 0 failed, 0
skipped, 11 expected RAS warnings in 1,187.29 seconds**, with 13.539 GB peak
sampled RSS. See the [full-suite report](full-pytest-2026-10-05-toads-ordering-carry.md).

## Scope and remaining validation

The persistent Java `ordering` field is now carried through thinning, growth,
and later iterations. The updated integrated input changes the first pass from
the previous cap-limited run to threshold convergence at iteration 39. The
second pass converges at iteration 25. The outputs and hashes are a synthetic
Python integration reference for future regression checks; do not label them
Java-equivalent. Next, compare source-matched intermediate/final base outputs
against the Java/MATLAB pipeline, then test the selected route on the supplied
real MRI. No human MRI was used in this run.
