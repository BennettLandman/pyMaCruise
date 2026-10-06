# Compiled TOADS queue in the brain-resolution synthetic base workflow — 2026-10-06

## Result

The six-stage restartable brain-resolution synthetic base workflow completed
using the Numba-compiled production TOADS thinning path. Both configured TOADS
passes converged before the 50-iteration cap:

| Pass | Iterations | Terminal energy ratio | Stop reason |
|---|---:|---:|---|
| First TOADS | 39 | `6.00817439e-05` | Energy threshold reached |
| Second TOADS | 25 | `8.62507077e-05` | Energy threshold reached |

N4, first TOADS, Dura, second TOADS, CRUISE core (ACE/GVF/TGDM), and MaSeg
refinement stage checks all passed. The run manifest reports `status=complete`.
Peak sampled RSS was **8.314 GB** under a **16 GB** watchdog.

The outputs from this run were compared by stage and output filename against
the earlier brain-resolution run using the Python reference queue at
`build/vandy-brain-resolution-ordering-carry/base-run/run.json`. All **30 of 30
NIfTI outputs have identical SHA-256 hashes**. The `provenance.json` and
`SEG_refine-summary.json` metadata files differ, as expected from distinct
run directories and recorded provenance; the NIfTI results do not.

A separate sequential read of all **30 NIfTI outputs** verified finite voxel
values, RAS axis codes, and manifest shapes. The final `inner.nii.gz`,
`central.nii.gz`, `outer.nii.gz`, and `SEG_refine.nii.gz` are on the native
`(172,216,164)` grid. A second invocation with `--resume` exited 0 and verified
and skipped all six stages.

This is strong exact-regression evidence that the compiled production queue
preserves the prior Python path at full synthetic brain resolution. It is not
a whole-workflow Java/MATLAB output comparison and uses no human MRI.

## Reproduction

From the repository root:

```bash
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/toads-numba-compiled/base-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Append `--resume` to validate cached stage hashes and skip completed work. The
completed resume reported every stage as verified/skipped, exited 0, and peaked
at **0.124 GB** sampled RSS.

## Representative output hashes

| Output | SHA-256 |
|---|---|
| First TOADS classification | `f1e83afda2104ac1738d25457ad3d7ddd1b53886d1f6f6fffa840b6a9b2b3e02` |
| Second TOADS classification | `82afdcc4569329bed6520f32534c049e2cec90621cf0fdf9c636ede47bcbc39a` |
| Inner level set | `2a704d8739afd9425f4602ac8ac93064c584063263322a12a7f7a2067ad0c381` |
| Central level set | `a9caf2a388c2254ad927d4a17a9e4622e3818ec1836ad7fd8145081af2e050d9` |
| Outer level set | `21066e399f92bf6674b6c9e5ac0894d444854d49afb5a7be8b2759216f125540` |
| Refined MaSeg labels | `ac4492009662430d835dcb7c22d7e014e6f55cbe818c8e164c7448bef80baf7f` |

The final CLI invocation printed `Base run complete`; the RSS watchdog exited 0
with `terminated_for_rss_limit=false`. Full-suite and focused test details for
the committed implementation are in the [test report](full-pytest-2026-10-06-numba-toads.md).
