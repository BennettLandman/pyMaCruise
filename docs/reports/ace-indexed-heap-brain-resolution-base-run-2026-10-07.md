# Brain-resolution synthetic base run after ACE heap correction

## Result

The restartable six-stage synthetic base workflow completed on the current
source tree (`638d781`) at native RAS grid `(172, 216, 164)`. All six stages
passed their output checks, and the manifest contains 30 NIfTI outputs. A
sequential scan confirmed all 30 outputs are finite, retain the expected RAS
axis codes, and match the manifest shapes. `--resume` verified and skipped all
six completed stages.

This is synthetic integration evidence only. It is not whole-workflow Java or
MATLAB equivalence and does not validate the supplied human MRI.

## TOADS stopping behavior

The first TOADS pass completed its configured 50 iterations but did not reach
the energy threshold: terminal energy ratio `0.0316733867`. The second pass
converged at iteration 37 with energy ratio `1.19687047e-05`. Therefore, the
workflow status is `complete`, but the manifest correctly records
`converged=false` and
`result_interpretation=stage-run-complete-but-toads-did-not-meet-energy-threshold`.
The first-pass nonconvergence is an open integration blocker; successful stage
checks do not imply a converged base result.

## Memory and restart evidence

The full command ran under a 16 GB RSS watchdog and peaked at **7.977 GB**;
the watchdog did not terminate it. The second command below reported all six
stages as verified and skipped, then returned a completed run manifest.

Full run:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.1 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/ace-indexed-heap-2026-10-07/base-run \
  --registration-levels 1 --registration-iterations 1 \
  --post-refresh-registration-iterations 1 --registration-max-iterations 1 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Resume check: rerun the same command with `--resume`; observed output:

```text
Resume: verified n4; skipping
Resume: verified toads-first; skipping
Resume: verified dura; skipping
Resume: verified toads-second; skipping
Resume: verified cruise-core; skipping
Resume: verified maseg-refinement; skipping
Base run complete: .../run.json
```

## Reproducibility finding

Compared with the earlier run recorded in
[`toads-compiled-queue-base-run-2026-10-06.md`](toads-compiled-queue-base-run-2026-10-06.md),
29 of the 30 common NIfTI file hashes differ. The earlier run converged in both
passes (iterations 39 and 25); this run did not converge in the first pass.
The N4 output hash is unchanged, making the first TOADS pass the earliest
observed image-output divergence; the remaining 28 differences are downstream.
The current source tree contains subsequent source-order updates to TOADS
float32 arithmetic, registration energy, JAMA solving, approximate thinning
and growing, in addition to the ACE changes. This comparison therefore does
not isolate the ACE change as the cause. Establishing which intervening source
changes alter the brain-resolution trajectory, and whether this first-pass
nonconvergence reflects a remaining port defect or expected source behavior,
is the next integration investigation. The run that produced this report was
executed on `638d781`, before checkpoint signatures included implementation
code. A follow-up now fingerprints sorted Python package sources in the
manifest and each stage dependency signature, so a changed implementation
invalidates old checkpoints instead of silently reusing them. This source
fingerprint was exercised by the focused restartable-workflow suite; details
are recorded after the implementation commit. Comparing the old and current
first-pass logs, iteration 1 energy differs by about `1.05e-7`; by iteration 5,
the values are `0.192783281` versus `0.230697453`. This shows early trajectory
drift, but does not identify which source-order change is responsible.

| Iteration | Prior run energy ratio | Current run energy ratio |
|---:|---:|---:|
| 1 | 0.161843687 | 0.161843792 |
| 2 | 1.09107924 | 1.09157264 |
| 3 | 0.774802864 | 0.772367716 |
| 4 | 1.09418404 | 1.09190392 |
| 5 | 0.192783281 | 0.230697453 |
| 6 | 0.105536781 | 0.0854224488 |
| 7 | 0.0742151216 | 0.0651180595 |
| 8 | 0.0348341092 | 0.0317154676 |
| 9 | 0.0233107414 | 0.0190226417 |
| 10 | 0.00618879683 | 0.00255699526 |

## Evidence boundaries

- Focused ACE plus integration suite on implementation commit `999744c`:
  **25 passed in 40.15 s**. The docs/report commit is `638d781`; no full suite
  was run for this base workflow check.
- After adding source-change invalidation, focused restartable workflow tests:
  **3 passed in 15.66 s** on commit `cce9f59`. The test simulates a changed
  implementation fingerprint, verifies all six stages rerun once, verifies a
  same-version resume skips all six, and checks selective downstream rerun
  after output corruption. It is not a full suite.
- Latest full repository suite remains **1,220 passed, 0 failed, 0 skipped,
  11 expected RAS warnings in 1,463.83 s** on `7ae1dee`, peak sampled RSS
  **12.448 GB** under a 24 GB watchdog. That result predates the current ACE
  source comparison changes and is not represented as validation of this
  workflow rerun.
- Newly mapped source-line delta for this integration rerun: **0**.
- No real imaging data was used.
