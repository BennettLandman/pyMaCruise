# Brain-resolution synthetic base run after ACE heap correction

## Result

The initial run on `638d781` is retained below as historical evidence. The
corrected replay on implementation commit `9986ac0` completed the restartable
six-stage synthetic base workflow at native RAS grid `(172, 216, 164)`. All
six stages passed their output checks; all 30 manifest NIfTI outputs are finite
and RAS; and the manifest records `status=complete`, `converged=true`, and
`source_equivalence_claim=false`. A second `--resume` verified and skipped all
six stages.

This is synthetic integration evidence only. It is not whole-workflow Java or
MATLAB equivalence and does not validate the supplied human MRI.

## TOADS stopping behavior

With the copied atlas conversion audit restored, the first TOADS pass
converged at iteration 12 with terminal ratio `6.96809075e-05`; the second
converged at iteration 33 with ratio `8.35972751e-05`. Both are below the
configured `1e-4` threshold and stopped before the 50-iteration cap. The old
first-pass ratio and hashes in the historical sections below are superseded:
the former materialized atlas had been interpreted in the wrong voxel order.
Convergence is synthetic integration evidence, not whole-workflow Java or
MATLAB equivalence.

## Memory and restart evidence

The corrected full command ran under a 16 GB RSS watchdog and peaked at
**8.526 GB**; the watchdog did not terminate it. The resume command reported
all six stages as verified and skipped and peaked at **0.127 GB**. The original
attempt peaked at 7.977 GB and is preserved in its history.

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
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1 \
  --resume
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

The measurements in this section describe the **pre-audit** run on `638d781`
and are retained for the investigation history. The follow-up below identifies
the missing atlas conversion audit as the cause of this workflow-path
discrepancy; these measurements must not be used to characterize the corrected
run.

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

The classifications were also compared directly on the shared synthetic grid.
First-pass maps differ at **17,021 of 6,092,928 voxels**. Per-label Dice is
below 0.97 for labels 16 (`0.956931`), 22 (`0.906954`), and 24 (`0.952799`).
Second-pass maps differ at **48,505 of 6,092,928 voxels**; labels 5, 6, 14,
16, 18, 22, and 24 have Dice below 0.97, with a minimum of `0.879295` for
label 14. These are old-Python versus new-Python regression measurements
only. Neither run is the Java/MATLAB reference, so they are not source-
equivalence scores or acceptance decisions.

The complete `638d781` first-pass trace took **8,313.4 seconds (2 h 18 min)**.
It shows recurring decreases followed by energy rebounds; selected adjacent
low/rebound pairs are:

| Low iteration / ratio | Next iteration / ratio |
|---|---|
| 10 / `0.00255699526` | 11 / `0.21966733` |
| 16 / `0.00219647912` | 17 / `0.235409945` |
| 22 / `0.00119246391` | 23 / `0.230245948` |
| 28 / `0.00097300089` | 29 / `0.23061122` |
| 34 / `0.000317806494` | 35 / `0.230662525` |
| 40 / `0.00205375929` | 41 / `0.231708989` |
| 46 / `0.000686858548` | 47 / `0.231917679` |

The prior pre-audit trajectory crossed the stop threshold at iteration 39
(`6.00817439e-05`), while the `638d781` run did not cross it by iteration 50.
That apparent convergence regression resulted from the missing conversion
audit in workflow atlas staging. The audit-aware run converged in both passes;
however, source comparison beyond the first two iterations and end-to-end
Java/MATLAB equivalence remain open.

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

## Full-grid selected Java comparison through two iterations — 2026-10-07

After the full-run discrepancy above was found, I ran the ignored
`VandyCoreSyntheticRunExport` harness on the same masked N4 image used by the
base workflow, with a two-iteration limit and a 6 GB Java heap cap. The Java
input was serialized x-fastest (Fortran order) from
`stages/n4/attempt-01/t1-n4.nii.gz`; its SHA-256 is
`0d08d8419e31a1034423062e5d26f55e0bd34d3ae551dab127bd99d17359f28a`. Both
implementations processed the internal `(167,207,160)` grid and produced the
native `(172,216,164)` result. The Python invocation used a 16 GB RSS watchdog
and peaked at **6.288 GB**.

| Iteration | Java energy ratio | Python energy ratio |
|---:|---:|---:|
| 1 | `1.6191568` | `1.61915684` |
| 2 | `0.27991822` | `0.279918224` |

The two native-grid classifications are voxel-exact: **0 changed voxels**,
Dice `1.0` for all 11 labels, and matching x-fastest SHA-256
`f8d74d62ac9f285b11f102e15314a9769643322dbd17799dde45cf6d37885406`. This is
selected VandyCRUISE source equivalence for two iterations on the
brain-resolution synthetic fixture. It is not 50-iteration convergence or
whole MaCRUISE equivalence.

This trace does not match the first two energy ratios in the earlier full base
workflow log above (`0.161843792`, `1.09157264`). The two-iteration Java and
Python CLI runs used the same N4 file and explicit alignment settings, so the
earlier full-workflow trace initially appeared unreconciled. The root cause
and repair are described below.

## Root cause: materialized atlas lost its conversion audit

The base workflow copies atlas files into the run directory so the run is
self-contained. That copy included the descriptor and NIfTI assets but omitted
`data/atlases/CONVERSION_AUDIT.json`. `load_toads_atlas_data` searches for this
audit near the descriptor; when it was missing, the loader treated the
materialized RAS NIfTI arrays as if their native array order were the original
Java voxel order. The standalone CLI used the repository descriptor and found
the audit, which explains why the same full-grid input matched Java when run
directly and differed when run through the base workflow.

`base_workflow.py` now copies the audit alongside the materialized descriptor,
tracks the source and copied audit hashes in the run manifest and stage
dependencies, and migrates older manifests so their stages are invalidated
and recomputed. A direct loader comparison confirmed that the materialized
atlas topology and all 11 shape arrays now match the source descriptor's
restored arrays exactly. The focused restartable-workflow suite passed:

```text
./.venv/bin/pytest -q tests/test_base_workflow_restartable.py
3 passed in 25.22s
```

The audit-aware standalone `toads-loop` run reproduced the full-grid Java
two-iteration ratios and classification exactly (see above), including when
using the materialized descriptor and writing Dura outputs. The complete
six-stage replay with audit-aware materialization and checkpoint migration has
now completed. It produced 30 finite RAS NIfTIs and passed all stage output
checks. The first and second TOADS passes converged at iterations 12 and 33;
rerunning with `--resume` verified and skipped all six stages. The earlier
50-iteration nonconvergence and output-hash comparison describe the incorrect
pre-audit materialization and are superseded for current workflow assessment.
The corrected run still does not establish whole MaCRUISE equivalence.

## Corrected full replay evidence — 2026-10-07

Implementation commit: `9986ac0`. The command above was run with `--resume`
after manifest migration invalidated all affected checkpoints. The manifest
records six complete stages with attempts N4=3, first TOADS=3, Dura=2, second
TOADS=2, CRUISE core=2, and MaSeg refinement=2. Stage checks passed for 38
recorded products (30 NIfTIs, two JSON sidecars, and six stage logs). A
sequential scan loaded all 30 manifest NIfTIs and verified finite values and
RAS axis codes. The final output list contains 16 products: 12 NIfTIs, two
JSON sidecars, and two stage logs.

The final resume invocation reported `Resume: verified ...; skipping` for
`n4`, `toads-first`, `dura`, `toads-second`, `cruise-core`, and
`maseg-refinement`, then exited 0. Its peak sampled RSS was 0.127 GB. Focused
restart regression after the code fix:

```text
./.venv/bin/pytest -q tests/test_base_workflow_restartable.py
3 passed in 25.22s
```

The selected Java comparison is limited to two full-grid TOADS iterations,
where ratios and classification are exact. The complete six-stage run has not
been compared end to end against Java/MATLAB. No human MRI was used.
