# Brain-resolution synthetic base run and atlas-audit correction

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

## Full-grid selected Java comparison through convergence — 2026-10-07

After the full-run discrepancy above was found, I ran the ignored
`VandyCoreSyntheticRunExport` harness on the same N4 image used by the base
workflow, with the selected `1e-4` stop threshold, a 12-iteration cap, and a
6 GB Java heap cap. The source NIfTI SHA-256 is
`0d08d8419e31a1034423062e5d26f55e0bd34d3ae551dab127bd99d17359f28a`; the
x-fastest float32 little-endian Java input serialization SHA-256 is
`5729e75b1002609e12f3031d0cac23054a3b25628435a229c7d47c2638c5442b`. Both
implementations processed the internal `(167,207,160)` grid and produced the
native `(172,216,164)` result. The 12-iteration cap reached the same natural
stop condition in Java and Python.

| Iteration | Java energy ratio | Python energy ratio |
|---:|---:|---:|
| 1 | `1.6191568` | `1.61915684` |
| 2 | `0.27991822` | `0.279918224` |
| 3 | `0.012432367` | `0.0124323666` |
| 4 | `0.0014986269` | `0.00149862689` |
| 5 | `0.0010588061` | `0.00105880608` |
| 6 | `0.004032372` | `0.00403237203` |
| 7 | `0.021509102` | `0.0215091016` |
| 8 | `0.0030093095` | `0.00300930953` |
| 9 | `4.1210206e-4` | `0.000412102061` |
| 10 | `1.4254554e-4` | `0.000142545541` |
| 11 | `1.0079462e-4` | `0.000100794619` |
| 12 | `6.968091e-5` | `6.96809075e-5` |

All 12 ratios have identical float32 bit patterns in Java and Python,
including the iteration 6–7 rebound. The final native-grid classifications
are voxel-exact: **0 changed voxels out of 6,092,928**, Dice `1.0` for all 11
labels, and matching x-fastest payload SHA-256
`a83c7e94d28825846ef0011f8bf81cd35fe4ab9d4258c40f44c0d35a3b736ebd`. The
compressed Java-reference NIfTI and provenance are committed as
[`toads-first-java-12iter.nii.gz`](../../tests/data/synthetic/base_workflow_brain_resolution/toads-first-java-12iter.nii.gz)
and
[`toads-first-java-12iter.json`](../../tests/data/synthetic/base_workflow_brain_resolution/toads-first-java-12iter.json).
This is selected VandyCRUISE source equivalence for the first complete
brain-grid TOADS pass; it is not whole MaCRUISE equivalence.

The 12-ratio Java trace also matches the corrected workflow log. The earlier
pre-audit ratios (`0.161843792`, `1.09157264`) came from the missing atlas
conversion audit described below, not from the selected loop implementation.

## Second-pass TOADS comparison through convergence — 2026-10-07

The full-grid second TOADS pass consumes the Dura-stripped T1, matching brain
mask, and audited source-order atlas from the corrected six-stage workflow.
The first Java attempt stopped at iteration 18 because custom diagnostic
instrumentation wrote about 10.1 GB of full-volume snapshots and exhausted
available disk space. This was a harness I/O failure, not a Python/Java
algorithm difference. I removed those disposable generated snapshots and
recompiled the ignored comparison harness with full-array snapshot writes
disabled while retaining its per-iteration Java energy and membership-hash
trace. The selected segmentation, parameter values, and registration behavior
were preserved; Java heap remained capped at 6 GB.

The complete selected Java run and Python workflow both naturally converged at
iteration 33 (threshold `1e-4`). All 33 energy ratios have identical float32
representations:

| Iteration | Java ratio | Python ratio |
|---:|---:|---:|
| 1 | `0.55921155` | `0.559211552` |
| 2 | `0.050367113` | `0.0503671132` |
| 3 | `0.09037805` | `0.0903780535` |
| 4 | `0.06370069` | `0.0637006909` |
| 5 | `0.026589138` | `0.0265891384` |
| 6 | `0.0066712564` | `0.00667125639` |
| 7 | `0.0019709794` | `0.00197097939` |
| 8 | `0.0063308454` | `0.00633084541` |
| 9 | `0.0090437075` | `0.00904370751` |
| 10 | `0.007371022` | `0.0073710219` |
| 11 | `0.0076013017` | `0.00760130165` |
| 12 | `0.005956981` | `0.00595698087` |
| 13 | `0.007843109` | `0.00784310885` |
| 14 | `0.011103654` | `0.0111036543` |
| 15 | `0.012672974` | `0.0126729738` |
| 16 | `0.009928998` | `0.00992899761` |
| 17 | `0.007628008` | `0.00762800779` |
| 18 | `0.007457276` | `0.00745727587` |
| 19 | `0.00448313` | `0.00448312983` |
| 20 | `0.005107686` | `0.00510768592` |
| 21 | `0.0074348385` | `0.00743483845` |
| 22 | `0.0058335043` | `0.00583350426` |
| 23 | `0.0048457794` | `0.00484577939` |
| 24 | `0.0028443993` | `0.00284439931` |
| 25 | `0.0028866357` | `0.00288663572` |
| 26 | `0.0015331173` | `0.00153311726` |
| 27 | `0.00093698537` | `0.000936985365` |
| 28 | `0.00021373917` | `0.000213739171` |
| 29 | `0.0006168019` | `0.000616801903` |
| 30 | `0.0007654224` | `0.000765422417` |
| 31 | `0.0013008404` | `0.00130084041` |
| 32 | `0.00026988387` | `0.000269883865` |
| 33 | `8.3597275e-05` | `8.35972751e-05` |

The native-grid classifications are voxel-exact: **0 changed voxels of
6,092,928**, Dice **1.0** for each of labels 1, 5, 6, 14, 15, 16, 17, 18, 22,
24, and 25. Both x-fastest classification payloads have SHA-256
`9db0db03df770a613c875865295ef56272874043643cefabeb2ee3d5245f35fe`. The compressed Java
reference and provenance are committed as
[`toads-second-java-33iter.nii.gz`](../../tests/data/synthetic/base_workflow_brain_resolution/toads-second-java-33iter.nii.gz)
and
[`toads-second-java-33iter.json`](../../tests/data/synthetic/base_workflow_brain_resolution/toads-second-java-33iter.json).
The reference uses the matching Python NIfTI header to preserve RAS geometry;
its voxel payload is the Java result. Together with the first-pass comparison
above, both TOADS passes in the synthetic base workflow now match selected
VandyCRUISE Java through convergence. This still does not establish Java
equivalence for Dura, MaACE, GVF, TGDM, refinement, or the complete workflow,
and no real MRI was used.

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
12-iteration ratios and classification exactly (see above), including when
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

The selected Java comparison covers one complete first-pass TOADS trajectory
through its natural stop, where all 12 ratios and the final classification
match exactly. The complete six-stage run has not been compared end to end
against Java/MATLAB. No human MRI was used.
