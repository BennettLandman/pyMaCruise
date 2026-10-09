# Corrected restartable synthetic base MaCRUISE run

## Run and evidence

This is a synthetic/native-grid workflow integration result on the corrected
implementation commit `b785bff` (`Match retained TOADS registration center`).
It is not a complete MATLAB/JIST workflow comparison and does not validate a
real MRI.

The run completed all seven stages using the bundled atlas and committed
brain-resolution synthetic inputs. Both TOADS passes reached their natural
energy thresholds rather than the iteration cap:

| Stage | Iterations | Terminal energy ratio |
|---|---:|---:|
| First TOADS | 11 | `9.77438758e-05` |
| Second TOADS | 43 | `6.88079308e-05` |

The second pass took longer than the held-fixed replay (43 rather than 24
iterations) because its T1 is the output of the corrected first TOADS and Dura
stages. This is an integrated run, not the same held-fixed input used for the
Java source comparison.

The final manifest records `status=complete`, `converged=true`, seven complete
stages, 46 stage-output records, and 24 final-output records. An independent
audit verified all 46 output SHA-256 hashes. All 31 NIfTI outputs are finite
and RAS-oriented. The four algorithmic VTK mesh products declare RAS
coordinates in millimetres; their vertex counts are 115,518 (inner), 120,726
(central), and 123,828 (outer), with 231,032 / 241,448 / 247,652 faces.
Surface products are retained because the selected thickness calculation
consumes the meshes; visualization remains outside the port.

The same-settings restart command verified and skipped all seven completed
stages. The RSS watchdog observed a 0.121 GB peak for resume under a 4 GB cap.
The full workflow used a 16 GB cap, recorded an 8.500 GB peak, and took about
2 h 49 min from N4-stage creation to the final manifest write. The synthetic
brain-resolution run therefore fits below the user's 16 GB target, with
substantial runtime dominated by the second TOADS pass.

## Commands

Full run:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/pymacruise base-run \
  --t1 tests/data/synthetic/base_workflow_brain_resolution/t1.nii.gz \
  --brain-mask tests/data/synthetic/base_workflow_brain_resolution/brain-mask.nii.gz \
  --maseg tests/data/synthetic/base_workflow_brain_resolution/maseg.nii.gz \
  --atlas-descriptor data/atlases/cruise-atlas-10obj-toads2008.txt \
  --output-dir build/base-workflow-2026-10-09-retained-center \
  --registration-levels 4 --registration-iterations 50 \
  --post-refresh-registration-iterations 5 --registration-max-iterations 5 \
  --toads-iterations 50 --dura-iterations 1 --crop-edge-width 2 \
  --ace-diffusion-iterations 1 --gvf-iterations 200 --tgdm-iterations 1
```

Resume used the identical command and options with `--resume`; the runner
reported `Resume: verified ...; skipping` for N4, both TOADS stages, Dura,
CRUISE core, MaSeg refinement, and surface-thickness. It exited successfully
with a 0.121 GB watchdog peak.

The complete suite was run on commit `b785bff`:

```sh
./.venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 16 --interval-seconds 0.5 --report-interval-seconds 60 -- \
  ./.venv/bin/pytest -q
```

Result: **1,503 passed, 0 failed, 0 skipped, 11 expected RAS-orientation
warnings in 1,764.57 s (29:24); 13.787 GB peak sampled RSS.** The watchdog did
not terminate the test process.

## Interpretation and next work

This closes the corrected synthetic base-workflow integration milestone. The
retained-center fix also restores exact held-fixed Java/Python output for both
TOADS passes, as documented in the [source-fix report](toads-retained-registration-center-source-fix-2026-10-09.md).
The combined workflow run exercises the corrected behavior through all seven
stages, but does not by itself establish whole-workflow numerical equivalence.
Continue with the folded plus workflow and remaining selected-source and
composite checks. The current selected TOADS source-method estimate remains
3,511–4,911 of 13,651 lines (about 26–36%); this integration cycle maps no
additional selected Java method lines.
