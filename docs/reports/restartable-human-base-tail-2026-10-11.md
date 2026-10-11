# Restartable human MaCRUISE base tail — 2026-10-11

## Result and scope

The new `base-resume-from-toads` command imported the completed human base
run's `toads-second` outputs after validating the source `run.json`, its T1,
brain mask, MaSeg and atlas hashes, the source-stage command/settings, output
hashes, and RAS/grid/finite properties. The imported `toads-second` pass was
marked converged after 27 iterations with terminal energy ratio
`2.59382614e-05`; its iteration limit was not reached.

**Important source-profile limitation:** the imported upstream run used a bounded tail-development profile recorded in its manifest (`Dura=1`, ACE
diffusion 1, and legacy TGDM count 1). Its TOADS/registration controls were
levels/iterations 4/50, post-refresh 5, registration cap 5, and TOADS max 50.
In particular, it did not execute the selected five-pass Dura removal. The imported tail in this earlier checkpoint also used non-source core values
(WM iso 0.6 and ACE updated-CSF mode 0). Therefore its upstream stages and
tail did not use the selected V3 settings. A separate full-chain human
run is documented in [the follow-up report](full-chain-human-base-2026-10-11.md);
that run used selected high-level stage parameters and TGDM outer count 7, but followed the chosen native-RAS path rather than BET/SPECTRE; lower-level Dura parameter parity remains under review.

The four-stage tail then completed its stages in a separate manifest:
MaSeg-to-TOADS preparation, ACE/GVF/TGDM core, source-faithful segmentation
refinement, and surface/thickness. This tail manifest is linked to the completed
seven-stage upstream run manifest by its hash-verified `toads-second` inputs.
It has 43 output records;
all 43 hashes validate. All 30 NIfTIs are finite and RAS. The run peaked at
**3.101 GiB sampled RSS** under a **12 GiB** watchdog, which did not terminate
it. A same-settings `--resume` verified and reused all four attempt-1 stages;
it completed under a **4 GiB** watchdog at **0.792 GiB** peak RSS. The inner,
central, and outer level-set arrays and source-faithful refined segmentation
are exactly equal to the earlier corrected held-fixed human tail products.

This is a restartable corrected tail linked to a completed upstream seven-stage
run, not a fresh end-to-end execution in this checkpoint and not full-pipeline
Java/MATLAB equivalence. The tail reuses the completed source run rather than
re-executing N4, registration, first TOADS, Dura, or second TOADS. The low raw
MaSeg containment remains; the support fractions are not Dice against an
unavailable historical Java/MATLAB final segmentation. See the
[support localization](in-vivo-maseg-support-localization-2026-10-11.md),
[held-fixed TGDM comparison](in-vivo-fused-input-tgdm-held-fixed-java-comparison-2026-10-11.md),
[held-fixed topology comparison](in-vivo-fused-core-topology-held-fixed-java-comparison-2026-10-11.md),
and [surface-extraction comparison](in-vivo-jist-surface-extraction-held-fixed-2026-10-11.md).

## Reproduction

From the `pymacruise` repository root, the bounded human tail run was:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 12 --interval-seconds 0.5 --report-interval-seconds 30 -- \
  .venv/bin/python -m pymacruise.cli base-resume-from-toads \
  --source-run build/invivo-base-2026-10-10-ras-grid/run.json \
  --maseg build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --output-dir build/invivo-base-2026-10-10-restartable-human-tail-2026-10-11
```

Resume verification:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/python scripts/run_with_rss_watchdog.py \
  --limit-gb 4 --interval-seconds 0.2 --report-interval-seconds 15 -- \
  .venv/bin/python -m pymacruise.cli base-resume-from-toads \
  --source-run build/invivo-base-2026-10-10-ras-grid/run.json \
  --maseg build/invivo-base-2026-10-10-inputs/ras-grid/orig_target_seg_ras_uint16.nii.gz \
  --output-dir build/invivo-base-2026-10-10-restartable-human-tail-2026-10-11 \
  --resume
```

Focused tests on implementation commit `e37a49a`:

```sh
PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl \
  .venv/bin/pytest -q --tb=short tests/test_base_resume_workflow.py
```

Result: **4 passed in 1.00 s**. CLI `--help` and `git diff --check` also
passed. No full test suite was run for this workflow addition. All human image
inputs and outputs remain in ignored `build/`; no human images were committed.
