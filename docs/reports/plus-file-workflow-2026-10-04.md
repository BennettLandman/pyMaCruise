# Restartable plus continuation from persisted P39 products

`run_plus_file_workflow` and `pymacruise plus-run` now connect persisted Dura-channel NIfTIs and prealigned segmentation/atlas/brain-mask NIfTIs through three independently verified stages:

1. P39 channel transfer to the segmentation grid and shared multi-atlas-GM crop.
2. Plus boundary/tissue/topology preparation, ACE, GVF, and nested TGDM.
3. TGDM uncrop, inner/central/outer RAS-mm surfaces, and volume/surface thickness.

The run manifest hashes inputs, settings, and stage outputs. `--resume` skips a stage only after its output hashes and RAS/finite NIfTI contracts pass; a corrupted final product invalidates and rebuilds only that stage while valid predecessors are skipped. Each retry writes a new attempt directory, preserving prior files for diagnosis. Images are scanned for nonfinite values in slabs to avoid making full-volume validation copies.

## Synthetic test

Command:

```sh
./.venv/bin/pytest -q tests/test_plus_file_workflow.py tests/test_plus_file_workflow_cli.py
```

The committed 25³ folded-cortex fixture supplies the aligned segmentation, atlas maps, brain mask, and synthetic P39 channel inputs. The test checks all three stage records; RAS full-grid `.nii.gz` level sets and thickness; RAS-mm VTK surfaces; thickness metrics; full valid-resume skipping; changed-settings rejection; and recovery from a corrupted final-stage output while reusing the first two stages. Focused result before commit: **2 passed in 16.91s**; post-commit rerun on `4ae3b55`: **2 passed in 17.08s**; CLI help lists every required input and supported setting.

This is file-based synthetic integration from P39 Dura products and prealigned maps. It does not run FLAIR registration, atlas registration, TOADS, or LesionTOADS. It does not establish NiftyReg equivalence, whole MaCRUISEplus equivalence, or real-data performance. Marching cubes retains its tracked JIST `IsoSurfaceOnGrid` comparison risk.


A full repository pytest run under the 24 GB RSS watchdog passed on tested implementation commit `4ae3b55`: **1,158 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 844.58s**; peak sampled RSS **14.536 GB**. See [full-suite report](full-pytest-2026-10-04-plus-file-workflow.md).


## Manifest interruption hardening

After the full-suite run, the runner was tightened to persist workflow/stage `running` status before validation or computation. This prevents a failed resume from retaining a misleading top-level `complete` value. The synthetic test now sets a previously complete tail stage to `running` to simulate interruption, verifies that valid predecessors still resume, confirms the tail reruns, and parses generated VTK meshes to check RAS+ millimetre coordinates. Focused `tests/test_plus_file_workflow.py tests/test_plus_file_workflow_cli.py`: **2 passed in 17.95s** after commit `371a807`. The repository-wide 1,158-test run above predates this narrow status/validation hardening and remains tied to `4ae3b55`.


## T1-grid contract

The workflow now rejects a segmentation or supplied brain mask whose shape/affine differs from the T1 reference; this makes the requested same-grid brain-mask contract explicit before any crop or segmentation stage. The synthetic file-workflow test includes a shape-mismatched brain mask rejection. Focused `tests/test_plus_file_workflow.py tests/test_plus_file_workflow_cli.py`: **2 passed in 19.53s** on the current working tree; commit-bound rerun pending.
