# Persisted LesionTOADS Dura-removal outputs

**Implementation commit:** `2b1c92ceddc7400016f1d0ccb70a32b2d9740419`

P41r extends `pymacruise plus-full-run` after P41m. It reads the persisted 17-channel membership array, class-index segmentation, and source-ordered 13 class names, then uses the existing source-checked `generate_lesion_dura_removal_outputs` adapter. It writes float32 CSF, GM, WM-fill, and hard-WM-mask `.nii.gz` images on the registered recrop grid with an input/output-hash manifest and restartable repair. Native T1-grid transfer is not part of this checkpoint.

## Validation

- Focused command: `PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl ./.venv/bin/pytest -q tests/test_plus_full_workflow.py` — **1 passed, 0 failed, 0 skipped in 18.37s** on implementation commit `2b1c92ceddc7400016f1d0ccb70a32b2d9740419`.
- The committed compact synthetic plus atlas case runs the `plus-full-run` CLI to terminal LesionTOADS outputs. Each saved channel equals a direct call to the selected output adapter on persisted P41m state, voxel-for-voxel. The test also verifies class/channel order, float32 dtype, RAS affine and dimensions, hashes, regular resume, corrupt P41m terminal-output repair, corrupt Dura-channel repair, and unchanged valid P41m checkpoint hash.
- Memory-watched full suite: **1,193 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,000.66s**. Peak sampled RSS was **14.374 GB** under the 24 GB watchdog; the run was not terminated. Command: `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`.

## Scope

This validates synthetic Python file-workflow composition around a selected source-checked adapter. It is not Java composite equivalence or real-MRI validation. The outputs remain on the registered recrop grid; transfer to native T1/MaSeg space and downstream Dura-to-CRUISE surfaces are the next integration work. No human MRI was used.

See [matching full-suite evidence](full-pytest-2026-10-05-lesion-dura-output.md).
