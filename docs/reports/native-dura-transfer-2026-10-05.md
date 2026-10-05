# Native T1-grid transfer of LesionTOADS Dura channels

**Implementation commit:** `2c32fd89ae62d75ffa3e445e8f7079e59b620d10`

P41s extends `pymacruise plus-full-run` after P41r. It transfers the selected GM, WM-fill, and hard-WM-mask Dura channels from the registered recrop grid to the original T1 reference grid, then writes float32 compressed NIfTIs in RAS. The transfer uses the existing `transfer_to_reference_space` physical-affine, nearest-neighbor adapter and has its own manifest, input signature, output hashes, resume validation, and corrupt-output repair.

## Validation

- Focused command: `PYTHONPATH=src MPLCONFIGDIR=/private/tmp/pymacruise-mpl ./.venv/bin/pytest -q tests/test_plus_full_workflow.py` — **1 passed, 0 failed, 0 skipped in 18.67s** on implementation commit `2c32fd89ae62d75ffa3e445e8f7079e59b620d10`.
- The committed synthetic plus fixture runs the CLI through P41s. All three stored native-grid outputs match direct calls to `transfer_to_reference_space` voxel-for-voxel. Checks cover native T1 shape/affine, RAS orientation, float32 type, source call order, hash metadata, valid resume, and repair of a corrupted native output without changing the prior Dura-output manifest.
- Memory-watched full suite: **1,193 passed, 0 failed, 0 skipped, 11 expected RAS warnings in 1,026.32s**. Peak sampled RSS was **14.479 GB** under the 24 GB watchdog; the run was not terminated. Command: `PYTHONPATH=src ./.venv/bin/python scripts/run_with_rss_watchdog.py --limit-gb 24 --interval-seconds 0.1 --report-interval-seconds 120 -- env PYTHONPATH=src ./.venv/bin/pytest -q`.

## Scope

This is synthetic validation of the Python physical-affine transfer adapter. It is not a JIST `ReorientVolume` numerical comparison for differing origins/spacings, full plus execution, Java composite equivalence, or real-MRI validation. MaSeg-derived masks and the downstream boundary/ACE/GVF/TGDM/surface/thickness stages remain to be composed from the raw-input workflow. No human MRI was used.

See [matching full-suite evidence](full-pytest-2026-10-05-native-dura-transfer.md).
