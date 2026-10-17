# AtlasEM compressed NIfTI payload loader

## Scope

Commit `dec82b5f0fc449de3fc83db86a604201e7303975` adds `load_atlas_em_nifti_payloads`. It reads a 3D integer topology map and each listed 3D shape prior from `.nii.gz`, normalizes all maps to RAS with the shared warning behavior, checks descriptor dimensions, rejects non-finite values and inconsistent grids, and reports descriptor-versus-header voxel-size differences while retaining descriptor resolutions for algorithm geometry. Shape arrays are float32; topology labels retain their integer representation.

The source Java implementation reads raw signed-byte topology and little-endian float32 shape files. This loader deliberately targets the repository's `.nii.gz` storage contract. Synthetic NIfTI loading and axis reorientation are validated, but no same-data Java raw-file comparison is claimed. Header affine consistency is used to ensure all converted maps occupy the same grid; actual `ibsr8obj` assets are still missing.

## Validation

Command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 .venv/bin/pytest -q tests/test_atlas_em_descriptor.py
```

The descriptor-focused suite passed **10 tests in 0.64 s**. The synthetic cases verify axis-flip voxel preservation with the expected RAS warning, topology and shape arrays, rejection of mismatched grids and invalid topology values, alongside the descriptor metadata checks. A broader AtlasEM regression across descriptor loading, statistics, memberships, preprocessing, crop/recrop, registration-pyramid preparation, output mapping, CLI handoff, and synthetic Dura integration passed **33 tests in 6.36 s** on tested commit `cfd9929762d6786dbd60f42248c20d60320cc166`, with three expected orientation warnings. The payload loader implementation is commit `dec82b5f0fc449de3fc83db86a604201e7303975`. The full suite and Java raw-file payload comparison were not run.
