# In-vivo base crop identity-formula comparison — 2026-10-11

## Result

The corrected V3 preparation on the T1 grid matches the pinned
`PlusLesionTransferCropRasIdentityReferenceHarness` formula exactly. The
full-grid inputs are RAS `208 × 309 × 305`; with the source-selected 5-voxel
edge, Java and saved preparation provenance both give inclusive bounds
`[[30,170],[58,230],[105,235]]` and crop size `141 × 173 × 131`.

The Java formula output matched the saved Python crop arrays exactly for
atlas GM, atlas WM, atlas mask, atlas brain-tissue mask, TOADS GM, TOADS WM
fill, and TOADS mask. It also matched a direct slice of the MaSeg labels used
as the eighth input. Each array has **3,195,483 voxels**, with zero changed
voxels and zero maximum absolute error. Each saved TOADS crop also equals the
corresponding slice from its already-transferred full-grid map: all three
checks had zero changed voxels.

The Java harness ran with `-Xmx768m`; peak sampled RSS was **0.667 GiB** under
the 2 GiB watchdog.

## Reproduction

From the `pymacruise` repository root:

```sh
PYTHONPATH=src .venv/bin/python scripts/compare_base_preparation_crop_identity_java.py \
  --preparation-dir build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11 \
  --output-dir build/invivo-base-2026-10-10-ras-grid/java-crop-formula-v3-t1-grid-2026-10-11
```

The runner checks the pinned harness SHA-256
`a3fb7afa0365609e513d0cb480e169643637144e3010909bd0409befec266c94`, packs
the eight full-grid arrays, executes the harness, compares all returned
arrays, and writes `comparison.json` and temporary arrays under the requested
ignored `build/` directory. The inputs include all three transferred TOADS
maps; no upstream preparation or CRUISE stage is rerun.

## Scope

This is an **identity-grid formula comparison**, not execution of JIST
`CropBasedOnGM` or `CubicVolumeCropper`. It supports exact crop-bound and
voxel-copy agreement for these held-fixed RAS arrays. The TOADS crop checks
verify copying from already-transferred same-grid maps; they do not validate
the earlier oblique-grid resampling against MIPAV `ReorientVolume`. Cortical
separation, that transfer, and full-workflow equivalence remain unvalidated by
this report.
