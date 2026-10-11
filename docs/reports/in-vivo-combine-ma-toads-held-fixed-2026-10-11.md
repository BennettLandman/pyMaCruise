# In-vivo `CombineMAandTOADS` held-fixed comparison — 2026-10-11

## Result

The selected V3 Java `CombineMAandTOADS.execute()` output matches all three
saved Python fusion outputs bit-for-bit on the corrected human preparation.
Each comparison covers **3,195,483 voxels** on the shared RAS crop grid
`141 × 173 × 131`:

| Output | Changed voxels | Maximum absolute error | Mean absolute error |
|---|---:|---:|---:|
| GM | 0 | 0 | 0 |
| WM fill | 0 | 0 | 0 |
| Mask | 0 | 0 | 0 |

The Java harness completed with **0.245 GiB peak sampled RSS**, below its
2 GiB watchdog limit. It ran with the selected defaults: 3D center-plus-six-
neighbor WM erosion and `GMconfident=0.001`.

The held-fixed inputs and Python outputs are from
`build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11/`.
This is the MaSeg file recorded by the full run manifest, regridded by the
input-preparation step onto the supplied T1 and brain-mask grid
`208 × 309 × 305` with RAS axis codes. The second-pass TOADS maps also have
that same grid, so this human case exercises the identity-grid transfer.
The preparation provenance records `remove_cerebellum=false` and the V3 order:
MaSeg cortical separation, nearest-neighbor transfer to the MaSeg grid,
shared crop from positive MaSeg GM, then fusion. The checked inputs were the
cropped atlas GM, atlas WM, atlas brain-tissue mask, TOADS GM, TOADS WM fill,
and TOADS mask. This is **source comparison of fusion on held-fixed Python
preparation products**; it does not establish parity of the earlier stages or
of the whole workflow.

## Reproduction

From the `pymacruise` repository root, run:

```sh
PYTHONPATH=src .venv/bin/python scripts/compare_base_tissue_preparation_saved_java.py \
  --preparation-dir build/invivo-base-2026-10-10-ras-grid/source-selected-maseg-toads-preparation-v3-t1-grid-2026-10-11 \
  --output-dir build/invivo-base-2026-10-10-ras-grid/java-fusion-v3-t1-grid-2026-10-11
```

The command reuses the preserved source and headless harness:
[`CombineMAandTOADS.java`](../../tests/java_reference/upstream/edu/vanderbilt/masi/plugins/CRUISE/MaCRUISE/CombineMAandTOADS.java),
[`MASImorph.java`](../../tests/java_reference/upstream/edu/vanderbilt/masi/plugins/CRUISE/utilities/MASImorph.java),
and [`PlusTissueUpdateSourceHarness.java`](../../tests/java_reference/plus_tissue_update/PlusTissueUpdateSourceHarness.java).
It writes packed arrays, Java outputs, and a machine-readable comparison report
under the requested output directory. The human input hashes, preparation
provenance hash, source settings, output metrics, and watchdog result are
recorded in that generated `comparison.json`.

## Remaining evidence gap

This comparison validates fusion only on held-fixed Python-prepared inputs.
The human transfer was an identity-grid operation: MaSeg and all TOADS inputs
had matching shape, affine, and RAS orientation. The selected Java adapter is
orientation-code based and preserves extents; nonidentity axis permutations,
flips, and other grids still need a direct source comparison. MaSeg cortical
separation and `CropBasedOnGM` / `CubicVolumeCropper` bounds and voxel copying
also lack a held-fixed Java comparison for this human run. Full base-workflow
parity and final ribbon quality remain separate checks.
