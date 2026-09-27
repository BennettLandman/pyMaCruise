# Image and Coordinate Conventions

**Status:** P06 started; source review recorded, synthetic landmark validation pending.
**Scope:** Establish the spatial rules PyMaCRUISE will use when loading, transforming, cropping, combining, and writing NIfTI volumes. All persisted image files use compressed NIfTI (`.nii.gz`); this is a working contract until the synthetic checks in P06 pass.

## Source behavior observed

| Topic | Evidence from the selected MaCRUISE workflow | Port implication |
|---|---|---|
| Input arrays and headers | `prepare_raw.m` copies/compresses T1 and MaSeg files without applying a reorientation. The pipeline's `load_untouch_nii_gz` / `save_untouch_nii_gz` helpers are used for direct array edits. | This historical behavior does not set the PyMaCRUISE output convention. On ingest, axis-reorient non-RAS images to RAS and emit a visible warning identifying the input orientation and conversion. |
| Toads-to-reference transfer | `ToadsCruiseSpaceToOriginalSpace` calls `ReorientVolume.solve(target, reference, Resolution.UNCHANGED, NEAREST_NEIGHBOR)` and then assigns the reference header. The companion `ReorientVolume` implementation maps target axis orientation to reference orientation; this is an orientation operation, not registration to align unrelated anatomy. | Keep orientation conversion explicit and separate from registration. Use nearest-neighbor for labels/masks in this stage. The exact MIPAV handling of oblique directions and header fields still needs an oracle/synthetic comparison. |
| Reference-space output helper | `run_TransferToOriginSpace_matlab.m` copies the target voxel array into the reference NIfTI structure and saves it; its commented-out flips are inactive. | This MATLAB helper assumes corresponding voxel indexing and does not itself resample. The Python port must validate shape/index compatibility before reproducing this behavior rather than silently assigning a new header. |
| Cropping | `CropBasedOnGM` computes crop parameters from GM and applies the same parameters to the companion image. The configured border is 5 voxels. | Store crop start/stop indices and update the affine origin using the crop start; apply one crop transform consistently to aligned companions. |
| Physical coordinate system | NIfTI files supply index-to-world metadata; ANTs uses ITK conventions internally, while NiBabel exposes its selected affine in RAS+ world coordinates. | All reported world coordinates and persisted volume/surface outputs use RAS+ millimetres. Convert explicitly at library boundaries that use LPS. Never infer a left/right flip from array shape alone. |

## Working Python contract

1. **RAS output orientation:** All output images use voxel axes ordered Right-Anterior-Superior (RAS); surfaces and reported coordinates use RAS+ millimetres. For an input whose axis codes are not RAS, permute/flip the voxel axes while preserving the physical mapping, emit a user-visible warning, and record the original orientation and applied transform in run metadata. This normalization applies to images, label maps, and surface geometry before downstream use.
2. **Index origin and voxel centers:** Integer indices refer to voxel centers. The homogeneous affine maps `[i, j, k, 1]` to the world coordinate of that voxel center; no half-voxel offset is added by PyMaCRUISE.
3. **World frame and units:** Public geometry helpers and every output use RAS+ coordinates in millimetres. If a library boundary uses LPS+, convert with `diag(-1, -1, 1, 1)` and test round trips.
4. **Input normalization and warnings:** Inspect orientation/frame metadata on every spatial input. If an image, surface, transform, or other spatial input is not in RAS, normalize/reorient it to RAS before use and display a warning naming the source orientation/frame and conversion. Preserve the transform in the run manifest. For inputs with no trustworthy coordinate-frame metadata, require an explicit frame instead of guessing.
5. **Affine selection:** For NIfTI input, use NiBabel's selected image affine (sform when coded, otherwise qform when coded, otherwise its documented fallback). Record both raw qform/sform and codes when validating; if both are coded but materially disagree, flag the image for explicit handling rather than silently choosing one in pipeline stages.
6. **Interpolation:** Use nearest-neighbor for discrete labels and masks. For continuous-valued MRI/level-set data, choose the interpolation method per source stage and document any ANTs/MIPAV boundary conversion. No universal interpolation default is assumed.
7. **Headers:** Persist NIfTI outputs in RAS orientation. On unchanged-grid voxelwise operations, preserve the normalized reference affine and relevant metadata while setting an output datatype compatible with the resulting values. On crop, update affine translation for the crop start. On resampling, construct output geometry from an RAS reference grid. Do not copy a reference header onto data with a different shape or voxel-to-world mapping without an explicit validated rule.
8. **Geometry checks:** Before voxelwise operations, require matching shape and affine within a declared tolerance, or perform a separately requested resampling step. Matching shape alone does not imply alignment.

## Pending decisions and validation

- Compare RAS normalization and warnings to NiBabel orientation handling and MIPAV `ReorientVolume` on axis-aligned, permuted-axis, flipped-axis, oblique, and LPS-coordinate synthetic images/surfaces.
- Confirm exact crop bounds, edge clipping/padding, and affine update against `CubicVolumeCropper` using a small asymmetric GM mask and landmark companion image.
- Inspect NiBabel and ANTsPy runtime behavior after the P12 environment is created; verify transform direction and physical units with a hand-calculated translation/rotation phantom.
- Add synthetic files and automated checks under `tests/data/synthetic/` and `tests/` at the planned fixture/test steps. This document does not claim those checks have passed.
- Resolve conflicting coded qform/sform behavior against the legacy runtime before treating the affine-selection policy as a final compatibility guarantee.

## Source references

- MaCRUISE `pipeline/prepare_raw.m`
- MaCRUISE `pipeline/run_TransferToOriginSpace_matlab.m`
- MaCRUISE `pipeline/run_TransferToOriginSpace.m` and `run_TransferToOriginSpace_jist.m`
- Companion `toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/utilities/ToadsCruiseSpaceToOriginalSpace.java`
- Companion/Vandy `edu/jhu/ece/iacl/algorithms/volume/ReorientVolume.java`
- Companion `toads-cruise/src/edu/vanderbilt/masi/plugins/CRUISE/utilities/CropBasedOnGM.java`
- MaCRUISE crop and output-space mappings in [source-map.md](source-map.md)
