# Workflow Input and Output Contract

**Status:** P07 validated on synthetic inputs; downstream pipeline behavior remains to be implemented.
**Reference:** MaCRUISE SVN revision 9360, `run_maCRUISE_v3.m` and `pipeline/maCRUISEpipeline.m`; MaCRUISEplus uses `run_maCRUISEplus_v1.m` and `pipeline/maCRUISEpluspipeline.m`.

## Supported workflow profiles

| Profile | Required input volumes | Additional behavior |
|---|---|---|
| Base V3 | T1-weighted MRI, MaSeg multi-atlas segmentation, and a user-supplied brain mask | N4 preprocessing, TOADS, V3 topology correction and segmentation refinement. The brain mask is used instead of BET/SPECTRE skull stripping. FLAIR is not an input to this runner. |
| Plus V1 | T1-weighted MRI, brain mask, FLAIR, and MaSeg multi-atlas segmentation | Rigid FLAIR-to-T1 registration followed by lesion-TOADS and the full MaSeg-informed plus surface workflow. The brain mask is supplied by the user; BET/SPECTRE skull stripping is not run. |

The two profiles must be explicit configuration choices. Do not infer plus mode from the presence of an arbitrary file, nor silently use the legacy no-MaSeg or Origin variants.

## Input requirements

### T1 and FLAIR

- Persisted workflow inputs, intermediates, and outputs use 3D, single-volume compressed NIfTI (`.nii.gz`) images with finite, scalar voxel data and a valid spatial affine. If an external image is supplied in another format, normalize it to `.nii.gz` before storing or reusing it in the workflow.
- T1 is required for both profiles. FLAIR is required only for the plus profile and must be registered to the T1 reference grid by the explicit rigid-registration stage.
- A brain mask is required for both profiles. Supply a binary 3D `.nii.gz` mask whose shape and physical affine exactly match T1 after RAS orientation normalization; accepted foreground encodings are 1 or 255, background is 0. The mask must be nonempty. It is a user-provided input, not an internally generated BET or SPECTRE result. Both restartable full-workflow runners regrid oblique inputs by default to an orthogonal RAS grid; `--no-regrid-inputs` disables that step and requires an already axis-aligned RAS grid for TOADS.
- Preserve the stored voxel order and header on read. Keep physical geometry explicit. Input grids are validated independently; operations that combine arrays require matching grids or a documented resampling step. The full runners' default regrid uses linear interpolation for T1 and nearest-neighbor interpolation for the brain mask and MaSeg, all on one T1-derived grid at the source T1 voxel spacing. It reports the action and records the derived files and hashes in the run manifest. `--no-regrid-inputs` preserves an already axis-aligned RAS grid without interpolation; oblique TOADS inputs are rejected in that mode. FLAIR remains on its source grid and is registered to the regridded T1 grid by the plus preprocessing stage.
- The baseline T1 acquisition default is `T1_MPRAGE`; acquisition metadata is not inferred from the filename. N4 remains the selected bias-correction stage. BET and SPECTRE are excluded; no FSL or FreeSurfer runtime is needed.

### MaSeg

- MaSeg is a 3D, single-volume NIfTI label map. Values are compared by exact equality in the MATLAB source; the port should reject non-finite or non-integer-valued labels rather than silently round them.
- The source does not remap the MaSeg labels. It derives binary `uint8` masks whose selected voxels are 255 and all others are 0.
- Active cortical-separation sets in `cortical_seperation.m` are:
  - **Atlas_GM:** `100–109`, `112–125`, `128–129`, `132–157`, `160–187`, and `190–207`; labels `188` and `189` are excluded by the source's explicit list.
  - **Atlas_WM and Atlas_Mask:** `4, 23, 30, 36, 37, 44, 45, 51, 52, 55–62` (the listed integer labels; `Atlas_Mask` intentionally uses this WM set alone).
  - **Atlas_BrainTissueMask:** GM labels plus WM labels plus `31, 32, 47, 48, 49, 50, 75, 76`.
- For base and plus segmentation, MaSeg is the spatial reference for cortical-separation products and for transferring TOADS results back to the segmentation space. Geometry differences must be resolved explicitly; same array shape alone is not sufficient evidence of alignment.

### `pymacruise plus-run` fusion mask

The standalone `plus-run` continuation accepts `--multi-atlas-mask` as its
fusion support mask. Supply the MaSeg-derived `Atlas_BrainTissueMask.nii.gz`
(GM labels, WM labels, and the additional brain-tissue labels listed above).
Do not pass `Atlas_Mask.nii.gz` to this option: that separate output contains
only the WM label set and is not the fusion support mask. The CLI flag name is
retained for compatibility with existing invocations.

## Output contract

Outputs are written beneath a caller-selected output directory. Every volume and surface output is reported in RAS (RAS-oriented NIfTI voxel axes; RAS+ physical coordinates in millimetres for surface geometry). Non-RAS input images/surfaces are reoriented before use with a visible warning, and the source frame/orientation plus conversion are recorded. The Python package should provide a deterministic manifest with profile, input paths and hashes, software/dependency versions, configuration, asset paths, coordinate frames/transforms, and output paths.

### Required scientific products

- Refined segmentation `SEG_refine.nii.gz` on the MaSeg reference grid, retaining the input MaSeg's label values where the MATLAB refinement logic preserves them.
- Cortical-separation masks `Atlas_GM.nii.gz`, `Atlas_WM.nii.gz`, `Atlas_Mask.nii.gz`, and `Atlas_BrainTissueMask.nii.gz`.
- TOADS or lesion-TOADS tissue products and their scientifically consumed intermediates.
- Final inner, central, and outer CRUISE level sets and surfaces, plus required inverse-cropped/origin-space products.
- Thickness volume and summary products generated from the final inner/outer pair.
- For full legacy-output mode, MNI-space level sets/surfaces. These require the exact MNI reference and mask, which are currently absent from the supplied source trees. Until obtained and validated, expose this output as unavailable/disabled with a clear reason; do not synthesize a replacement template.

### Excluded legacy presentation products

- Do not port FreeView batch scripts, screenshots/figures, MATLAB figures, rendering, or viewer behavior.
- Do not create FreeSurfer/FreeView-specific ASC presentation exports or subject-specific FreeSurfer comparisons as part of this port.
- Installing or using PyMaCRUISE must not require a FreeSurfer installation, license, environment configuration, dataset, or generated output.
- Preserve scientific surface meshes, level sets, coordinate transforms, and thickness products required by the segmentation workflow. The later P47 discussion will consider modern Python visualization tools as a separate follow-on effort.

## Failure behavior and provenance

- Fail early with the input path and reason for missing files, invalid dimensions/data, invalid label values, unusable geometry, incompatible grids, or missing required bundled atlas/configuration files.
- Never treat an existing stale output as valid solely because its filename exists. Use input/configuration/asset signatures in the run manifest before resuming stages.
- Do not write partial output as if the full profile succeeded. Record stage state and distinguish optional omitted products from failed required stages.
- Human subject input images are not copied into tests or committed. Synthetic fixtures are generated and checked in during the planned fixture steps.
- Every emitted volume is RAS-oriented; every mesh and coordinate report is in RAS+ millimetres. A conversion warning is emitted whenever an input is found in another orientation/frame.

## Pending synthetic acceptance cases

- Valid base profile; valid plus profile; absent FLAIR accepted only for base and rejected for plus.
- Missing T1/MaSeg, 2D/4D volumes, non-finite intensities, malformed affine, non-integer MaSeg labels, and out-of-contract label values produce clear diagnostics.
- Missing, empty, nonbinary, shape-mismatched, or affine-mismatched brain masks produce clear diagnostics; valid non-RAS masks warn and normalize to the T1 RAS grid.
- MaSeg labels map to the exact expected synthetic masks, including labels near each set boundary and labels intentionally excluded.
- Shape-equal but affine-shifted inputs are rejected or explicitly resampled; FLAIR output is on the T1 reference grid after rigid registration.
- Permuted/flipped-axis and LPS image/surface inputs are normalized to RAS with visible warnings and recorded transforms; generated volume and mesh outputs pass RAS orientation/frame checks.
- Output manifest identifies required scientific outputs and MNI outputs unavailable due to missing MNI assets; it contains no visualization-product requirement.

The MaSeg-to-mask label sets are implemented in `src/pymacruise/cortical_separation.py` and checked against the exact source lists, including labels at range boundaries and exclusions. The remaining cases above are requirements, not claims of completed implementation or validation. They will be assigned to the deterministic fixture and test roadmap steps.
