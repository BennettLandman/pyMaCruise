# Original MaCRUISE reference-run packet

This packet is prepared for someone who has access to the original MATLAB,
MIPAV, and JIST runtime. It contains the exact `.nii.gz` inputs used for the
current human BASE and PLUS Python runs, with SHA-256 hashes in
[`SHA256SUMS`](SHA256SUMS). It does not contain source code or runtime images.

## Input sets

These are two separate input sets. Do not mix files across them.

| Run | Input directory | Files | Grid and orientation |
| --- | --- | --- | --- |
| BASE | [`inputs/base-root`](inputs/base-root) | `T1.nii.gz`, `brain_mask.nii.gz`, `MaSeg.nii.gz` | All `(184, 256, 256)`, RAS; mask and MaSeg match the T1 grid. |
| PLUS | [`inputs/100260`](inputs/100260) | `T1.nii.gz`, `FLAIR.nii.gz`, `brain_mask.nii.gz`, `MaSeg.nii.gz` | T1, mask, and MaSeg `(176, 240, 256)`, RAS; FLAIR `(192, 512, 512)`, RAS, to be rigidly registered to T1. |

The BASE `T1.nii.gz` and the PLUS `100260/T1.nii.gz` are different scans and
have different dimensions. The file hashes are the definitive identity check.

## Source runs requested

Use MaCRUISE **v3 or newer**; v3.2.0 is preferred when available. Versions
before v3 had known bugs and are not an acceptable reference. Preserve the
complete per-stage output tree, logs, and parameters. Use separate output
directories under [`outputs`](outputs) so the BASE and PLUS products stay
distinct.

For both runs, use the supplied brain mask directly, retain the input voxel
grid, and only apply axis permutations/flips needed to report RAS. Do not run
BET or SPECTRE as a replacement for these supplied masks. Do not resample an
oblique image in the no-regrid profile. The included NIfTIs are RAS already.
The reference outputs used in comparison must be saved as `.nii.gz`; keep
source meshes in their original VTK form as well.

The checked-in legacy MATLAB entry points do not expose this exact contract:
`run_maCRUISE_v3.m` defaults to BET, and `run_maCRUISEplus.m` calls its legacy
skull-strip path. Please use a small run-only driver around the original
algorithm functions to supply the mask and bypass skull stripping. Do not
modify the original algorithm methods. The entry functions and settings below
identify the selected source workflow. FreeSurfer/MNI surface export, figures,
and viewer output are outside this run request; they are not needed for the
scientific image and mesh comparisons.

The source trees already provided to the project are:

- MaCRUISE MATLAB: `/Users/landmaba/dev/deleteme/masimatlab/trunk/xnatspiders/matlab/MaCRUISE/`
- MIPAV: `/Users/landmaba/dev/mipav/mipav-src/`
- JIST: `/Users/landmaba/dev/mipav/jist-scm-2026-09-26/`
- VandyCRUISE algorithms: `/Users/landmaba/dev/mipav/toads-cruise-scm-2026-09-26/VandyCRUISE/`
- MaCRUISE-specific JIST adapters and `MaGenericTGDM`: `/Users/landmaba/dev/mipav/toads-cruise-scm-2026-09-26/toads-cruise/`

Record the actual revision or source-tree hashes used by the runner. The
paths above describe the current machine; the staged inputs and atlas hashes
make the request portable if the runner uses equivalent source checkouts.

## BASE run driver

The packet includes run-only MATLAB drivers for BASE and PLUS and a NiBabel
mask helper at [`matlab/apply_brain_mask_nifti.py`](matlab/apply_brain_mask_nifti.py).
They follow the selected V3 stage order but call source stages directly to
accept the user mask without BET or SPECTRE. BASE applies the mask after
unmasked source N4, then runs both TOADS passes, Dura removal, the original
crop/boundary/tissue-fusion/ACE/GVF/TGDM stages, and refinement. PLUS calls the
original `FlairRegistration`/NiftyReg rigid wrapper, then source N4, applies
the supplied mask, and runs AtlasEM, Dura, LesionTOADS, tissue fusion,
ACE/GVF/TGDM, and refinement. This PLUS source reference uses Aladin while the
Python workflow uses ANTsPy; compare and report the resulting registration
transform separately, and do not interpret the end-to-end comparison as
registration-backend equivalence. Both drivers omit MNI export, FreeSurfer
conversion, figures, viewer files, and cleanup/deletion. Timestamped output
directories are created under `outputs/`; these are Git-ignored, and each
driver refuses to reuse an existing run directory.

From the repository root, with MATLAB, the selected MIPAV/JIST install, MASI
fusion checkout, ANTs N4 executable, and Python with `nibabel` available, run:

```matlab
addpath('source-reference-run/matlab');
run_base_selected_mask_reference( ...
  '/path/to/MaCRUISE', ...
  '/path/to/MIPAV-install', ...
  '/path/to/full-multi-atlas', ...
  '/path/to/ANTs/bin/N4BiasFieldCorrection', ...
  'python3');
```

The Python mask helper has two committed synthetic tests covering mask values,
RAS geometry, qform/sform preservation, and mismatched-grid rejection. The MATLAB
driver has not been executed in this environment because MATLAB and the selected
MIPAV/JIST application runtime are unavailable. Before using its
outputs as evidence, confirm the logged effective options, source hashes,
TOADS stop iterations, generated outputs, and mask/grid identity. Its stages
are the original MATLAB/JIST functions; the driver itself is new orchestration
that should be treated as an external-run recipe until someone completes and
returns this run.

### PLUS run driver

[`matlab/run_plus_selected_mask_reference.m`](matlab/run_plus_selected_mask_reference.m)
implements the corresponding selected-mask PLUS orchestration. It retains the
original `run_maCRUISEplus.m` rigid FLAIR registration and the
`maCRUISEpluspipeline.m` stage order, but replaces SPECTRE-derived masking with
the supplied T1-grid mask. It additionally applies the Dura mask to registered
FLAIR before LesionTOADS. From the repository root, call it with the same five
arguments as the BASE driver:

```matlab
addpath('source-reference-run/matlab');
run_plus_selected_mask_reference( ...
  '/path/to/MaCRUISE', ...
  '/path/to/MIPAV-install', ...
  '/path/to/full-multi-atlas', ...
  '/path/to/ANTs/bin/N4BiasFieldCorrection', ...
  'python3');
```

The driver and its stage order have only static source-review checks here; it
has not been parsed or executed by MATLAB. Original Aladin and the selected
MIPAV/JIST application must be present. Its selected-mask path is a run-only
adaptation, not an unchanged invocation of the historical wrapper. Return all
outputs and logs even if a stage fails so any driver mismatch can be corrected.

### BASE

- Inputs: `inputs/base-root/{T1.nii.gz,brain_mask.nii.gz,MaSeg.nii.gz}`.
- Atlas descriptor: `../data/atlases/cruise-atlas-10obj-toads2008.txt` from
  the repository root; SHA-256 is recorded in `SHA256SUMS`.
- Selected source entry: the v3 BASE sequence in `maCRUISEpipeline.m`, followed
  by `run_SegmentationRefine.m`. Retain outputs at every stage, including N4,
  both TOADS passes, Dura removal, TGDM inputs/level sets, and refinement.
- Match the current BASE run: no regridding; unmasked ANTs N4 with intensity
  rescaling off; V3 Dura regularization WM `0.05` / GM `0.02`; selected 10-class
  TOADS atlas; source-faithful segmentation refinement. The Python run's
  exact resolved options and command are in
  [`in-vivo-base-v3-no-mask-n4-source-profile-2026-10-14.md`](../docs/reports/in-vivo-base-v3-no-mask-n4-source-profile-2026-10-14.md).

### PLUS

- Inputs: `inputs/100260/{T1.nii.gz,FLAIR.nii.gz,brain_mask.nii.gz,MaSeg.nii.gz}`.
- Atlas descriptor and topology: `../data/atlases/lesionToads-atlas-2012/`
  (the 12-object LesionTOADS descriptor and its matching topology NIfTI); hashes
  are in `SHA256SUMS`.
- Selected source entry: `FlairRegistration` → `maCRUISEpluspipeline.m` →
  `run_SegmentationRefine.m`. Retain all nine workflow stages and their
  intermediates, especially the FLAIR-to-T1 transform, N4 products, TOADS/Dura
  handoffs, ACE/GVF/TGDM fields, refined segmentation, and surfaces/thickness.
- Use the axis-only `--no-regrid-inputs` geometry profile to match the current
  historical-environment comparison. Record the rigid FLAIR-to-T1 transform,
  interpolation, N4 options, lesion convergence settings, and Dura transfer
  behavior actually used. The Python reference run manifest and command are
  detailed in
  [`in-vivo-100260-plus-full-run-2026-10-14.md`](../docs/reports/in-vivo-100260-plus-full-run-2026-10-14.md).
- Record the exact compiled v3+ runtime version. The supplied VandyCRUISE
  source passes plugin `firstLimit=3.0` through `AlgorithmLesionToads` to
  `LesionToadSegmentation`, where it becomes the effective limit. Include
  iteration-level thinning event counts if available. The earlier zero-limit
  continuation was a sensitivity experiment, not the source wrapper profile;
  see [`p41z-source-wrapper-limit0-continuation-2026-10-15.md`](../docs/reports/p41z-source-wrapper-limit0-continuation-2026-10-15.md).

If the source driver cannot accept the supplied mask without invoking BET or
SPECTRE, return that run as an **unmatched historical-default reference** and
record the exact mask produced by the wrapper. It can inform diagnosis but
cannot close the matched-source acceptance row. Please do not silently use a
different mask, atlas, orientation, or registration transform.

## Files and metadata to return

For each run, return the complete output directory (a compressed archive is
fine), including:

1. Every intermediate scientific NIfTI and final segmentation, surface, and
   thickness output. Preserve RAS geometry and use `.nii.gz` for volume images.
2. The full console, MATLAB, MIPAV/JIST, ANTs, and scheduler logs; do not strip
   error messages or failed-stage logs.
3. A `RUN_INFO.md` or text file containing the exact source revision/tag,
   MIPAV/JIST/plugin versions, MATLAB/MCR version, Java version, ANTs version,
   operating system, launch command, effective options, atlas paths, and
   whether the run used the supplied mask or a generated mask.
4. SHA-256 hashes for all four inputs (three for BASE), atlas descriptor and
   referenced atlas files, and returned outputs. Include image dimensions,
   voxel sizes, qform/sform, axis codes, and the FLAIR registration transform
   for PLUS.
5. The stopping iteration, energy ratio, stop reason, and whether stopping was
   natural convergence or an iteration limit.

Do not include FreeSurfer surface conversion, MNI export, figures, or viewer
files in the archive. Retaining them is unnecessary for this port comparison.

After the output archive is available here, I can map its intermediates to the
Python manifests and run the voxel, Dice, transform, and surface comparisons.
The original runtime remains the only item this host cannot execute itself.

Before comparison, audit the returned directory with the repository tool:

```sh
PYTHONPATH=src .venv/bin/python scripts/audit_source_reference_run.py \
  /path/to/returned/base-run --require-run-info --require-ras
```

The audit writes `source-run-manifest.json` with SHA-256 for every regular file,
and NIfTI shape, affine, qform/sform, voxel size, orientation, and finite-value
checks. `--require-ras` fails on non-RAS image outputs; omit it only when
deliberately collecting a legacy orientation for later conversion/comparison.
It rejects uncompressed `.nii` images and symlinks but never modifies source
outputs. Mesh and log files are hashed without interpretation.
