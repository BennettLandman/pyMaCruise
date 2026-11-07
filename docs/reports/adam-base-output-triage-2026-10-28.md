# Adam's BASE output comparison, October 28, 2026

## Received-file audit

The follow-up supplied `/Users/landmaba/Downloads/PyMaCRUISE/`, containing
Python inputs/outputs and an independent MaCRUISE v3.2.0 Output tree. The
findings below supersede the screenshot-only interpretation retained later
in this report. The original files were read without modification.

All **68 Python stage-output hashes match** the received manifest after
relocating its Linux paths. All three input hashes match. Original `Raw/Seg`
preserves the supplied MaSeg values and affine; original `Raw/T1` preserves
the supplied T1 values after float32 conversion and its affine. This supports
scan identity, not equality of preprocessing or algorithms.

### Numerical settings and completion

The supplied minimal command used the smoke-test defaults. `run.json` records
`status=complete`, `converged=false`, and
`stage-run-complete-but-toads-did-not-meet-energy-threshold`. Both TOADS passes
hit their two-iteration limit; neither reached the energy threshold. ACE
diffusion used one iteration, and TGDM initial/inner/central each used one.
Registration levels and iteration caps were also one. WM iso-level was 0.6
and updated-CSF mode was 0.

This is a confirmed configuration problem for an MRI agreement evaluation.
Its contribution to the observed quality deficit must be measured with a
new run; increasing limits is not proof of original/Python equivalence.

### SEG_refine datatype and grid

| Product | Dimensions | Stored datatype | Fractional voxels |
| --- | --- | --- | --- |
| Supplied MaSeg | 256×256×150 | float64 | 0 |
| Original SEG_refine | 256×256×150 | float64 | 0 |
| Python working MaSeg | 257×257×151 | float64 | 0 |
| Python SEG_refine | 257×257×151 | int32 | 0 |

All four use identity scaling (slope 1, intercept 0). The reported fractional
Python SEG_refine is therefore not present in the supplied file. A viewer's
interpolation setting or selection of another product remains to be checked.
There is no evidence here for rounding or recasting this Python output.

Python's default regridding created an axis-aligned working T1/MaSeg grid,
explaining the extra voxel on each axis. Original SEG_refine retains the
input's slightly oblique affine. These grids need a defined spatial mapping
for comparison, not copied headers.

### Surface-coordinate verification

Original `Surfaces/target_image_GMimg_centralSurf.vtk` stores zero-based
**voxel coordinates**. Re-extracting the supplied original central field in
voxel coordinates produces the same vertex/triangle counts; nearest-vertex
differences are at most **0.0000203 voxel**, consistent with VTK rounding.
Applying original `Raw/Seg`'s affine to those VTK vertices reproduces the
supplied `Surfaces_FreeView` ASC vertices with a maximum difference of
**0.00000863 mm**. This verifies the mapping against a supplied independent
export, rather than inferring it from mesh bounds.

Python's VTK is explicitly RAS+ millimetres. The original in-plane spacing
is approximately **0.890076 mm**; displaying original voxel numbers as
millimetres exaggerates its in-plane size by roughly 12.35% and also omits
its translation/rotation. The existing surface extractor can reconstruct the
original field closely when held fixed; that does not validate Python's
upstream field generation.

### Disagreement after coordinate normalization

Mapping the original mesh through its verified affine leaves a symmetric
nearest-vertex mean difference of **1.57647 mm**, with maximum **11.4088 mm**.
These are vertex-nearest diagnostics, not vertex-to-triangle distances or a
continuous Hausdorff bound. No best-fit registration was applied.

The original central mesh has one connected component; Python has 140.
Python's largest component contains 304,984 of 307,656 vertices (99.13%);
the other 139 components contain 2,672 vertices. Both meshes have zero open
and nonmanifold edges. Closed triangles do not establish matching topology
or anatomical quality.

After nearest-neighbor physical resampling of Python SEG_refine to original
SEG_refine's grid, **585,542 voxels differ**, or **32.22%** of the foreground
union; **102/134** non-background labels have Dice below 0.97. Foreground Dice
is **0.84948**. These are diagnostic metrics for different numerical profiles,
including label resampling; they are not a controlled port-acceptance result.
The file-based common-grid montage shows more peripheral clearing in the
received Python labels. This replaces the tentative screenshot-only reading
of which result filled more sulci.

Machine-readable results and reproducible analysis:

- [Received-output audit](../../build/adam-output-triage-2026-10-28/received-output-audit.json)
- [Audit script](../../build/adam-output-triage-2026-10-28/audit_received.py)
- [Common-grid montage](../../build/adam-output-triage-2026-10-28/received-segmentation-common-grid.png)
- [Original central mesh in verified RAS mm](../../build/adam-output-triage-2026-10-28/original_central_RAS_mm_DIAGNOSTIC.vtk)
- [Python central mesh in RAS mm](../../build/adam-output-triage-2026-10-28/python_central_RAS_mm_DIAGNOSTIC.vtk)

The derived mesh/label exports preserve existing output shapes or explicitly
resample labels for comparison. They do not repair the reconstruction.

### Prepared next run

Update after the owner's defaults/grid correction: version 0.1.1 now defaults
to the full numerical profile and source-compatible no-regrid working geometry,
then exports final products back to the original input grid. The linked rerun
script now requires the native-export module and selects this corrected profile.
See [the implementation and validation report](base-full-defaults-native-output-2026-10-28.md).

[The rerun script](../../build/adam-output-triage-2026-10-28/run-base-mri-profile.sh)
uses the existing installed atlas, the same three relative input paths,
explicit MRI limits from the selected BASE profile, and a new
`outputs_mri_profile` directory. Run it from Adam's test directory with his
updated environment active. It uses `--no-regrid-inputs` for source-compatible
processing and restores native geometry in the final export.
Do not resume the old `outputs/` run with changed settings. Repeating the new
script with `--resume` is valid for the same new run.

Shell syntax and CLI argument parsing passed. **The new human MRI run has not
been executed.** Existing differences in supplied-mask/skull-stripping
preprocessing and full original/Python equivalence remain open. No runtime
algorithms, defaults, received image headers or received outputs were changed.

## Initial triage before the output folder was supplied

Adam reports a completed run with differences in SEG_refine, image dimensions
and affines, apparent noninteger segmentation output, and central-surface
shape/alignment. The supplied screenshots show apparently more filled sulci
in Python's segmentation and a visibly rougher Python central surface.
Different slice positions, camera settings, display scaling and coordinate
frames limit conclusions from screenshots. No header tables or underlying
output files from this run have been supplied. This report does not diagnose
the cause or claim a correction.

## First check: resolved BASE settings

`BaseRunOptions` and the BASE CLI default to short synthetic smoke-test caps:
TOADS 2, ACE diffusion 1, and TGDM 1 iteration for each surface unless explicit
surface counts override the compatibility count. Registration defaults are
also short. Completing these stages does not establish convergence or
agreement with original MaCRUISE.

The explicit MRI comparison recipe in
[`python-base-checks.md`](../testing/python-base-checks.md) instead selects
TOADS 50, ACE diffusion 25, TGDM initial/inner/central 6/4/7, registration
levels/iterations 4/50, post-refresh/max registration counts 5/5, Dura 5,
WM iso-level 0.5, and updated-CSF mode 2. Its no-regrid geometry profile is
specific to the supplied comparison case and discards oblique directions;
it is not a general physical-space alignment solution.

The local October 27 Python-only bundle contains `docs/tools-reference.md`
and PLUS instructions, but omits `docs/testing/python-base-checks.md`. Its
BASE note tells testers to retain their chosen settings without giving the
full MRI recipe. A minimal BASE command could therefore use smoke-test
settings. Adam's actual command and `run.json` are needed to determine whether
this happened; appearance alone is insufficient.

## Distributed-wheel refinement check

The tested wheel is
`build/python-only-tester-bundle/wheels/pymacruise-0.1.0-py3-none-any.whl`,
SHA-256 `b76ad934c4aa9cbbd04aa8aa93a6e9f2f4c8bdde4bb7f105af10cd0368bd21bb`.
It was unpacked to an ignored build directory and imported directly from
that unpacked wheel with the existing development dependencies. This is a
focused wheel-code check, not a fresh installation or Adam's environment.

Using `tests/data/synthetic/refinement/` level sets and segmentation values:

| Input label storage | Written SEG_refine storage | Values | Geometry | Scaling |
| --- | --- | --- | --- | --- |
| uint16 | uint16 | All integers | Input shape/affine retained | Slope 1, intercept 0 |
| float32, integer values | int32 | All integers | Input shape/affine retained | Slope 1, intercept 0 |

Both products were written through the wheel's `save_nifti_ras` and reloaded
from disk. The wheel's refinement function rejects fractional input labels.
This did not reproduce Adam's noninteger output. Check the exact active
SEG_refine file, stored datatype, scaled values, slope/intercept, package
location and wheel identity before choosing an export fix. Viewer
interpolation alone cannot establish the stored datatype or fractional labels.

## Geometry and surface interpretation

BASE's default regrid path resamples oblique inputs to an axis-aligned RAS
grid and can increase image dimensions. Its no-regrid compatibility path
permutes/flips axes and discards oblique directions while preserving the
image center. The working T1 must accompany outputs for overlays. Compare
qform/sform, voxel sizes and working-input geometry before any resampling.
Do not copy headers between different grids to force apparent agreement.

Python meshes use full level-set affines to produce RAS+ millimetre vertices.
The original surface's actual frame/export stage must be identified before
comparing placement or size. Camera appearance is not a mesh-scale measurement.

Existing held-fixed Java/Python extraction comparisons concern identical
input level sets. They do not validate the upstream fields or Adam's run.
The rough surface and altered segmentation warrant comparing upstream
GM/WM maps, ACE and TGDM after verifying inputs, settings and geometry;
cosmetic smoothing would not establish agreement.

## Initially requested material (now supplied in the follow-up)

- Exact launch command, `run.json`, stage logs and installed package location.
- Both SEG_refine NIfTIs, both central meshes and corresponding central level
  sets, retaining original files and their directory identities.
- Input T1/MaSeg/mask identities and Python's generated working T1, if present.
- Original MaCRUISE version/command and the viewer/mesh import convention.

No runtime algorithms, defaults or image headers were changed during triage.
