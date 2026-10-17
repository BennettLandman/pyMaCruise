# PyMaCRUISE command-line input/output reference

This page is the expanded `--help` contract for the supported inputs and
outputs. The command line exposes preparation, individual scientific stages,
restartable base and plus workflows, workflow configuration/preflight, and
standalone algorithm commands. Synthetic full-workflow validation is available;
source equivalence on complete workflows and matched real-MRI inputs remains in progress. Anatomical quality audits are diagnostic, not independent port-completion gates; see [the validation policy](validation-tolerances.md#port-acceptance-objective--owner-clarification-2026-10-14).

## `pymacruise workflow-preflight`

```text
pymacruise workflow-preflight --config workflow.json [--overwrite]
```

This command is a **dry-run preflight**, not a segmentation workflow. It checks
the JSON configuration, compressed NIfTI inputs, profile-specific bundled
atlas descriptor and referenced NIfTI headers, core dependencies, and the
optional ANTsPy dependency needed by the selected N4 correction. Non-RAS
inputs are reoriented for validation with a warning; original and normalized
axis codes are recorded. Atlas NIfTI payloads are not loaded.

The command writes `run.json` in the configured output directory. A directory
without a manifest must be empty. If `run.json` already exists, `--overwrite`
replaces only a regular manifest file and leaves other files untouched;
symlinked manifests are rejected. The manifest
sets `dry_run` to `true`, records `workflow_executed: false`, and leaves
`stages.pipeline` as `not_started`. No N4 correction or segmentation stage is
run. If ANTsPy is missing, install it with
`python -m pip install 'pymacruise[ants]'` before retrying.

Configuration paths are resolved relative to the JSON file. Version 1 uses
this schema for the base profile:

```json
{
  "schema_version": 1,
  "profile": "base",
  "t1_path": "../inputs/T1.nii.gz",
  "maseg_path": "../inputs/MaSeg.nii.gz",
  "brain_mask_path": "../inputs/BrainMask.nii.gz",
  "output_dir": "run-output",
  "bias_correction": "N4"
}
```

The plus profile additionally requires `flair_path` and may provide an
`atlasem_descriptor_path` and a `plus_options` object using the field names in
`PlusPipelineOptions`. The optional AtlasEM descriptor enables the selected
AtlasEM/Dura stage before LesionTOADS cropping; every referenced payload must
exist as a same-grid, 3D RAS `.nii.gz` image. Preflight checks descriptor and
payload hashes and NIfTI headers without loading voxel arrays. The bundled
LesionTOADS descriptor and the base TOADS descriptor remain profile-selected.
Either profile may provide `run_options`: base fields are names in `BaseRunOptions`;
plus fields are nested `registration`, `n4`, and `plus` objects matching the
corresponding runner configuration fields. Unknown fields are rejected. These
let synthetic runs lower iteration counts without changing algorithm defaults
in the implementation. The plus profile uses the bundled LesionTOADS atlas
descriptor. The base profile uses the bundled default TOADS descriptor. Only
the separate PLUS AtlasEM/Dura handoff descriptor is selectable in the config.

## `pymacruise workflow-run`

```text
pymacruise workflow-run --config workflow.json [--resume]
```

Executes the restartable workflow selected by the same version-1 JSON
configuration. Base runs use the bundled TOADS descriptor; plus runs use the
bundled LesionTOADS descriptor, optionally run the configured AtlasEM/Dura
handoff, and pass `plus_options` through to the existing MaCRUISEplus stages.
The output directory must be new or contain a compatible
restartable run. If it contains the `run.json` produced by
`workflow-preflight`, that dry-run manifest is preserved as `preflight.json`
before execution begins. `--resume` verifies existing checkpoints and resumes
the configured workflow. Non-RAS inputs are normalized by the runners with
visible warnings; persisted image outputs are RAS `.nii.gz`.

The base profile executes eight checked stages: N4, first TOADS, Dura removal,
second TOADS, MaSeg/TOADS tissue preparation, CRUISE/ACE/GVF/TGDM core, MaSeg
refinement, and scientific surface/thickness products. The last stage writes
inner, central, and outer RAS-mm VTK meshes, a measured inner VTK mesh, a
level-set-subtraction thickness NIfTI, a thickness summary JSON, and
provenance. Surface extraction now uses the source-derived connectivity-18
case templates for JIST `IsoSurfaceOnGrid`; selected same-level-set Java
comparisons pass the surface position gates. The selected surface
closest-point and level-set thickness operations are described in
[`thickness-port.md`](thickness-port.md).

## `pymacruise base-run`

```text
pymacruise base-run --t1 T1.nii.gz --brain-mask BrainMask.nii.gz \
  --maseg MaSeg.nii.gz --atlas-descriptor ATLAS.txt --output-dir RUN_DIR \
  [--regrid-inputs | --no-regrid-inputs] [--resume] \
  [--registration-levels N] [--registration-iterations N] \
  [--post-refresh-registration-iterations N] [--registration-max-iterations N] \
  [--toads-iterations N] [--dura-iterations N] \
  [--dura-membership-level VALUE] [--crop-edge-width N] \
  [--ace-diffusion-iterations N] [--ace-use-updated-csf-in-ace {0,1,2,3}] \
  [--gvf-iterations N] [--tgdm-iterations N] \
  [--tgdm-initial-iterations N] [--tgdm-inner-iterations N] \
  [--tgdm-central-iterations N] [--wm-iso-level VALUE] \
  [--source-faithful-refinement | --no-source-faithful-refinement]
```

Regridding is recommended and enabled by default. It resamples an oblique RAS
grid to an orthogonal RAS grid at the T1 voxel spacing. T1 uses linear
interpolation; the binary mask and MaSeg labels use nearest-neighbor
interpolation. The shared regridded inputs and hashes are recorded under the
run's `resources` directory and manifest. Axis-aligned RAS inputs are left
unchanged. Use `--no-regrid-inputs` for source-style axis-only geometry: each
input is reoriented to the closest RAS axis order using only axis permutations
and flips, while oblique direction cosines are ignored. No interpolation is
performed; the manifest records this geometry mode. The run uses the required
user-provided T1-grid mask and does not invoke BET or SPECTRE.

## `pymacruise base-surface-thickness`

```text
pymacruise base-surface-thickness --inner-levelset INNER.nii.gz \
  --central-levelset CENTRAL.nii.gz --outer-levelset OUTER.nii.gz \
  --output-dir OUTPUT_DIR [--provenance RUN.json]
```

This stage adapter accepts finite 3D `.nii.gz` level sets, reorients non-RAS
inputs to RAS with visible warnings, and requires all three images to share
one physical grid. It extracts all three scientific meshes, computes the
selected absolute level-set subtraction thickness volume, and writes the
selected closest-point thickness values on the inner surface. VTK point
coordinates are RAS+ millimetres; the thickness volume is RAS `.nii.gz`.
`thickness_summary.json` records means, sample standard deviations,
crossing/nonzero counts, and the source-derived JIST connectivity-18 surface
method. Undefined summary values are represented as JSON `null`. This
command does not render surfaces or create viewer-specific files.

## `pymacruise topology-correct`

```text
pymacruise topology-correct --volume INPUT.nii.gz --paint-mask MASK.nii.gz \
  --output CORRECTED.nii.gz [--multi-atlas-wm WM.nii.gz \
  --atlas-refined-output REFINED.nii.gz] [--paint-threshold VALUE] \
  [--lowest-level VALUE] [--highest-level VALUE] [--regularize VALUE] \
  [--provenance run.json]
```

This command delegates to the existing source-matched `(6,18)` topology
correction API. It accepts a 3D scalar image and paint mask on the same
physical grid. The optional multi-atlas WM mask and refined output must be
provided together; this enables the selected V3 WM post-processing stage.
Inputs are normalized to RAS with a visible warning for non-RAS orientations.
The corrected output (and optional refined output) is compressed NIfTI in
RAS. Optional provenance records inputs, parameters, outputs, and the RAS
convention.

## `pymacruise segmentation-refine`

```text
pymacruise segmentation-refine --segmentation MaSeg.nii.gz \
  --inner-levelset INNER.nii.gz --outer-levelset OUTER.nii.gz \
  --output REFINED.nii.gz --summary rules.json \
  [--inner-coefficient-mm VALUE] [--outer-coefficient-mm VALUE] \
  [--inner-central-coefficient-mm VALUE] \
  [--outer-central-coefficient-mm VALUE] [--neighborhood-radius VOXELS] \
  [--max-sweeps COUNT] [--deal-with-central | --no-deal-with-central] \
  [--source-faithful | --no-source-faithful] \
  [--provenance run.json]
```

The three inputs must be 3D `.nii.gz` images on the same physical grid after
RAS reorientation. `--summary` records the selected mode and per-rule
candidate, changed, unresolved, sweep, and threshold counts. The default
`source_faithful` behavior matches the V3 MATLAB positive-integer `tabulate`
result: `maxind` equals the actual modal label. `--no-source-faithful` remains
accepted for compatibility and produces the same labels for the supported
candidate groups. The label-refinement output is
written as RAS `.nii.gz`; each non-RAS input produces a visible warning.

## `pymacruise prepare-raw`

```text
pymacruise prepare-raw --t1 T1.nii.gz --maseg MaSeg.nii.gz \
  --brain-mask BrainMask.nii.gz --output-dir OUTPUT_DIR [--overwrite]
```

### Inputs

| Option | Required | Description and assumptions |
|---|---:|---|
| `--t1 PATH` | Yes | 3D, single-volume T1-weighted NIfTI image (`.nii` or `.nii.gz`) with finite scalar values and a usable affine. Outputs are persisted as `.nii.gz`. |
| `--maseg PATH` | Yes | 3D MaSeg NIfTI label map with finite integer-valued voxels. It must match T1 shape and physical affine after RAS orientation normalization. Label values are preserved. |
| `--brain-mask PATH` | Yes | User-supplied 3D binary NIfTI mask for the T1 image. Values must be 0/1 or 0/255, with at least one foreground voxel. Shape and physical affine must match T1 after RAS normalization. The mask is not resampled. |
| `--output-dir PATH` | Yes | Directory for prepared inputs. Existing outputs cause an error unless `--overwrite` is supplied. |
| `--overwrite` | No | Allow replacing `T1.nii.gz`, `Seg.nii.gz`, and `BrainMask.nii.gz` in the output directory. |

Each non-RAS input is reoriented to RAS and produces a visible warning. The
brain mask is the selected skull-mask source: this command does not run BET or
SPECTRE and does not require FSL or FreeSurfer.

### Outputs

| File | Description and properties |
|---|---|
| `T1.nii.gz` | Prepared T1 in RAS voxel orientation, preserving physical voxel locations and source scalar values. |
| `Seg.nii.gz` | Prepared MaSeg labels on the T1 RAS grid, preserving integer labels and dtype. |
| `BrainMask.nii.gz` | Prepared binary brain mask on the exact T1 RAS grid, preserving 0/1 or 0/255 encoding. |

The Python API returns source orientation codes for the three inputs. The full
workflow manifest will also record source orientations and conversions when
P34/P42a workflow orchestration is implemented.

## `pymacruise atlasem-atlas-inventory`

```text
pymacruise atlasem-atlas-inventory --descriptor atlas.txt
```

This reads descriptor metadata without loading image voxels. The report lists
the structure order, referenced topology/shape/intensity files, and whether
each referenced image is present as `.nii.gz`. It is useful for checking an
atlas bundle before a run; it does not establish that the selected AtlasEM
atlas is available or validate payload image contents.

## `pymacruise atlasem-prelesion-dura`

```text
pymacruise atlasem-prelesion-dura --t1 T1.nii.gz \
  --memberships MEMBERSHIPS.nii.gz --brain-mask BrainMask.nii.gz \
  --class-names '["BG","CSF","GM","WM"]' --output-dir OUTPUT_DIR \
  [--max-iterations 5] [--provenance run.json]
```

### Inputs

| Option | Required | Description and assumptions |
|---|---:|---|
| `--t1 PATH` | Yes | 3D scalar T1 NIfTI; this is the source of the output grid and voxel values. |
| `--memberships PATH` | Yes | 4D float membership NIfTI with classes on the final axis. Its spatial grid must match the T1 physical grid after RAS orientation normalization. |
| `--brain-mask PATH` | Yes | 3D binary mask matching the T1 physical grid. The mask is required and is not resampled. |
| `--class-names JSON` | Yes | JSON string array in membership-channel order. Names control selected GM/WM/Dura output mapping. |
| `--output-dir PATH` | Yes | Destination for compressed RAS NIfTI products and the restart manifest. |
| `--max-iterations COUNT` | No | Maximum selected Dura outer passes; default is 5. |
| `--provenance PATH` | No | Optional JSON sidecar recording inputs, effective parameters, and outputs. |

The command consumes existing AtlasEM memberships. It does not load an atlas,
estimate registration, or run the AtlasEM EM loop. Non-RAS images are
reoriented with a warning. The output grid matches the input T1 in RAS.

### Outputs

The command writes `prelesion_dura_inputs.nii.gz`,
`prelesion_dura_segmentation.nii.gz`, `dura_stripped_t1.nii.gz`,
`dura_corrected_white_matter.nii.gz`, `dura_corrected_gray_matter.nii.gz`,
`dura_brain_mask_level_set.nii.gz`, `dura_brain_mask.nii.gz`,
`dura_wm_topology.nii.gz`, and `dura_gm_topology.nii.gz`. All are compressed
NIfTI on the native T1 grid in RAS. `run.json` records hash-verified AtlasEM
output and Dura stages for resume.

## `pymacruise atlasem-run`

```text
pymacruise atlasem-run --t1 T1.nii.gz --brain-mask BrainMask.nii.gz \
  --descriptor atlas.txt --output-dir OUTPUT_DIR [--modality T1_MPRAGE] \
  [--background-ratio 0.0] [--crop-border 10] [--levels 4] \
  [--initial-registration-iterations 50] \
  [--post-recrop-registration-iterations 1] [--max-em-iterations 99] \
  [--max-membership-change 0.01] [--smoothing 10.0] \
  [--max-dura-iterations 5] [--resume | --overwrite] [--provenance run.json]
```

### Inputs and assumptions

| Option | Required | Description and assumptions |
|---|---:|---|
| `--t1 PATH` | Yes | 3D scalar T1 NIfTI. AtlasEM currently supports one T1 channel. |
| `--brain-mask PATH` | Yes | Binary brain mask on the T1 physical grid; required in place of BET/SPECTRE and never resampled. |
| `--descriptor PATH` | Yes | AtlasEM structure descriptor. Every referenced topology and shape payload must be available as compressed NIfTI; descriptor class order determines membership order. |
| `--output-dir PATH` | Yes | Destination for AtlasEM and pre-LesionTOADS Dura products and the stage manifest. |
| `--modality NAME` | No | Intensity-atlas modality key; default is `T1_MPRAGE`. |
| `--background-ratio VALUE` | No | Source preprocessing background ratio; default is `0.0`. |
| `--crop-border VOXELS` | No | Padding added to the selected AtlasEM crop; default is 10 voxels. |
| `--levels COUNT` | No | Number of registration pyramid levels; default is 4. |
| `--initial-registration-iterations COUNT` | No | Iteration cap per initial pyramid level; default is 50. |
| `--post-recrop-registration-iterations COUNT` | No | Iteration cap for the full-resolution recropped rigid stage; default is 1. |
| `--max-em-iterations COUNT` | No | Fixed-prior EM iteration cap; default is 99. |
| `--max-membership-change VALUE` | No | EM stop threshold for maximum membership change; default is 0.01. |
| `--smoothing VALUE` | No | Selected membership smoothing parameter; default is 10.0. |
| `--max-dura-iterations COUNT` | No | Maximum pre-LesionTOADS Dura outer passes; default is 5. |
| `--resume` | No | Reuse hash-verified stages when inputs, settings, assets, and implementation match. |
| `--overwrite` | No | Start a new run by replacing the run manifest and expected product files in the output directory. |
| `--provenance PATH` | No | Optional JSON sidecar recording input hashes, settings, stage counts, and outputs. |

All spatial outputs are written as native-grid RAS `.nii.gz`. Non-RAS input
orientations produce visible warnings and are normalized to RAS. The command
executes the supported single-channel descriptor-driven AtlasEM registration,
fixed-prior EM, output generation, and pre-LesionTOADS Dura handoff. Current
execution uses synthetic descriptor cases; the selected `ibsr8obj` descriptor
and payload are unavailable, and direct Java optimizer parity remains open.

### Outputs

AtlasEM products are `atlasem_memberships.nii.gz`,
`atlasem_segmentation.nii.gz`, `atlasem_internal_mask.nii.gz`, and
`atlasem_dura_inputs.nii.gz`. The Dura stage writes the nine products listed
under `atlasem-prelesion-dura`. The command also writes a hash-checked
`run.json` stage manifest; `--resume` verifies and skips valid AtlasEM/Dura
stages independently.

## `pymacruise calculate`

```text
pymacruise calculate --input-a A.nii.gz --input-b B.nii.gz \
  --operation OP --output RESULT.nii.gz [--provenance run.json]
```

Applies a voxel-wise N-dimensional array operation after RAS axis orientation.
Inputs must have identical shapes and physical affines; this command does not
register or resample. Supported operations are `ADD`, `SUBTRACT`, `MULTIPLY`,
`DIVIDE`, `MIN`, and `MAX`. Arithmetic is not clamped or intensity-scaled.
The result is compressed RAS NIfTI, preserving image A's header and promoting
integer inputs to a floating type to avoid overflow. Division follows NumPy
IEEE behavior for zero denominators.

## `pymacruise cortical-separation`

```text
pymacruise cortical-separation --maseg MaSeg.nii.gz --output-dir OUTPUT_DIR
```

Reads an integer-valued 3D MaSeg label image, reorients non-RAS inputs with a
warning, and writes the four selected masks `Atlas_GM.nii.gz`,
`Atlas_WM.nii.gz`, `Atlas_Mask.nii.gz`, and `Atlas_BrainTissueMask.nii.gz` as
uint8 0/255 RAS images. The selected source intentionally defines
`Atlas_Mask` from the WM label set.

## `pymacruise base-resume-from-toads`

```text
pymacruise base-resume-from-toads --source-run SOURCE/run.json \
  --maseg MaSeg.nii.gz --output-dir DISTINCT_TAIL_DIR [--resume]
```

Imports a completed, hash-verified source run through its second TOADS stage,
then resumes the selected BASE tail using the supplied same-grid MaSeg. The
source run is read-only; tail outputs and their manifest go to the distinct
output directory. `--resume` reuses compatible verified tail checkpoints.
Images written by downstream stages use compressed RAS NIfTI. This command
does not establish a new source TOADS run.

## Other currently exposed commands

- `pymacruise --help` lists the current command groups and their required
  options.
- `pymacruise --version` reports the installed package version.
The base profile's restart behavior and final-output manifest are exercised by
the synthetic workflow tests; plus-profile source equivalence remains open.
See
[`io-contract.md`](io-contract.md) for profile-level input/output semantics and
[`toads-port-plan.md`](toads-port-plan.md) for the TOADS blocker plan.

## Supported standalone scientific tools

Image inputs and outputs for these commands use `.nii.gz`. The image APIs
normalize inputs to RAS and warn when they reorient an input. Surface
coordinates are RAS+ millimetres and use tagged ASCII legacy VTK. These
commands expose selected standalone stages; they do not assemble the complete
base MaCRUISE or LesionTOADS pipelines.

Each of the eight standalone commands also accepts optional
`--provenance RUN.json`. After successful processing it writes a JSON sidecar
with a schema version, command name, PyMaCRUISE version, resolved input and
output paths, the `RAS` spatial convention, and effective parameters
(including API defaults used by that command). Output directories are
expanded into the individual output image paths. Optional inputs that were
not supplied are omitted.

### `pymacruise n4`

```text
pymacruise n4 --input T1.nii.gz --output T1-n4.nii.gz \
  [--shrink-factor N] [--iterations N [N ...]] \
  [--convergence-tolerance VALUE] [--spline-mesh X Y Z] \
  [--rescale-intensities | --no-rescale-intensities] [--provenance RUN.json]
```

Reads a finite, non-negative 3D intensity image and writes its N4-corrected
image as `.nii.gz` in RAS orientation. The defaults match the selected
wrapper's unadorned ANTs N4 command: shrink factor 4, four 50-iteration
levels, zero convergence tolerance, a 1×1×1 spline mesh, and intensity
rescaling off because the source wrapper supplies no mask (enable explicitly
with `--rescale-intensities`). ANTsPy is optional; install PyMaCRUISE with its `ants` extra
to run this command. The N4 estimation itself is provided by ANTsPy/ANTs.

### `pymacruise ace`

```text
pymacruise ace --gray GM.nii.gz --white WM.nii.gz \
  --enhanced GM-enhanced.nii.gz --skeleton ACE-skeleton.nii.gz \
  --thinned-skeleton ACE-skeleton-thinned.nii.gz [--atlas-csf CSF.nii.gz] \
  [--diffusion-iterations 25] [--provenance RUN.json]
```

Inputs are aligned 3D gray/white matter membership or intensity images; the
optional CSF atlas image must share their physical grid. Writes enhanced gray
matter and both ACE skeleton products. `--diffusion-iterations` defaults to
25. This exposes the current ACE stage, not all multi-atlas MaACE options or
the complete MaCRUISE workflow.

### `pymacruise gvf`

```text
pymacruise gvf --input IMAGE.nii.gz --output GVF.nii.gz \
  [--lambda VALUE] [--max-iterations N] [--solver-tolerance VALUE] \
  [--provenance RUN.json]
```

Reads a finite 3D scalar image and writes a RAS-oriented, component-last
three-vector NIfTI image. Defaults are lambda 1.0, 200 source multigrid
V-cycles, and absolute maximum-residual threshold 1e-6.

### `pymacruise tgdm`

```text
pymacruise tgdm --initial INITIAL.nii.gz --white WM.nii.gz \
  --gray GM.nii.gz --ace-gray GM-ACE.nii.gz --output-dir DIR \
  [--brain-mask MASK.nii.gz] [--ace-skeleton SKELETON.nii.gz] \
  [--gvf GVF.nii.gz] [--provenance RUN.json]
```

The four required finite 3D scalar images and each supplied optional image
must share one physical grid. Optional GVF is a component-last vector image.
Writes `initial_level_set.nii.gz`, `tissue_sign.nii.gz`,
`initial_inner.nii.gz`, `inner.nii.gz`, `central.nii.gz`, and `outer.nii.gz`,
all on the RAS grid. This is the nested TGDM stage, not full TOADS or MaCRUISE.

### `pymacruise dura`

```text
pymacruise dura --t1 T1.nii.gz --white WM.nii.gz --gray GM.nii.gz \
  --white-mask WM-mask.nii.gz --output-dir DIR [--max-iterations N] \
  [--membership-level 0.3] [--wm-topology-regularization 0.05] \
  [--gm-topology-regularization 0.02] [--provenance RUN.json]
```

The T1, white/gray membership, and original white-matter mask must be finite
3D images on one physical grid. `--max-iterations` defaults to 5. Writes
RAS-oriented `stripped_t1.nii.gz`, `corrected_white_matter.nii.gz`,
`corrected_gray_matter.nii.gz`, `brain_mask_level_set.nii.gz`,
`brain_mask.nii.gz`, `wm_topology.nii.gz`, and `gm_topology.nii.gz`.

### `pymacruise surface-extract`

```text
pymacruise surface-extract --level-set LEVELSET.nii.gz --output SURFACE.vtk \
  [--level VALUE] [--provenance RUN.json]
```

Reads a finite 3D scalar level set and writes a triangular RAS+ millimetre
surface. The isovalue defaults to zero. The VTK header declares coordinate
frame and units.

### `pymacruise thickness-levelset`

```text
pymacruise thickness-levelset --inner INNER.nii.gz --outer OUTER.nii.gz \
  --output THICKNESS.nii.gz [--iso-level VALUE] [--provenance RUN.json]
```

Inputs must share one physical RAS grid. Writes the voxelwise absolute
level-set difference in millimetres and prints mean, sample standard
deviation, and crossing-voxel count for the selected iso-level (default 0).

### `pymacruise thickness-surface`

```text
pymacruise thickness-surface --inner INNER.vtk --outer OUTER.vtk \
  --inner-output INNER-MEASURED.vtk [--central CENTRAL.vtk \
  --central-output CENTRAL-MEASURED.vtk] [--provenance RUN.json]
```

Inputs must be triangular ASCII legacy VTK with a declared RAS/LPS frame.
Writes closest-surface point measurements into output meshes and prints
summary statistics. Supply `--central` and `--central-output` together.

## Not currently exposed as complete tools

The `toads-loop` command below runs the selected registered TOADS loop, but it
does not complete the two-pass MaCRUISE TOADS/dura workflow. There is no full
LesionTOADS command. The user-supplied mask and native-T1 RAS working space are
selected; historical SPECTRE preprocessing remains an explicit working-space
deviation documented in the port plan.


## `pymacruise toads-loop`

```text
pymacruise toads-loop --t1-n4 T1_N4.nii.gz --brain-mask BrainMask.nii.gz \
  --output Classification.nii.gz [--atlas-descriptor ATLAS.txt] \
  [--modality T1_MPRAGE] [--registration-levels 4] \
  [--registration-iterations 50] [--post-refresh-registration-iterations 5] \
  [--registration-max-iterations 5] \
  [--padding 10] [--maximum-iterations 50] \
  [--convergence-threshold 1e-4] [--relation-limit 3] \
  [--maximum-distance 1e-4] [--category-exponents CSF GM WM] \
  [--source-faithful-nonfinite-registration-trials] \
  [--preserve-input-axis-order] \
  [--dura-output-dir DURA_INPUTS_DIR] [--provenance run.json]
```

### Inputs and assumptions

| Option | Required | Description and assumptions |
|---|---:|---|
| `--t1-n4 PATH` | Yes | 3D single-channel T1 already corrected by the standalone `n4` command; compressed `.nii.gz`. This command does not run N4. |
| `--brain-mask PATH` | Yes | User-provided binary brain mask on the T1 physical grid. Non-RAS T1 or mask images are reoriented with a visible warning. |
| `--atlas-descriptor PATH` | No | Defaults to the bundled base TOADS descriptor. Referenced atlas images must remain beside the descriptor using its declared relative paths and `.nii.gz` format. |
| `--output PATH` | Yes | Native-grid terminal classification, compressed `.nii.gz`; existing output is not overwritten. |
| `--dura-output-dir PATH` | No | Optional directory for the classification and four source-scaled float32 Dura-removal inputs (CSF, GM, WM fill, WM mask), all restored to native RAS. Existing output collisions are rejected. |
| Registration and loop options | No | Selected source controls with defaults shown above. Post-refresh alignment defaults to the Java `alignIter=1`; it is separate from startup-pyramid registration iterations and per-loop registration iterations. CSF, GM, and WM exponents default to the atlas `Optimized factors` row (or Java's all-ones defaults when omitted); `--category-exponents` overrides those three values for the single T1 channel. |
| `--source-faithful-nonfinite-registration-trials` | No | Enable Java-style edge interpolation for NaN registration trials; degenerate inputs can propagate a nonfinite transform. |
| `--preserve-input-axis-order` | No | Reconstruct axial RAS/LPS voxel order internally for source-compatibility while writing images in RAS. |

### Output

The command writes one uint8 classification map on the original T1 grid with
RAS orientation and affine. Labels are restored from the registered crop using
the selected TOADS padding and x-fast uncrop semantics. The optional Dura-removal files are inputs for a later stage; this command
does not execute dura removal, run the second MaCRUISE TOADS pass, create final
CRUISE output products, or execute a complete base/plus workflow by itself.
The separate `base-run` and `plus-full-run` commands assemble restartable
end-to-end workflows; synthetic workflow execution is validated, while
whole-workflow Java/MATLAB equivalence remains open.

For the selected synthetic handoff, use the three exported maps required by
the current Dura CLI input contract and retain the original T1 as its reference:

```text
pymacruise toads-loop --t1-n4 T1_N4.nii.gz --brain-mask BrainMask.nii.gz \
  --output Classification.nii.gz --dura-output-dir toads_dura_inputs
pymacruise dura --t1 T1.nii.gz \
  --white toads_dura_inputs/T1_N4_wmfill.nii.gz \
  --gray toads_dura_inputs/T1_N4_gm.nii.gz \
  --white-mask toads_dura_inputs/T1_N4_wmmask.nii.gz \
  --output-dir dura_outputs --max-iterations 1
```

The synthetic integration test compares the seven Dura outputs with the direct
API after one evolution pass. This example does not perform the second TOADS
pass.

The selected second-pass handoff can be exercised with the Dura stripped T1
and corrected brain mask:

```text
pymacruise toads-loop --t1-n4 dura_outputs/stripped_t1.nii.gz \
  --brain-mask dura_outputs/brain_mask.nii.gz \
  --output second_pass/Classification.nii.gz \
  --dura-output-dir second_pass/dura_inputs
```

This uses the current native-RAS working space. The historical SPECTRE
registration and inverse transfer are not implemented in this handoff.

## `pymacruise remove-cerebellum`

```text
pymacruise remove-cerebellum --segmentation MaSeg.nii.gz \
  --gm Toads2_gm.nii.gz --wm-fill Toads2_wmfill.nii.gz \
  --wm-mask Toads2_wmmask.nii.gz --output-dir cerebrum_maps \
  [--provenance run.json]
```

The CLI nearest-neighbor transfers the integer MaSeg label map to each TOADS
map grid, then zeros labels 35, 38–41, and 71–73. It preserves each input
map's scalar dtype and writes three RAS `.nii.gz` files. Physical-affine
transfer has only been compared with the selected JIST operation on identical
RAS grids; arbitrary origin/spacing changes remain unvalidated.

## `pymacruise crop-tissue-maps`

```text
pymacruise crop-tissue-maps --gm Toads2_gm.nii.gz \
  --wm-fill Toads2_wmfill.nii.gz --wm-mask Toads2_wmmask.nii.gz \
  --output-dir cropped_maps [--edge-width 5] [--provenance run.json]
```

Positive GM voxels determine one inclusive crop box. The configured edge is
added on all sides, with zero padding beyond the source grid; the same bounds
are applied to all three maps. Outputs retain their input types and are saved
as RAS `.nii.gz` files. The default edge width is the selected wrapper value
of five voxels.

## `pymacruise prepare-cruise-inputs`

```text
pymacruise prepare-cruise-inputs --gm Toads2_gm.nii.gz \
  --wm-fill Toads2_wmfill.nii.gz --wm-mask Toads2_wmmask.nii.gz \
  [--segmentation MaSeg.nii.gz --remove-cerebellum] \
  --output-dir cruise_inputs [--edge-width 5] [--provenance run.json]
```

This composes the selected postpass operations in order: optional MaSeg
brainstem/cerebellum removal, one positive-GM-derived crop shared by all tissue
maps, then (6,18) topology correction of cropped WM-fill painted by cropped
WM-mask. It writes the three cropped maps and corrected WM-fill on native
cropped RAS geometry. Base workflow defaults are used; multi-atlas topology
refinement is not enabled by this command.

## `pymacruise base-tissue-preparation`

```text
pymacruise base-tissue-preparation --maseg MaSeg.nii.gz \
  --gm Toads2_gm.nii.gz --wm-fill Toads2_wmfill.nii.gz \
  --wm-mask Toads2_wmmask.nii.gz --output-dir prepared_maps \
  --provenance run.json [--edge-width 5] \
  [--gm-confident-level 0.001] [--remove-cerebellum]
```

Derives MaSeg cortical tissue maps, transfers the three second-pass TOADS maps
to the MaSeg grid with nearest-neighbor interpolation, crops all maps from
positive MaSeg GM support, then fuses them with the selected
`CombineMAandTOADS` settings. It writes the intermediate atlas/TOADS maps and
the fused GM, WM-fill, and mask as RAS `.nii.gz` files. Optional
brainstem/cerebellum removal is disabled in the selected V3 base workflow.
This fusion mask is derived from MaSeg and is distinct from the external
brain mask used to skull-strip TOADS.

## `pymacruise toads-cruise-core`

```text
pymacruise toads-cruise-core --gm Toads2_gm.nii.gz \
  --wm-fill Toads2_wmfill.nii.gz --wm-mask Toads2_wmmask.nii.gz \
  [--segmentation MaSeg.nii.gz --remove-cerebellum] \
  --output-dir cruise_core [--edge-width 5] \
  [--ace-diffusion-iterations 25] [--gvf-lambda 1] \
  [--gvf-iterations 200] [--gvf-solver-tolerance 1e-6] \
  [--boundary-csf CSF.nii.gz] [--use-updated-csf-in-ace {0,1,2,3}] \
  [--tgdm-initial-iterations 6] \
  [--tgdm-inner-iterations 4] [--tgdm-central-iterations 7] \
  [--base-preparation-provenance run.json] \
  [--provenance run.json]
```

This continues from second-pass maps through cleanup/crop/topology,
base-parameter ACE, generalized GVF, nested TGDM, and native-grid inverse
crop. When `--base-preparation-provenance` is supplied, the command instead
accepts already-fused, cropped maps, verifies their grid against the recorded
MaSeg crop, performs the selected WM-fill topology correction, and avoids
repeating cleanup/crop. It writes cropped tissue maps, topology-corrected
WM-fill, ACE outputs, the component-last GVF, and native-RAS
inner/central/outer level sets. Base ACE does not merge updated multi-atlas
CSF. The selected MATLAB wrapper passes the central-surface iteration option
to the outer TGDM stage as well, so the CLI preserves that mapping. Surfaces
and thickness outputs are not generated.


## `pymacruise plus-run`

```text
pymacruise plus-run --t1-reference T1.nii.gz --segmentation MaSeg.nii.gz \
  --dura-csf CSF.nii.gz --dura-gm GM.nii.gz --dura-wm-fill WMFill.nii.gz \
  --dura-wm-mask WMMask.nii.gz --multi-atlas-gm AtlasGM.nii.gz \
  --multi-atlas-wm AtlasWM.nii.gz \
  --multi-atlas-mask Atlas_BrainTissueMask.nii.gz \
  --brain-mask BrainMask.nii.gz --output-dir plus-run [--resume] \
  [--dura-transfer-mode physical-affine|jist-unchanged] \
  [--crop-edge-width 5] [--gvf-lambda 0.2] [--wm-iso-level 0.5] \
  [--refine-topology-by-maseg | --no-refine-topology-by-maseg]
```

This restartable continuation starts from the persisted four-channel P39 Dura
products and already aligned segmentation/atlas/brain-mask inputs. It records
and verifies three stages: channel transfer and crop; plus topology, ACE, GVF,
and TGDM; and full-grid restoration, surfaces, and thickness. `--resume`
requires the same input hashes and settings. Verified stages are skipped; an
invalid stage is recomputed in a new attempt directory. Outputs include RAS
`.nii.gz` level sets and thickness, RAS+ millimetre VTK surfaces, thickness
metrics, and `run.json`. The default `physical-affine` mode expects Dura maps
on the full T1 grid. `jist-unchanged` accepts common smaller Dura extents that
carry the T1 reference header, then applies the MaSeg crop bounds with padding.

`--multi-atlas-mask` must be the MaSeg-derived `Atlas_BrainTissueMask.nii.gz`
(GM, WM, and additional brain-tissue labels). It must not be the WM-only
`Atlas_Mask.nii.gz`; that WM-only product remains distinct in
`pymacruise cortical-separation` outputs.

The command does not run FLAIR registration, atlas registration, TOADS, or
LesionTOADS. Use `pymacruise plus-register-flair` for the separate ANTsPy rigid
FLAIR-to-T1 checkpoint. These commands are workflow slices, not a complete
MaCRUISEplus run or a claim of equivalence to NiftyReg/Aladin or the Java
workflow.


## `pymacruise plus-register-flair`

```text
pymacruise plus-register-flair --t1 T1.nii.gz --flair FLAIR.nii.gz \
  --output-dir registration-checkpoint [--resume] \
  [--metric mattes|meansquares|GC] [--metric-bins 32] \
  [--sampling-rate 0.2] [--iterations 1000 500 250 0] \
  [--shrink-factors 6 4 2 1] [--smoothing-sigmas 3 2 1 0]
```

T1 is the fixed image and FLAIR is the moving image. The command normalizes
inputs to RAS with warnings for reorientation, then uses the existing rigid
ANTsPy adapter to write `flair_in_t1.nii.gz`, forward/inverse transform files,
and `flair-registration.json`. Resume verifies input/settings signatures,
output hashes, and the warped image grid before skipping registration.
ANTsPy is an explicit Python backend substitute; NiftyReg/Aladin numerical
equivalence is unvalidated. This checkpoint alone does not run N4, TOADS,
LesionTOADS, or the full plus workflow.

## `pymacruise plus-prepare-lesion`

```text
pymacruise plus-prepare-lesion --t1 T1.nii.gz --flair FLAIR.nii.gz \
  --brain-mask BrainMask.nii.gz --output-dir lesion-inputs [--resume] \
  [--metric mattes|meansquares|GC] [--metric-bins 32] \
  [--sampling-rate 0.2] [--registration-iterations 1000 500 250 0] \
  [--shrink-factors 6 4 2 1] [--smoothing-sigmas 3 2 1 0] \
  [--n4-shrink-factor 4] [--n4-iterations 50 50 50 50] \
  [--n4-convergence-tolerance 0] [--n4-spline-mesh 1 1 1]
```

This restartable input-preparation checkpoint follows the selected plus order:
rigidly register FLAIR to T1 with ANTsPy, N4-correct T1 and registered FLAIR,
then apply the supplied T1-grid brain mask and set negative values to zero. It
writes `t1_brain.nii.gz`, `flair_brain.nii.gz`, a 4-D channel-last
`modalities_t1_flair.nii.gz`, and a hash-checked manifest. Inputs must be
compressed 3-D NIfTI; T1 and mask must share a grid, while FLAIR may start on
another grid. Outputs use RAS and modality order is exactly `T1, FLAIR`. The
ANTsPy registration is not numerically validated against NiftyReg/Aladin. This
command prepares LesionTOADS inputs; it does not execute LesionTOADS, the Dura
workflow, or complete MaCRUISEplus.

## `pymacruise plus-crop-lesion-inputs`

```text
pymacruise plus-crop-lesion-inputs --modalities T1Flair.nii.gz \
  --topology LesionTopology.nii.gz --output-dir lesion-crop \
  [--padding 10] [--resume]
```

Consumes the 4-D channel-last T1-then-FLAIR image from
`plus-prepare-lesion` and a 3-D LesionTOADS topology atlas. It applies the
selected source normalization, unions nonzero image support with transformed
topology support, crops and zero-pads both channels, then records inclusive
source bounds and the post-crop rigid translation. It writes whole-grid
normalized and cropped RAS `.nii.gz` images plus a hash-checked manifest.
`--padding` defaults to the source-selected ten voxels. Resume verifies both
input/settings signature and output files before skipping. This is input
preparation for initial atlas alignment; it does not execute the iterative
LesionTOADS workflow.

## `pymacruise plus-lesion-initial-align`

```text
pymacruise plus-lesion-initial-align --cropped-modalities T1FlairCrop.nii.gz \
  --crop-manifest crop.json --atlas-descriptor LesionTOADS-atlas.txt \
  --output-dir initial-alignment [--init-align-iterations 10] \
  [--pyramid-levels 2] [--resume]
```

Consumes the crop checkpoint and its manifest, then initializes source-ordered
13-class T1/FLAIR memberships and the selected rigid atlas alignment. The
cropped image must be RAS, 4-D channel-last `.nii.gz`, with T1 then FLAIR; the
manifest must match its hash and describe the same crop geometry. The atlas
descriptor and its referenced data must be available. Outputs are a
channel-last 4-D RAS membership image, a transform/call-order summary, and a
hash-checked restart manifest. The defaults use ten initial alignment
iterations and two pyramid levels; these remain explicit settings rather than
a claim of full LesionTOADS equivalence. This checkpoint does not run the
iterative LesionTOADS evolution or complete MaCRUISEplus.

## `pymacruise plus-lesion-registered-recrop`

```text
pymacruise plus-lesion-registered-recrop \
  --normalized-modalities T1FlairNormalized.nii.gz \
  --crop-manifest crop.json \
  --initial-alignment-manifest initial-alignment.json \
  --atlas-descriptor LesionTOADS-atlas.txt --output-dir registered-recrop \
  [--resume]
```

Consumes P41f's full-grid normalized image and P41f/P41g manifests. It verifies
the crop inputs, initial alignment's 13-class membership image and transform
summary, then follows the source order for uncropping with the preprocessing
transform, finding/cropping the new image bounds, and sampling the transformed
topology and shape maps with the optimized atlas transform. It writes
`recropped_modalities_t1_flair.nii.gz`, `transformed_topology.nii.gz`,
`transformed_shapes_13.nii.gz`, a summary, and a hash-checked restart manifest.
All image outputs are RAS `.nii.gz` on one crop grid. Resume verifies the
inputs and outputs before skipping the recrop. The committed small synthetic
atlas fixture keeps routine validation bounded; full LesionTOADS iterations
and real-image behavior remain open.

## `pymacruise plus-lesion-first-iteration`

```text
pymacruise plus-lesion-first-iteration \
  --registered-recrop-manifest registered-recrops.json \
  --atlas-descriptor LesionTOADS-atlas.txt --output-dir first-iteration \
  [--iterations 99] [--first-thinning-limit 3] \
  [--first-thinning-volume-ratio 0.5] [--finite-zero-variance] [--resume]
```

Consumes the registered-recrop outputs and computes the selected standard-13
initial membership pass followed by loop-entry setup and first thinning. It
writes 17-channel memberships (13 classes plus lesion, black hole, outlier,
and lesion relation), first-thinning segmentation/skeleton/availability/
ordering images, a diagnostic summary, and a hash-checked manifest. The
source-faithful default preserves zero-variance IEEE NaNs; if these make the
initial objective nonfinite, the manifest records a failed checkpoint and the
run stops. `--finite-zero-variance` enables an explicit stabilization: matching
centroids with zero deviation have distance zero and nonmatching ones infinite
distance. The summary records zero-deviation channels and nonfinite membership
count so this numerical difference remains visible. All outputs use the
registered RAS grid. This step does not execute later LesionTOADS updates.

## `pymacruise plus-lesion-registered-branch`

```text
pymacruise plus-lesion-registered-branch \
  --registered-recrop-manifest registered-recrops.json \
  --atlas-descriptor LesionTOADS-atlas.txt --output-dir registered-branch \
  [--iterations 99] [--first-thinning-limit 3] \
  [--first-thinning-volume-ratio 0.5] [--max-distance 0.0001] \
  [--first-iteration-manifest first-iteration/lesion-first-iteration.json] \
  [--finite-zero-variance] [--resume]
```

This checkpoint runs the compact selected-plus path from initial registered
membership and first thinning through the first post-thinning statistics and
membership update, atlas registration, and one wrapper branch. It writes the
updated 17-channel membership, segmentation, skeleton, availability, ordering,
relation map, summary, and hash-checked manifest as RAS `.nii.gz`/JSON. The
checkpoint also saves membership-energy diagnostics, initial distance maps,
refreshed centroids, energy-ratio intermediates, and registration bounds for
the next loop pass. The
source-faithful zero-variance behavior remains the default and records a failed
manifest with its variance and nonfinite-energy diagnostics. Use
`--finite-zero-variance` to opt into the explicitly stabilized synthetic
path. This is one selected synthetic branch; later iterations, full plus
execution, Java composite equivalence, and real-image validation remain open.
Passing `--first-iteration-manifest` reuses P41i's persisted memberships,
diagnostic energies, distance maps, relations, and first-thinning state; the
checkpoint must match the supplied recrop manifest and numeric settings.

## `pymacruise plus-lesion-registered-loop`

```text
pymacruise plus-lesion-registered-loop \
  --registered-checkpoint-manifest registered-branch/lesion-registered-branch.json \
  --atlas-descriptor LesionTOADS-atlas.txt --output-dir registered-loop-2 \
  [--iterations 99] [--thinning-limit 3] [--max-passes 1] \
  [--correction-fields fields.nii.gz] [--correction-polynomial-degree {1,2,3}] \
  [--initial-correction-degree N] \
  [--finite-zero-variance] [--resume]
```

Continue from a registered branch checkpoint or a previous registered-loop
checkpoint. `--max-passes` bounds work for this invocation; reaching that
boundary is recorded as a checkpoint, not as source convergence. The manifest
stores the iteration number, RAS outputs, atlas/parent hashes, and selected
settings so a later invocation can continue without replaying the completed
segment. Source-faithful zero-variance behavior remains the default. This is
synthetic selected-loop integration. Optional correction fields must be a
channel-last RAS `.nii.gz` on the registered T1/FLAIR grid. An adaptive degree
ceiling requires fields; each checkpoint writes the updated field volume and
current degree, and the next segment carries both forward automatically.
When resuming the same checkpoint output directory, repeat the original
`--correction-fields` argument if it was used to create that checkpoint.
Natural convergence, complete Java equivalence, and real-image validation
remain open.

## `pymacruise plus-lesion-loop-run`

```text
pymacruise plus-lesion-loop-run \
  --registered-recrop-manifest recrop/registered-recrop.json \
  --atlas-descriptor LesionTOADS-atlas.txt --output-dir lesion-loop \
  [--iterations 99] [--first-thinning-limit 3] \
  [--first-thinning-volume-ratio 0.5] [--max-distance 0.0001] \
  [--loop-convergence-threshold N] [--n-iter-max 99] \
  [--thinning-limit 3] [--max-passes 99] \
  [--correction-fields fields.nii.gz] [--correction-polynomial-degree {1,2,3}] \
  [--initial-correction-degree N] \
  [--finite-zero-variance] [--resume]
```

Compose the persisted P41i initial membership/thinning, P41k registered
branch, and P41m registered-loop iterations from a completed P41h
recrop manifest. A top-level JSON manifest links each child manifest by hash
and records stage state before execution. `--resume` asks each child checkpoint
to verify and reuse its outputs; corrupt loop outputs are recomputed. Source-
faithful zero-variance behavior remains the default. `--finite-zero-variance`
selects the explicit finite synthetic path. An optional
`--loop-convergence-threshold` overrides only P41m's convergence threshold;
`--max-distance` continues to control the P41k registered branch and defaults
to the selected source tolerance `1e-4`. The workflow allows 99 P41m passes by
default, stopping at convergence or `--n-iter-max`; an explicitly smaller
`--max-passes` writes a continuing checkpoint instead. A terminal result
writes classification and lesion-subtype NIfTIs. This
does not establish complete MaCRUISEplus, Java composite equivalence, or
real-image validation.

Optional correction fields must be channel-last, RAS, compressed NIfTI on the
registered recrop T1/FLAIR grid. An adaptive degree ceiling requires these
fields. The registered-loop checkpoint writes updated fields and their degree,
so later segments resume them; `--resume` with the same workflow also requires
the original input field path if it was supplied.

## `pymacruise plus-full-run`

Run the restartable raw-input plus path from T1, FLAIR, a user-provided T1-grid brain mask, MaSeg labels, and a descriptor-owned atlas through preprocessing, registered LesionTOADS, Dura-channel transfer, and terminal CRUISE surfaces and thickness.

```bash
pymacruise plus-full-run --t1 T1.nii.gz --flair FLAIR.nii.gz \
  --brain-mask brain-mask.nii.gz --maseg MaSeg.nii.gz \
  --atlas-descriptor atlas.txt \
  [--atlasem-descriptor AtlasEM.txt] \
  --output-dir plus-output \
  --registration-iterations 1000 500 250 0 --registration-sampling-rate 0.25 \
  --shrink-factors 6 4 2 1 --smoothing-sigmas 3 2 1 0 \
  --n4-iterations 50 50 50 50 \
  [--regrid-inputs | --no-regrid-inputs] \
  [--dura-transfer-mode {physical-affine,jist-unchanged}] [--resume]
```

The full plus workflow accepts the same `[--regrid-inputs | --no-regrid-inputs]`
option, recommended and enabled by default. When T1 is oblique, it creates a
shared orthogonal RAS grid for T1, brain mask, and MaSeg under
`resources/regridded-inputs`; T1 uses linear interpolation and the discrete
inputs use nearest neighbor. FLAIR stays on its input grid and is registered
by the plus preprocessing stage. With `--no-regrid-inputs`, T1, mask, MaSeg,
and FLAIR are each reoriented by axis permutations/flips to closest RAS;
oblique direction cosines are ignored and no regridding interpolation is done.

Native Dura transfer defaults to `--dura-transfer-mode physical-affine`,
which maps source values to the full T1 grid using physical RAS coordinates.
Use `--dura-transfer-mode jist-unchanged` only to reproduce selected JIST
behavior: preserve source extents and values, assign the reference header, and
apply full-grid crop bounds with zero padding. The selected synthetic case
shows substantial downstream map and surface changes under this mode; it is
recorded in run signatures and manifests.

The required T1-grid mask replaces BET/SPECTRE skull stripping. MaSeg and the brain mask must match the T1 grid after RAS orientation normalization; a mismatch is rejected before preprocessing. The FLAIR may have a different physical grid; rigid ANTsPy registration resamples it onto T1. The atlas descriptor resolves its own compressed topology and shape maps. The top-level `plus-full-workflow.json` records ordered child manifests, input/settings signatures, stage state, and hashes. `--resume` verifies completed outputs and delegates repair of invalid checkpoints to each stage. The defaults follow the selected wrapper route; the ANTsPy registration backend is a documented substitute for NiftyReg/Aladin, not source-equivalent registration.

`--atlasem-descriptor` is an optional, separate AtlasEM structure descriptor.
When supplied, the runner executes the selected pre-LesionTOADS AtlasEM and
Dura handoff after T1/FLAIR registration, N4 correction, and user-mask
application. It passes Dura-stripped T1 and Dura-masked registered FLAIR into
LesionTOADS cropping. Every referenced AtlasEM topology/shape image must be a
compressed `.nii.gz` file beside the descriptor using its relative paths.
Omitting this option retains the prior LesionTOADS-only chain and does not
claim that the selected MATLAB PLUS source sequence is complete. The actual
selected `ibsr8obj` descriptor and payload are not bundled yet; the generated
AtlasEM descriptors in tests validate mechanics only and are not a substitute
for those selected resources.

Images are compressed NIfTI in RAS. Outputs include preprocessed T1/FLAIR and channel-last modalities; with `--atlasem-descriptor`, AtlasEM memberships, Dura products, and Dura-conditioned T1/FLAIR; atlas-aligned/recropped state; registered LesionTOADS products; Dura channels; and native-T1 transfer products. The auxiliary native CSF channel is retained for the plus-file API; selected GM, WM-fill, and WM-mask inputs use nearest-neighbor physical-affine transfer. This is a geometry-aware Python substitute for the selected JIST orientation-transfer boundary; arbitrary-grid numerical equivalence has not been established. MaSeg creates the selected cortical masks, then the plus continuation writes inner/central/outer level sets, RAS-mm VTK surfaces, volume thickness, and JSON metrics. If a synthetic case has no valid volume-thickness crossings, unavailable summary statistics are JSON `null` with `available: false`; this is distinct from a measured zero thickness.

Options include `--crop-padding`, `--init-align-iterations`, `--pyramid-levels`, `--lesion-iterations`, `--first-thinning-limit`, `--first-thinning-volume-ratio`, `--max-distance` (default `1e-4`), `--loop-convergence-threshold` (defaults to `--max-distance`), `--n-iter-max` (default `99`), `--loop-thinning-limit`, and `--max-passes` (default `99`). `plus-full-run` runs until terminal convergence or the configured iteration ceiling; use an explicit smaller `--max-passes` for a restartable bounded segment. `--finite-zero-variance` opts out of the source-faithful NaN behavior for zero-variance terms. `--dura-transfer-mode` defaults to `physical-affine`; select `jist-unchanged` to preserve Dura source extents/values, assign the reference header, and crop by full-grid MaSeg bounds with padding.

For a compact committed synthetic case, use `tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt` and the matching NIfTI topology/shape maps. `tests/test_plus_full_workflow.py` runs a 25³ folded-cortex synthetic T1/FLAIR/MaSeg/mask through the complete CLI composition with low-cost settings, including resume and corruption repair. This validates synthetic workflow composition, not composite Java equivalence or real-MRI behavior.
