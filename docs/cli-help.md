# PyMaCRUISE command-line input/output reference

This page is the expanded `--help` contract for the supported inputs and
outputs. The command line exposes preparation, individual scientific stages, and a
standalone registered TOADS loop. The complete two-pass base and plus workflow
commands remain under development.

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

The plus profile additionally requires `flair_path` and may provide a
`plus_options` object using the field names in `PlusPipelineOptions`. It uses
the bundled LesionTOADS atlas descriptor. The base profile uses the bundled
default TOADS descriptor. Custom atlas paths are not accepted by this
preflight interface.

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
  [--max-sweeps COUNT] [--no-deal-with-central] [--no-source-faithful] \
  [--provenance run.json]
```

The three inputs must be 3D `.nii.gz` images on the same physical grid after
RAS reorientation. `--summary` records the selected mode and per-rule
candidate, changed, unresolved, sweep, and threshold counts. The historical
`source_faithful` behavior is the default while owner sign-off is pending; it
preserves the V3 MATLAB `maxind` row-index result. Use `--no-source-faithful`
to opt into assigning the actual modal label. The label-refinement output is
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

## Other currently exposed commands

- `pymacruise --help` lists the current command groups and their required
  options.
- `pymacruise --version` reports the installed package version.
- `pymacruise cortical-separation --maseg MaSeg.nii.gz --output-dir DIR`
  writes four RAS-oriented binary `.nii.gz` masks derived from MaSeg labels.
- `pymacruise calculate --input-a A.nii.gz --input-b B.nii.gz --operation OP
  --output RESULT.nii.gz` performs the selected generic voxelwise operator;
  both images must share one physical grid. Supported operations are `ADD`,
  `SUBTRACT`, `MULTIPLY`, `DIVIDE`, `MIN`, and `MAX` (case-insensitive).
  Inputs are axis-reoriented to RAS with a warning when needed, and output is
  `.nii.gz`. Optional `--provenance RUN.json` records the effective operation,
  paths, software version, and RAS convention.

The full-pipeline configuration file, restart behavior, and final output
manifest will be specified when P34 and P42a are implemented. See
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
  [--rescale-intensities] [--provenance RUN.json]
```

Reads a finite, non-negative 3D intensity image and writes its N4-corrected
image as `.nii.gz` in RAS orientation. The defaults match the selected
wrapper's unadorned ANTs N4 command: shrink factor 4, four 50-iteration
levels, zero convergence tolerance, a 1×1×1 spline mesh, and intensity
rescaling off. ANTsPy is optional; install PyMaCRUISE with its `ants` extra
to run this command. The N4 estimation itself is provided by ANTsPy/ANTs.

### `pymacruise ace`

```text
pymacruise ace --gray GM.nii.gz --white WM.nii.gz \
  --enhanced GM-enhanced.nii.gz --skeleton ACE-skeleton.nii.gz \
  --thinned-skeleton ACE-skeleton-thinned.nii.gz [--atlas-csf CSF.nii.gz]
```

Inputs are aligned 3D gray/white matter membership or intensity images; the
optional CSF atlas image must share their physical grid. Writes enhanced gray
matter and both ACE skeleton products. `--diffusion-iterations` defaults to
25. This exposes the current ACE stage, not all multi-atlas MaACE options or
the complete MaCRUISE workflow.

### `pymacruise gvf`

```text
pymacruise gvf --input IMAGE.nii.gz --output GVF.nii.gz \
  [--lambda VALUE] [--max-iterations N] [--solver-tolerance VALUE]
```

Reads a finite 3D scalar image and writes a RAS-oriented, component-last
three-vector NIfTI image. Defaults are lambda 1.0, 200 iterations, and solver
tolerance 1e-6.

### `pymacruise tgdm`

```text
pymacruise tgdm --initial INITIAL.nii.gz --white WM.nii.gz \
  --gray GM.nii.gz --ace-gray GM-ACE.nii.gz --output-dir DIR \
  [--brain-mask MASK.nii.gz] [--ace-skeleton SKELETON.nii.gz] [--gvf GVF.nii.gz]
```

The four required finite 3D scalar images and each supplied optional image
must share one physical grid. Optional GVF is a component-last vector image.
Writes `initial_level_set.nii.gz`, `tissue_sign.nii.gz`,
`initial_inner.nii.gz`, `inner.nii.gz`, `central.nii.gz`, and `outer.nii.gz`,
all on the RAS grid. This is the nested TGDM stage, not full TOADS or MaCRUISE.

### `pymacruise dura`

```text
pymacruise dura --t1 T1.nii.gz --white WM.nii.gz --gray GM.nii.gz \
  --white-mask WM-mask.nii.gz --output-dir DIR [--max-iterations N]
```

The T1, white/gray membership, and original white-matter mask must be finite
3D images on one physical grid. `--max-iterations` defaults to 5. Writes
RAS-oriented `stripped_t1.nii.gz`, `corrected_white_matter.nii.gz`,
`corrected_gray_matter.nii.gz`, `brain_mask_level_set.nii.gz`,
`brain_mask.nii.gz`, `wm_topology.nii.gz`, and `gm_topology.nii.gz`.

### `pymacruise surface-extract`

```text
pymacruise surface-extract --level-set LEVELSET.nii.gz --output SURFACE.vtk [--level VALUE]
```

Reads a finite 3D scalar level set and writes a triangular RAS+ millimetre
surface. The isovalue defaults to zero. The VTK header declares coordinate
frame and units.

### `pymacruise thickness-levelset`

```text
pymacruise thickness-levelset --inner INNER.nii.gz --outer OUTER.nii.gz \
  --output THICKNESS.nii.gz [--iso-level VALUE]
```

Inputs must share one physical RAS grid. Writes the voxelwise absolute
level-set difference in millimetres and prints mean, sample standard
deviation, and crossing-voxel count for the selected iso-level (default 0).

### `pymacruise thickness-surface`

```text
pymacruise thickness-surface --inner INNER.vtk --outer OUTER.vtk \
  --inner-output INNER-MEASURED.vtk [--central CENTRAL.vtk \
  --central-output CENTRAL-MEASURED.vtk]
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
  [--registration-iterations 50] [--registration-max-iterations 5] \
  [--padding 10] [--maximum-iterations 50] \
  [--convergence-threshold 1e-4] [--relation-limit 3] \
  [--maximum-distance 1e-4] [--category-exponents 2 2 2] \
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
| Registration and loop options | No | Selected source controls with defaults shown above. `--category-exponents` sets CSF, GM, and WM exponents for the single T1 channel. |

### Output

The command writes one uint8 classification map on the original T1 grid with
RAS orientation and affine. Labels are restored from the registered crop using
the selected TOADS padding and x-fast uncrop semantics. The optional Dura-removal files are inputs for a later stage; this command
does not execute dura removal, run the second MaCRUISE TOADS pass, create final
CRUISE output products, or execute the complete base/plus workflow. The broader workflow CLI
and synthetic full-workflow validation remain open.
