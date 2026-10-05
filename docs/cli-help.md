# PyMaCRUISE command-line input/output reference

This page is the expanded `--help` contract for the supported inputs and
outputs. The command line exposes preparation, individual scientific stages,
restartable base and plus workflows, workflow configuration/preflight, and
standalone algorithm commands. Synthetic full-workflow validation is available;
source-equivalence and real-MRI acceptance remain in progress.

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
`plus_options` object using the field names in `PlusPipelineOptions`. Either
profile may provide `run_options`: base fields are names in `BaseRunOptions`;
plus fields are nested `registration`, `n4`, and `plus` objects matching the
corresponding runner configuration fields. Unknown fields are rejected. These
let synthetic runs lower iteration counts without changing algorithm defaults
in the implementation. The plus profile uses
the bundled LesionTOADS atlas descriptor. The base profile uses the bundled
default TOADS descriptor. Custom atlas paths are not accepted by this
preflight interface.

## `pymacruise workflow-run`

```text
pymacruise workflow-run --config workflow.json [--resume]
```

Executes the restartable workflow selected by the same version-1 JSON
configuration. Base runs use the bundled TOADS descriptor; plus runs use the
bundled LesionTOADS descriptor and pass `plus_options` through to the existing
MaCRUISEplus stages. The output directory must be new or contain a compatible
restartable run. If it contains the `run.json` produced by
`workflow-preflight`, that dry-run manifest is preserved as `preflight.json`
before execution begins. `--resume` verifies existing checkpoints and resumes
the configured workflow. Non-RAS inputs are normalized by the runners with
visible warnings; persisted image outputs are RAS `.nii.gz`.

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
three-vector NIfTI image. Defaults are lambda 1.0, 200 source multigrid
V-cycles, and absolute maximum-residual threshold 1e-6.

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
  [--registration-iterations 50] [--post-refresh-registration-iterations 1] \
  [--registration-max-iterations 5] \
  [--padding 10] [--maximum-iterations 50] \
  [--convergence-threshold 1e-4] [--relation-limit 3] \
  [--maximum-distance 1e-4] [--category-exponents CSF GM WM] \
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

### Output

The command writes one uint8 classification map on the original T1 grid with
RAS orientation and affine. Labels are restored from the registered crop using
the selected TOADS padding and x-fast uncrop semantics. The optional Dura-removal files are inputs for a later stage; this command
does not execute dura removal, run the second MaCRUISE TOADS pass, create final
CRUISE output products, or execute the complete base/plus workflow. The broader workflow CLI
and synthetic full-workflow validation remain open.

For the selected synthetic handoff, use the three exported maps required by
the current Dura CLI input contract and retain the original T1 as its reference:

```text
pymacruise toads-loop --t1-n4 T1_N4.nii.gz --brain-mask BrainMask.nii.gz \
  --output Classification.nii.gz --dura-output-dir toads_dura_inputs
pymacruise dura --t1 T1.nii.gz \
  --wm-fill toads_dura_inputs/T1_N4_wmfill.nii.gz \
  --gm toads_dura_inputs/T1_N4_gm.nii.gz \
  --wm-mask toads_dura_inputs/T1_N4_wmmask.nii.gz \
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

## `pymacruise toads-cruise-core`

```text
pymacruise toads-cruise-core --gm Toads2_gm.nii.gz \
  --wm-fill Toads2_wmfill.nii.gz --wm-mask Toads2_wmmask.nii.gz \
  [--segmentation MaSeg.nii.gz --remove-cerebellum] \
  --output-dir cruise_core [--edge-width 5] \
  [--ace-diffusion-iterations 25] [--gvf-lambda 1] \
  [--gvf-iterations 200] [--tgdm-initial-iterations 6] \
  [--tgdm-inner-iterations 4] [--tgdm-central-iterations 7] \
  [--provenance run.json]
```

This continues from second-pass maps through cleanup/crop/topology,
base-parameter ACE, generalized GVF, nested TGDM, and native-grid inverse
crop. It writes cropped tissue maps, topology-corrected WM-fill, ACE outputs,
the component-last GVF, and native-RAS inner/central/outer level sets. Base
ACE does not merge updated multi-atlas CSF. The selected MATLAB wrapper passes
the central-surface iteration option to the outer TGDM stage as well, so the
CLI preserves that mapping. Surfaces and thickness outputs are not generated.


## `pymacruise plus-run`

```text
pymacruise plus-run --t1-reference T1.nii.gz --segmentation MaSeg.nii.gz \
  --dura-csf CSF.nii.gz --dura-gm GM.nii.gz --dura-wm-fill WMFill.nii.gz \
  --dura-wm-mask WMMask.nii.gz --multi-atlas-gm AtlasGM.nii.gz \
  --multi-atlas-wm AtlasWM.nii.gz --multi-atlas-mask AtlasMask.nii.gz \
  --brain-mask BrainMask.nii.gz --output-dir plus-run [--resume] \
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
metrics, and `run.json`.

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
run stops, including counts of nonfinite energy values and memberships plus
zero-standard-deviation channels. `--finite-zero-variance` enables an explicit stabilization: matching
centroids with zero deviation have distance zero and nonmatching ones infinite
distance. The summary records zero-deviation channels and nonfinite membership
count so this numerical difference remains visible. All outputs use the
registered RAS grid. This step does not execute later LesionTOADS updates.

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
run stops, including counts of nonfinite energy values and memberships plus
zero-standard-deviation channels. `--finite-zero-variance` enables an explicit stabilization: matching
centroids with zero deviation have distance zero and nonmatching ones infinite
distance. The summary records zero-deviation channels and nonfinite membership
count so this numerical difference remains visible. All outputs use the
registered RAS grid. This step does not execute later LesionTOADS updates.

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
  [--first-thinning-volume-ratio 0.5] [--max-distance 0] \
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
  [--finite-zero-variance] [--resume]
```

Continue from a registered branch checkpoint or a previous registered-loop
checkpoint. `--max-passes` bounds work for this invocation; reaching that
boundary is recorded as a checkpoint, not as source convergence. The manifest
stores the iteration number, RAS outputs, atlas/parent hashes, and selected
settings so a later invocation can continue without replaying the completed
segment. Source-faithful zero-variance behavior remains the default. This is
synthetic selected-loop integration; natural convergence, complete Java
equivalence, and real-image validation remain open.

## `pymacruise plus-lesion-loop-run`

```text
pymacruise plus-lesion-loop-run \
  --registered-recrop-manifest recrop/registered-recrop.json \
  --atlas-descriptor LesionTOADS-atlas.txt --output-dir lesion-loop \
  [--iterations 99] [--first-thinning-limit 3] \
  [--first-thinning-volume-ratio 0.5] [--max-distance 0] \
  [--loop-convergence-threshold N] [--n-iter-max 99] \
  [--thinning-limit 3] [--max-passes 1] \
  [--finite-zero-variance] [--resume]
```

Compose the persisted P41i initial membership/thinning, P41k registered
branch, and one bounded P41m registered-loop segment from a completed P41h
recrop manifest. A top-level JSON manifest links each child manifest by hash
and records stage state before execution. `--resume` asks each child checkpoint
to verify and reuse its outputs; corrupt loop outputs are recomputed. Source-
faithful zero-variance behavior remains the default. `--finite-zero-variance`
selects the explicit finite synthetic path. An optional
`--loop-convergence-threshold` overrides only P41m's convergence threshold;
`--max-distance` continues to control the P41k registered branch. A segment
that reaches convergence or the iteration cap writes terminal classification
and lesion-subtype NIfTIs; otherwise it writes a continuing checkpoint. This
does not establish complete MaCRUISEplus, Java composite equivalence, or
real-image validation.

## `pymacruise plus-full-run`

Run the restartable raw-input plus path from T1, FLAIR, a user-provided T1-grid brain mask, MaSeg labels, and a descriptor-owned atlas through preprocessing, registered LesionTOADS, Dura-channel transfer, and terminal CRUISE surfaces and thickness.

```bash
pymacruise plus-full-run --t1 T1.nii.gz --flair FLAIR.nii.gz \
  --brain-mask brain-mask.nii.gz --maseg MaSeg.nii.gz \
  --atlas-descriptor atlas.txt \
  --output-dir plus-output \
  --registration-iterations 1000 500 250 0 \
  --shrink-factors 6 4 2 1 --smoothing-sigmas 3 2 1 0 \
  --n4-iterations 50 50 50 50 --resume
```

The required T1-grid mask replaces BET/SPECTRE skull stripping. MaSeg and the brain mask must match the T1 grid after RAS orientation normalization; a mismatch is rejected before preprocessing. The FLAIR may have a different physical grid; rigid ANTsPy registration resamples it onto T1. The atlas descriptor resolves its own compressed topology and shape maps. The top-level `plus-full-workflow.json` records ordered child manifests, input/settings signatures, stage state, and hashes. `--resume` verifies completed outputs and delegates repair of invalid checkpoints to each stage. The defaults follow the selected wrapper route; the ANTsPy registration backend is a documented substitute for NiftyReg/Aladin, not source-equivalent registration.

Images are compressed NIfTI in RAS. Outputs include preprocessed T1/FLAIR and channel-last T1/FLAIR modalities, atlas-aligned/recropped state, registered LesionTOADS products, Dura channels, and native-T1 transfer products. The auxiliary native CSF channel is retained for the plus-file API; selected GM, WM-fill, and WM-mask inputs use nearest-neighbor physical-affine transfer. This is a geometry-aware Python substitute for the selected JIST orientation-transfer boundary; arbitrary-grid numerical equivalence has not been established. MaSeg creates the selected cortical masks, then the plus continuation writes inner/central/outer level sets, RAS-mm VTK surfaces, volume thickness, and JSON metrics. If a synthetic case has no valid volume-thickness crossings, unavailable summary statistics are JSON `null` with `available: false`; this is distinct from a measured zero thickness.

Options include `--crop-padding`, `--init-align-iterations`, `--pyramid-levels`, `--lesion-iterations`, `--first-thinning-limit`, `--first-thinning-volume-ratio`, `--max-distance`, `--loop-convergence-threshold`, `--n-iter-max`, `--loop-thinning-limit`, and `--max-passes`. `--finite-zero-variance` opts out of the source-faithful NaN behavior for zero-variance terms.

For a compact committed synthetic case, use `tests/data/synthetic/plus_lesion_atlas_init/mini-lesion-toads-atlas.txt` and the matching NIfTI topology/shape maps. `tests/test_plus_full_workflow.py` runs a 25³ folded-cortex synthetic T1/FLAIR/MaSeg/mask through the complete CLI composition with low-cost settings, including resume and corruption repair. This validates synthetic workflow composition, not composite Java equivalence or real-MRI behavior.
