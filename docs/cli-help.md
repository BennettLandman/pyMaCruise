# PyMaCRUISE command-line input/output reference

This page is the expanded `--help` contract for the supported inputs and
outputs. The command line currently exposes preparation and individual
scientific stages; the full base and plus workflow commands remain under
development.

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

There is no full TOADS command and no full LesionTOADS command. The repository
contains selected supporting stages and synthetic TOADS-like or
LesionTOADS-related components, but these do not constitute either complete
algorithm. The TOADS/SPECTRE/BET port-vs-dependency decision remains the
tracked blocker for the complete TOADS workflow. These standalone interfaces
are not substitutes for either full pipeline.
