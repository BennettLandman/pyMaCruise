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
P34/P42 orchestration is implemented.

## Other currently exposed commands

- `pymacruise --help` lists the current command groups and their required
  options.
- `pymacruise --version` reports the installed package version.
- `pymacruise cortical-separation --maseg MaSeg.nii.gz --output-dir DIR`
  writes four RAS-oriented binary `.nii.gz` masks derived from MaSeg labels.
- `pymacruise calculate --input-a A.nii.gz --input-b B.nii.gz --operation OP
  --output RESULT.nii.gz` performs the selected generic voxelwise operator;
  both images must share one physical grid.

The full-pipeline configuration file, restart behavior, and final output
manifest will be specified when P34 and P42 are implemented. See
[`io-contract.md`](io-contract.md) for profile-level input/output semantics and
[`toads-port-plan.md`](toads-port-plan.md) for the TOADS blocker plan.
