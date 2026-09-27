# Synthetic Fixture Corpus

Every fixture here is project-generated synthetic data; no human image data is included. Image volumes are stored as compressed NIfTI (`.nii.gz`). The mesh fixture stores only scientific geometry (RAS+ vertices in millimetres and triangle indices) in compressed NumPy format; it has no visualization dependency. `phantom.json` records the seed, recipe, labels, affine, spacing, coordinate frame, and fixture purpose. `SHA256SUMS` covers every generated file except itself.

## Files

| File | Contents / use | Geometry / values |
|---|---|---|
| `T1.nii.gz` | Main T1-like volume with seeded noise and smooth multiplicative bias. | 32×28×24; anisotropic 1.0×1.25×1.5 mm; RAS, 7° in-plane obliquity. |
| `FLAIR.nii.gz` | Multimodal companion with a bright synthetic lesion. | Same affine/grid as `T1.nii.gz`. |
| `MaSeg.nii.gz` | Synthetic segmentation labels. | Background 0, WM 4, GM 100; lesion remains a separate mask. |
| `WM_mask.nii.gz`, `GM_mask.nii.gz`, `brain_mask.nii.gz` | Exact tissue masks. | `uint8`, 0/1 values, same affine as T1. |
| `lesion_mask.nii.gz`, `cavity_mask.nii.gz` | Asymmetric lesion and interior background cavity. | Binary `uint8` masks, same affine as T1. |
| `bias_field.nii.gz` | Known smooth multiplicative bias. | `float32`, positive-valued, same affine as T1. |
| `non_ras_orientation.nii.gz` | Tiny LPS-oriented image with distinct corner markers. | 5×7×9; intentionally non-RAS for warning/reorientation checks. |
| `boundary_labels.nii.gz` | Adjacent regions with a one-voxel boundary and multiple labels. | 8×8×8 `uint16`. |
| `disconnected_labels.nii.gz` | Two separated components with the same GM label and a distinct WM component. | 16×16×16 `uint16`. |
| `brain_surface_ras.npz` | Closed ellipsoid triangle mesh. | Vertices are RAS+ millimetres; faces are triangle indices. |

Regenerate with `PYTHONPATH=src .venv/bin/python -m pymacruise.synthetic --output tests/data/synthetic --seed 20260927`.
