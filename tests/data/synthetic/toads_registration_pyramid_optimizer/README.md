# Synthetic pyramid registration images

This committed synthetic fixture exercises the coarse-to-fine rigid atlas
registration wrapper. It contains no clinical or subject data.

| File | Shape (x, y, z, class) | RAS voxel spacing | Contents |
|---|---:|---:|---|
| `memberships.nii.gz` | 17×15×13×3 | 0.8×1.2×1.5 mm | Float32 membership values from NumPy `default_rng(2)`, uniform on [0.05, 0.9). |
| `atlas_shapes.nii.gz` | 19×17×15×3 | 0.9×1.1×1.3 mm | Float32 synthetic class-prior values from the same generator immediately after memberships. |

Both NIfTI files use positive RAS affine axes and preserve class as the final
fourth dimension. Tests verify the RAS axis codes, shape, and every voxel
against the deterministic generator before running scales 4 and 2. The
different odd input and atlas dimensions exercise independent block truncation
and full-resolution center handling.
