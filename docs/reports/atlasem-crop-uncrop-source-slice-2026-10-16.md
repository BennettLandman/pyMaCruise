# AtlasEM crop and uncrop source slice — 2026-10-16

## Ported behavior

`atlas_em_preprocessing.py` now implements three source-shaped operations from `AtlasEMPreprocess`:

- `find_atlas_em_crop_bounds` uses the selected robust-span crop predicate across all image channels and returns inclusive lower/upper voxel indices. It also supports the source's already-normalized predicate.
- `crop_atlas_em_images` copies the inclusive crop, adds the configured zero border, and adjusts the six-parameter rigid transform using the source rotation matrix and signed voxel resolutions.
- `uncrop_atlas_em_images` restores the cropped image into the original grid using the source border expansion/clipping and zero-fills outside it. It applies the source inverse translation update and returns that transform alongside the image.
- `recrop_atlas_em_after_registration` follows the selected wrapper order after registration: install the optimizer's current transform, uncrop, recompute foreground bounds, then crop and update the transform again.

The functions use `(x,y,z,channel)` arrays. This preserves the source's voxel/channel values while matching the Python project's existing NIfTI array convention.

## Source comparison

The synthetic Java harness pins `AtlasEMPreprocess.java` SHA-256 `acec23f916b3134e53381e83ff3de7ca569f785e6b866bfe135fb409c2926d2b` and the copied selected `RotationMatrix.java` SHA-256 `144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0`. On a two-channel 13×11×9 case with disconnected foreground contributions, anisotropic signed spacing, a nonzero rigid transform, and a two-voxel border, Python matches the Java excerpt's inclusive bounds, cropped dimensions, complete crop and uncrop arrays, and float32 transform bits.

Focused validation:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_cropping.py
```

Result on the tested working tree: **1 passed in 0.97 s**. The test also applies a synthetic registration translation change and compares the subsequent crop bounds, image values, and transform bits with a Java source excerpt. This validates the recropping wrapper sequence, not the optimizer that produces that transform. It does not validate full AtlasEM setup or the unavailable selected atlas configuration.

## Remaining boundary

The selected atlas descriptor and payload are still missing. AtlasEM rigid registration, any source-used deformable registration, and full workflow assembly remain open under P44b.
