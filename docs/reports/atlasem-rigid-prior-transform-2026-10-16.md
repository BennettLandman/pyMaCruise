# AtlasEM rigid-prior transform source slice — 2026-10-16

## Ported behavior

`transform_atlas_em_rigid_priors` in `src/pymacruise/atlas_em_preprocessing.py` now connects AtlasEM's rigid atlas-prior path to the existing source-checked TOADS sampler. It computes each class's positive-support image bounds, samples the atlas shape maps with the selected rigid image-to-template mapping and closest-edge trilinear interpolation, restores the descriptor's shape slope, applies the source default outside-bounds values, and normalizes class priors within the supplied mask.

This is a fixed-transform adapter. It does not parse an AtlasEM descriptor, estimate the rigid registration, perform deformable registration, or compose the missing `ibsr8obj` initialization.

## Source comparison

The Java reference harness captures the selected `DemonToadDeformableAtlas` rigid transformed-shape and bounding-box behavior, using the pinned `ParametricTransform` and `ImageFunctions.linearClosestInterpolation` implementations, followed by AtlasEM's masked `computeAtlas` prior normalization. Source identities are:

- `DemonToadDeformableAtlas.java`: `96cb5321f9cc95a84d9115b6ea638e9d84aded422171c8792906fc2130bd415a`
- `ParametricTransform.java`: `22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac`
- `ImageFunctions.java`: `be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892`

On a deterministic two-class synthetic volume with anisotropic spacing, a non-identity rigid transform, cropped positive support, and a nontrivial mask, the complete normalized `(x,y,z,class)` float32 output hash matches the Java excerpt exactly. The test also checks the active voxel count and pins the three source hashes.

Focused command:

```sh
PYTHONPATH=.:src PYTHONNOUSERSITE=1 OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 \
  .venv/bin/pytest -q tests/test_atlas_em_rigid_priors.py tests/test_atlas_em_preprocessing.py \
  tests/test_atlas_em_prior_normalization.py tests/test_toads_shape_volume.py
```

Result: **7 passed in 2.16 s**. This is a selected Java method-excerpt comparison on synthetic data, not a configured AtlasEM run or whole-workflow equivalence result.

## Remaining boundary

The selected rigid transform must still be produced by AtlasEM's registration setup, and the source descriptor/payload remains unavailable. The full selected AtlasEM initialization and its integration before pre-LesionTOADS Dura therefore remain open under P44b.
