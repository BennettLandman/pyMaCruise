# Larger asymmetric plus workflow fixture

This committed synthetic fixture scales the folded 25³ phantom in
`plus_pre_tgdm_composition` to an asymmetric **36×45×34** RAS grid while
preserving its physical field of view. The binary tissue and brain-mask
supports use nearest-neighbor scaling. T1/FLAIR intensity values and the
MaSeg labels are deterministic; FLAIR has a 0.25 mm RAS-x origin offset so
the workflow exercises registration. `provenance.json` records hashes, grids,
data types, and affines for all four `.nii.gz` inputs.

Regenerate the committed images with:

```sh
./.venv/bin/python tests/data/synthetic/plus_larger_grid/generate_fixture.py
```

The fixture has 3.5 times as many voxels as the prior 25³ end-to-end plus
case. It is still a compact synthetic field of view, not a brain-resolution
phantom, Java/source equivalence result, or real-MRI validation.
