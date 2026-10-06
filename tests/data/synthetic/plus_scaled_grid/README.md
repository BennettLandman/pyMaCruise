# Scaled synthetic plus-workflow fixture

This synthetic fixture scales the folded tissue, white-matter, and brain-mask
supports from `plus_pre_tgdm_composition` to an asymmetric **64×80×64** RAS
grid, preserving the source fixture's physical field of view. It contains a
float32 T1, a float32 FLAIR with a 0.25 mm RAS-x origin shift, a uint8 brain
mask, and uint16 MaSeg labels. T1/FLAIR values are generated deterministically
from the tissue masks. `provenance.json` records image hashes, shapes, data
types, axis codes, and affines. No human imaging data is used.

Regenerate the committed `.nii.gz` files with:

```sh
./.venv/bin/python tests/data/synthetic/plus_scaled_grid/generate_fixture.py
```

The fixture has about **5.8 times** the voxels of
`plus_larger_grid` (64×80×64 versus 36×45×34). It is a resource-scaling
phantom, not brain-resolution validation or source-equivalence evidence.
