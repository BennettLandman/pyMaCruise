# Brain-resolution synthetic base workflow inputs

These committed `(172,216,164)` `.nii.gz` inputs are deterministic atlas-derived
synthetic data for the restartable base workflow. They are not clinical data
and do not establish source equivalence or anatomical accuracy. All three
images are on the same RAS grid.

The fixture was generated with the repository's bundled VandyCRUISE atlas
maps. To regenerate a copy for comparison:

```bash
./.venv/bin/python tests/data/synthetic/toads_downsampled_vandy_atlas/generate_fixture.py \
  --target-shape 172 216 164 \
  --output-dir build/vandy-brain-resolution-memory
```

Copy only `t1.nii.gz`, `brain-mask.nii.gz`, and `maseg.nii.gz` into this
directory when updating the committed inputs; the generated atlas maps remain
in the repository's existing atlas assets. SHA-256 values:

| Image | SHA-256 |
|---|---|
| `t1.nii.gz` | `ebf945d06cc2519114763c949af0fc1ec88d938677748a41126852ab6782f158` |
| `brain-mask.nii.gz` | `2d617064fbea8e8a5d807165d2567d87a9b40f0d4487078810b25a755ad87a00` |
| `maseg.nii.gz` | `06cb8f49e12d7966cff025601695b41de9dae5526358ab23b2be7cb4a2631005` |
