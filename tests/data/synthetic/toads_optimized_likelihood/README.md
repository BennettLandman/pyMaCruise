# Optimized TOADS tissue likelihood phantom

This fixture stores three 7×6×5×3 float32 RAS `.nii.gz` volumes: normalized
synthetic input channels, a multiplicative intensity field, and a PSF bias
field. All are generated deterministically by
`scripts/generate_toads_likelihood_fixture.py`, use a 1 mm RAS affine, and
are numeric phantoms rather than human images.

`java_reference.json` contains per-voxel CSF/GM/WM data-term energies from a
Java formula excerpt of selected `MultipleToadSegmentation.computeExactMemberships`
branches, calling the actual copied source `PowerTable` implementation. It
covers all 630 voxel/class energies with the optional field both disabled and
enabled; PSF bias and class/modality/variance weights are applied in both
cases. The fixture does not include smoothing, atlas priors, normalization,
membership updates, or topology evolution.

Regenerate the numeric inputs with `scripts/generate_toads_likelihood_fixture.py`
and the Java reference with
`scripts/generate_java_toads_optimized_likelihood_reference.sh`.
