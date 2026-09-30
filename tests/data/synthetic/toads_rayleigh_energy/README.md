# TOADS RAYLEIGH intensity-energy fixture

This asymmetric 11×9×7, three-channel synthetic fixture contains `images`,
`field`, and `psf_bias` NIfTI volumes with the same true-RAS affine and
anisotropic millimeter spacing. Its intensity field, bias and image values
produce both positive and nonpositive corrected intensities. The Java
reference covers all four combinations of field on/off and fixed bias
on/off. One class/channel scale is zero and another is negative to validate
IEEE results through the already separately-tested Rayleigh kernel.

Regenerate the fixture with:

```sh
python scripts/generate_toads_rayleigh_energy_fixture.py
```

Regenerate `java_reference.json` with:

```sh
scripts/generate_java_toads_rayleigh_energy_reference.sh
```

The Java harness composes the selected RAYLEIGH branch of
`MultipleToadSegmentation.computeExactMemberships()` and verifies the source
SHA-256 `63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Only interior class-energy voxels are evaluated; the output fixture uses
zero-valued six-face borders to represent untouched faces. Centroid/variance
estimation, PSF/field estimation and all membership updates are excluded.
