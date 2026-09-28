# Optimized TOADS class-variance phantom

This fixture reuses the 7×6×5×3 RAS, 1 mm isotropic synthetic image and
multiplicative-field volumes in `../toads_optimized_likelihood/`, and the
first three membership channels in `../toads_direct_membership/`. It adds an
int16 categorical segmentation volume, where each voxel is assigned by
`(x + 2*y + z) % 3`. Three initial centroids are recorded by the Java
formula harness.

`java_reference.json` stores the assigned-class standard deviations for each
modality and class, with and without the field, using the optimized tissue
branch of `MultipleToadSegmentation.computeExactVariances`. The selected
source file SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
Regenerate the categorical image with
`scripts/generate_toads_optimized_class_variance_fixture.py` and the Java
reference with `scripts/generate_java_toads_optimized_class_variance_reference.sh`.
