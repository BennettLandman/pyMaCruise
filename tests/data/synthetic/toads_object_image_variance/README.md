# Base TOADS ordinary-object image variance

The fixture is an asymmetric `4 x 3 x 3` RAS grid with two modalities, three
ordinary OBJECT membership channels, and a multiplicative field. All image
volumes are committed as `.nii.gz` with 1 mm identity affine and true RAS axis
codes. Every class has nonzero membership on at least one face voxel. The
first modality is exactly equal to every class centroid, so its numerator is
zero and exercises the source `1e10f` fallback. The second modality exercises
nonzero residuals. With the optional field enabled, modality 0 remains
zero-residual while modality 1 is field-corrected.

`ToadsObjectImageVarianceReferenceHarness.java` is a pinned formula excerpt
from `MultipleToadSegmentation.computeExactImageVariance`, selected source
SHA-256 `63e135f36dd9117a4f6fdda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
It emits float32 bits for raw numerator and denominator, each modality's
unnormalized inverse variance, their float32 sum, pre-normalization standard
deviation, and normalized `varImg`, with and without the field. All object
classes are ordinary OBJECTs. LESION/OUTLIER mixture branches, approximate
variance, and full membership estimation are excluded.

Regenerate `.nii.gz` inputs with:

```sh
.venv/bin/python tests/data/synthetic/toads_object_image_variance/generate_fixture.py
```

Regenerate the pinned method-excerpt reference with:

```sh
mkdir -p build/java-reference/toads-object-image-variance
javac -d build/java-reference/toads-object-image-variance tests/java_reference/ToadsObjectImageVarianceReferenceHarness.java
java -cp build/java-reference/toads-object-image-variance ToadsObjectImageVarianceReferenceHarness > tests/data/synthetic/toads_object_image_variance/java_reference.json
```
