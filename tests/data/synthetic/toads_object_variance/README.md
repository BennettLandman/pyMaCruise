# Base TOADS ordinary-object variance phantom

This asymmetric `4 x 3 x 3` synthetic phantom stores images, optional
multiplicative intensity field, four membership channels, and a categorical
segmentation as `.nii.gz`. Every NIfTI uses a 1 mm identity affine with RAS
axis codes. Three classes occur in the segmentation. Class 2 is assigned only
at two boundary corners, and class 3 has no assigned voxels. Class 1 has zero
memberships at all its assigned voxels, exercising the source `1e-30f` weight
floor. The two ordinary parity classes also reach the volume faces.

`ToadsObjectVarianceReferenceHarness.java` isolates the ordinary-object
formula from the selected `MultipleToadSegmentation.computeExactVariances`
source. It uses the source's x/y/z traversal over the entire image, optional
`field * image` correction, left-associated float32 expression
`(mem*mem + 1e-30f) * diff * diff`, float32 numerator and denominator sums,
and square root of their ratio. The fixture classes are ordinary OBJECTs;
the separate lesion sub-branch is intentionally outside this case. The source
file SHA-256 is
`63e135f36dd9117a4f6dda2ec78990230e5e44ad91234ce3745d20e5976cd884`.
The Java harness reports raw float bits in `java_reference.json`.

Regenerate input volumes with:

```sh
.venv/bin/python tests/data/synthetic/toads_object_variance/generate_fixture.py
```

Regenerate the formula reference with:

```sh
mkdir -p build/java-reference/toads-object-variance
javac -d build/java-reference/toads-object-variance tests/java_reference/ToadsObjectVarianceReferenceHarness.java
java -cp build/java-reference/toads-object-variance ToadsObjectVarianceReferenceHarness > tests/data/synthetic/toads_object_variance/java_reference.json
```
