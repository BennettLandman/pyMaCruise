# P18czr: post-registration TOADS template refresh phantom

This synthetic RAS fixture validates one bounded wrapper sequence from the selected VandyCRUISE AlgorithmMultipleToads.runAlgorithm() after its initial rigid registerShapesPyramid() call. It does **not** establish equivalence of the full TOADS workflow or validate the registration optimizer itself.

The sequence under test is:

1. MultipleToadPreprocess.uncropImages()
2. MultipleToadPreprocess.findCroppingBoundaries()
3. MultipleToadPreprocess.cropImages()
4. MultipleToadDeformableAtlas.setImageInfo(...)
5. MultipleToadDeformableAtlas.computeTransformedTemplate()
6. MultipleToadPreprocess.updateTransformedTemplate(...)

updateTransformedTemplate() only assigns the resulting labels and sets a state flag. The executable Python composition reuses the selected crop, uncrop, bounds, and topology-resampling kernels from pymacruise.toads_crop.

## Phantom

- native_t1_channels_ras.nii.gz: two deterministic synthetic intensity channels, shape 14×13×11×2, RAS affine with 0.8, 1.1, and 1.4 mm spacing.
- atlas_template_ras.nii.gz: asymmetric, folded multi-label topology, shape 17×16×15, RAS affine with 1.0, 1.2, and 1.5 mm spacing.
- Previous inclusive crop bounds: [1, 11, 1, 10, 1, 8], with two voxels of padding. The registered transform has nonzero rotations and translations; the original crop is deliberately offset from the native-grid center so uncrop and re-crop translation adjustments are exercised.
- Recomputed bounds extend beyond the native image. The selected crop code preserves these inclusive bounds and zero-fills out-of-grid image samples. This checks padded/border behavior instead of clamping the crop to the image.

## Reference

java_reference.json was emitted by tests/java_reference/ToadsPostRegistrationTemplateRefreshReferenceHarness.java. The harness implements method excerpts for the selected source sequence and records SHA-256 digests of every image/label output, exact transform float bits, the recomputed bounds, crop dimensions, and call order.

Pinned selected-source SHA-256 values:

- AlgorithmMultipleToads.java: b6ee1cd14dc0ad5eb7a3c0e9809ce19978c3381b233bfcc752ab070a50bbf01d
- MultipleToadPreprocess.java: 31296a4c154ac5366dcd78aa7e5da28bb1c56163f1e230cc451cb3e4feb071e7
- MultipleToadDeformableAtlas.java: e949505a3ccdd63745fea7ba8a1dd5a67d5c4c4eb9644553d1d63beb6217a9e2
- RotationMatrix.java: 144d16dc33eca12a4a3a117c362b7aa1e767fe6561950f64ef04a66297818af0
- ParametricTransform.java: 22bd0cf545bf1112f24096be54e7688638bbcdbd6224e69dc6c74ea240063fac
- ImageFunctions.java: be0eae41dd1487de200930e72ed22a3585402caf7a04484180bbb45d77468892
- Numerics.java: 13925b6bb7f5d58e15c3cd7d4c66bdba5b60de6762e688667babe5621844e6fb

Regenerate the NIfTI inputs with:

    PYTHONPATH=src .venv/bin/python tests/java_reference/generate_toads_post_registration_template_refresh_fixture.py

Regenerate the Java excerpt reference with:

    mkdir -p build/java-reference/p18czr
    javac -d build/java-reference/p18czr tests/java_reference/ToadsPostRegistrationTemplateRefreshReferenceHarness.java
    java -cp build/java-reference/p18czr ToadsPostRegistrationTemplateRefreshReferenceHarness > tests/data/synthetic/toads_post_registration_template_refresh/java_reference.json

Focused validation:

    PYTHONPATH=src .venv/bin/pytest -q tests/test_toads_post_registration_template_refresh.py

The test checks exact selected-Java bounds, image and label digests, transform float bits, dimensions, call order, RAS fixture affines, and input immutability. It does not validate the full atlas optimizer or subsequent TOADS iterations.
