# TOADS registration pyramid fixture

`java_reference.json` contains the output of the test-only Java harness
`tests/java_reference/ToadsPyramidSubsampleHarness.java`, which transcribes the
supplied VandyCRUISE `ImageFunctions.subsample(float[][][][])` block-mean loops.
The 9×8×7×3 deterministic input is sampled at factors 2 and 4. Both factors
discard incomplete blocks; factor 4 also proves that a trailing short z axis
uses floor dimensions rather than padding or clipping. Values are serialized
in x/y/z/class order.

This is an algorithmic numeric fixture, not a medical image. No imaging data
are included.
