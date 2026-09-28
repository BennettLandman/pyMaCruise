# TOADS atlas intensity-prior parsing reference

This checks parsing of the selected bundled base TOADS descriptor's intensity
table in its declared class order. It is a metadata test; it contains no
synthetic or subject image data. The expected T1 MPRAGE row was emitted by the
standalone Java harness that mirrors `MultipleToadDeformableAtlas.modalityId`
and its `Intensity Atlas` parser.

Regenerate the expected artifact with:

```sh
javac -d /tmp/pymacruise-atlas-priors tests/java_reference/ToadsIntensityPriorReferenceHarness.java
java -cp /tmp/pymacruise-atlas-priors ToadsIntensityPriorReferenceHarness \
  data/atlases/cruise-atlas-10obj-toads2008.txt T1_MPRAGE \
  > tests/data/synthetic/toads_atlas_priors/java_reference.json
```
