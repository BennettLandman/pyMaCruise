# TOADS intensity variance factor reference fixture

This tiny metadata-only descriptor separates the `Intensity Variance Factor`
table from the distinct `Intensity Atlas` intensity means. Its three
structures make class order visible. It supplies MPRAGE and FLAIR variance
rows, while T1/T2/PD and an unknown modality exercise Java's initialized
all-ones defaults. The image values are intentionally irrelevant: this tests
descriptor semantics, not segmentation or numerical workflow equivalence.

The expected output is emitted by the standalone Java harness that mirrors
`MultipleToadDeformableAtlas` in the supplied TOADS source:

```sh
javac -d /tmp/pymacruise-variance-reference \
  tests/java_reference/ToadsIntensityVarianceFactorReferenceHarness.java
java -cp /tmp/pymacruise-variance-reference \
  ToadsIntensityVarianceFactorReferenceHarness \
  tests/data/synthetic/toads_intensity_variance_factors/mini-atlas.txt \
  MPRAGE FLAIR T1 T1_SPGR T2 T1_RAW MPRAGE_RAW PD PDFSE DIR UNKNOWN_MODALITY \
  > tests/data/synthetic/toads_intensity_variance_factors/java_reference.json
```

The selected source references are `getIntensityVariancePriors` at lines
270-279, descriptor initialization of `intensityVariance` to `1.0f` at lines
939-942, variance-section parsing at lines 1114-1124, and modality aliases in
`modalityId` at lines 1262-1271 of
`toads-cruise/src/edu/jhmi/rad/medic/methods/MultipleToadDeformableAtlas.java`.
