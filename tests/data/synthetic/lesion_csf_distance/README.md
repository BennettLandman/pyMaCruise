# Selected LesionTOADS CSF distance fixture

This synthetic comparison covers the exact atlas-name selection in
`LesionToadSegmentation.computeCSFDistance()` and its call to the selected
byte `fastDistance(boolean[][][], byte, int, int, int)` overload. The selected
wrapper comments out the helper call, initializes `CSFDist` to zero, and does
not export the result. This fixture validates a dormant support method, not
active LesionTOADS behavior or a segmentation.

`seeded_atlas_indices.nii.gz` is an asymmetric `(9, 10, 11)` integer-label
volume. It contains every exact seed name, a duplicate `Caudate` table entry,
border and interior seeds, and case/punctuation/whitespace controls.
`near_miss_atlas_indices.nii.gz` contains only negative-control spellings and
therefore yields no seeds. Both volumes have positive RAS axis directions,
anisotropic spacing, and nonzero origin. The Java arrays and flattened JSON
use `(x, y, z)` order. Cases exercise cap 1's source extra propagation and the
field-initialization default cap of 5. The distance kernel measures six-
connected voxel steps and ignores physical spacing.

Exact selected names: `Caudate`, `Thalamus`, `Putamen`, `Subcortical-GM`,
`SubcorticalGM`, and `BasalGanglia`. Matching is case-, punctuation-, and
whitespace-sensitive; duplicated matching names seed both corresponding
atlas indices.

Provenance:

- Source: selected VandyCRUISE `LesionToadSegmentation.java`
- Source SHA-256: `545a61027d22efcc5d07d6bea960039c938808867707cbbd61811cad46cc3a99`
- Wrapper SHA-256: `556940249cb3a8191f53f23a76308f11823789f7b41f7a83b071d784ae4a20cb`
- Excerpt: `computeCSFDistance()` plus selected byte `fastDistance(...)`
- Standalone harness: `tests/java_reference/LesionToadsCsfDistanceReferenceHarness.java`
- Harness SHA-256: `94384a693850c706c10744e8dcf4816c28fe4f7c37f81420117e9ded874f759b`
- Generated Java reference SHA-256: `160836366079cf9bfc57a953fa582439161ecf833bd765cf24d6df6599274269`

Regenerate the reference and image phantoms from the repository root:

```sh
mkdir -p /tmp/pymacruise-csf-distance-reference
javac -d /tmp/pymacruise-csf-distance-reference \
  tests/java_reference/LesionToadsCsfDistanceReferenceHarness.java
java -cp /tmp/pymacruise-csf-distance-reference \
  LesionToadsCsfDistanceReferenceHarness \
  > tests/data/synthetic/lesion_csf_distance/java_reference.json
PYTHONPATH=src .venv/bin/python \
  tests/data/synthetic/lesion_csf_distance/generate_fixtures.py
```

The generator writes only the two NIfTI files in this directory. The Java
command overwrites only the named JSON reference.
